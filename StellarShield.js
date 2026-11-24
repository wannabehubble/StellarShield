// ============================================================================
// StellarShield.js  -  v1.3.2
//
// Restores galaxies and other non-stellar objects that a star-removal tool
// (StarXTerminator, StarNet, ...) removed by mistake.
//
// How it works
//   1. Everything the star remover took out = original - starless.
//   2. That residual is thresholded and split into connected blobs.
//   3. Each blob is cross-matched against Gaia (local XPSD database) and
//      measured (flux, size, ellipticity).
//   4. Blobs judged non-stellar are restored from the original pixels through
//      a dilated, feathered mask. The stars image is corrected to match.
//
// Nothing is computed until Preview is pressed. The preview then shows the
// whole image at 1/3 resolution; drag a box on it to see that area at full
// resolution. After that the preview follows the controls. Settings are
// remembered between runs.
//
// Requirements
//   - A plate-solved image (any one of the selected views).
//   - Gaia XPSD database files configured in the Gaia process. The full
//     DR3 / EDR3 set is strongly preferred; DR3/SP is much shallower.
//
// Workflow
//   A) Clone your image, run the star remover on the clone, then select
//      Starless + Original.            (most accurate)
//   B) Run the star remover in place with a stars image, then select
//      Starless + Stars.               (original is reconstructed)
//
// Changes to existing images are undoable (Ctrl+Z).
// ============================================================================

#feature-id    Utilities > StellarShield
#feature-info  Restores galaxies and other non-stellar objects removed by a \
               star-removal tool, using Gaia cross-matching and morphology.

#include <pjsr/Sizer.jsh>
#include <pjsr/TextAlign.jsh>
#include <pjsr/StdButton.jsh>
#include <pjsr/StdIcon.jsh>
#include <pjsr/NumericControl.jsh>
#include <pjsr/UndoFlag.jsh>
#include <pjsr/DataType.jsh>
#include <pjsr/ColorSpace.jsh>
#include <pjsr/SampleType.jsh>

#define TITLE   "StellarShield"
#define VERSION "1.3.2"
#define SETTINGS "StellarShield/"

// Blob classes
#define CLS_STAR     1
#define CLS_PROTECT  2
#define CLS_GUARD    3
#define CLS_TOOLARGE 4

// Faint unmatched policy
#define FAINT_MORPH   0
#define FAINT_PROTECT 1
#define FAINT_STAR    2

// Preview
#define OVERVIEW_STEP 3     // the whole-image preview shows every 3x3 block as one pixel

// ----------------------------------------------------------------------------
// Parameters
// ----------------------------------------------------------------------------

var P = {
   starlessView: null,
   originalView: null,
   starsView:    null,
   gaiaRelease:  0,       // Gaia database files: 0 auto (deepest available), 1 DR3, 2 DR3/SP
   blend:        0,       // how the stars image combines with the starless: 0 auto, 1 added, 2 screened (StarX Unscreen)
   linear:       true,    // enables the photometric tests

   kSigma:       2.0,     // detection threshold, in noise sigmas
   minArea:      6,       // px
   maxArea:      5000,    // px; larger blobs are left alone
   matchRadius:  2,       // px, Gaia position -> blob

   completeMag:  19.5,    // Gaia assumed complete brighter than this
   excessMag:    0.75,    // blob brighter than its Gaia stars by this -> non-stellar
   faintPolicy:  FAINT_MORPH,
   sizeSigma:    3.0,     // morphology outlier threshold
   minEll:       0.35,    // minimum ellipticity to call something elongated

   guardMag:     11.0,    // stars brighter than this get a guard zone
   guardGrow:    0.25,    // guard zone = star footprint radius grown by this fraction
   armGrow:      0.50,    // extra reach along diffraction spikes
   guardPad:     3,       // px added to every guard radius

   dilate:       3,       // px
   feather:      4,       // px, width of the soft edge
   maskOnly:     false,

   previewMode:  0,       // 0 result, 1 starless before, 2 mask, 3 removed
   view:         null     // null = whole image; else { x0, y0, x1, y1 } in image px
};

// Settings remembered between runs: b = boolean, d = number
var PERSIST = {
   gaiaRelease: 'd', blend: 'd', linear: 'b', kSigma: 'd', minArea: 'd', maxArea: 'd', completeMag: 'd',
   excessMag: 'd', faintPolicy: 'd', sizeSigma: 'd', minEll: 'd', guardMag: 'd',
   dilate: 'd', feather: 'd', maskOnly: 'b', previewMode: 'd'
};

// ----------------------------------------------------------------------------
// Small numeric helpers (plain arrays, PJSR-safe)
// ----------------------------------------------------------------------------

function medianOf( a )
{
   var n = a.length;
   if ( n == 0 )
      return 0;
   a.sort( function( x, y ) { return x - y; } );
   return (n & 1) ? a[n >> 1] : 0.5*(a[(n >> 1) - 1] + a[n >> 1]);
}

// k-th smallest element of a (Float32Array or Array; reordered in place).
function selectK( a, k )
{
   var lo = 0, hi = a.length - 1;
   while ( hi > lo )
   {
      var pv = a[(lo + hi) >> 1], i = lo, j = hi;
      while ( i <= j )
      {
         while ( a[i] < pv ) ++i;
         while ( a[j] > pv ) --j;
         if ( i <= j )
         {
            var t = a[i]; a[i] = a[j]; a[j] = t;
            ++i; --j;
         }
      }
      if ( k <= j )
         hi = j;
      else if ( k >= i )
         lo = i;
      else
         break;
   }
   return a[k];
}

// Fast median of a plain array of numbers (reordered in place).
function fastMedian( a )
{
   var n = a.length;
   if ( n == 0 )
      return 0;
   var f = (a instanceof Float32Array) ? a : Float32Array.from ? Float32Array.from( a ) : new Float32Array( a );
   var m = selectK( f, n >> 1 );
   if ( n & 1 )
      return m;
   var lowMax = f[0];
   for ( var i = 1; i < (n >> 1); ++i )
      if ( f[i] > lowMax )
         lowMax = f[i];
   return 0.5*(lowMax + m);
}

function sigmaMAD( a, med )
{
   var d = new Array( a.length );
   for ( var i = 0; i < a.length; ++i )
      d[i] = Math.abs( a[i] - med );
   return 1.4826*medianOf( d );
}

// Median / sigma of a quantity in bins of a key (magnitude).
function BinnedStat( step )
{
   this.step = step;
   this.raw = {};      // bin index -> array of values
   this.bins = [];     // finalized: { k, med, sig, n }

   this.add = function( key, value )
   {
      var k = Math.floor( key/this.step );
      if ( this.raw[k] === undefined )
         this.raw[k] = [];
      this.raw[k].push( value );
   };

   this.finalize = function( minCount, sigFloorFrac )
   {
      this.bins = [];
      for ( var k in this.raw )
      {
         var v = this.raw[k];
         if ( v.length < minCount )
            continue;
         var med = medianOf( v );
         var sig = sigmaMAD( v, med );
         var floor = Math.abs( med )*sigFloorFrac;
         if ( sig < floor )
            sig = floor;
         this.bins.push( { k: parseInt( k ), med: med, sig: sig, n: v.length } );
      }
      this.raw = {};
      return this.bins.length;
   };

   // Nearest populated bin, or null.
   this.lookup = function( key )
   {
      var k = Math.floor( key/this.step ), best = null, bd = 1e30;
      for ( var i = 0; i < this.bins.length; ++i )
      {
         var d = Math.abs( this.bins[i].k - k );
         if ( d < bd )
         {
            bd = d;
            best = this.bins[i];
         }
      }
      return best;
   };
}

// Robust noise sigma from differences of horizontally adjacent pixels.
function estimateNoise( a, w, h )
{
   var N = w*h;
   var step = Math.max( 1, Math.floor( N/200000 ) ) | 1;
   var d = [];
   for ( var i = 0; i + 1 < N; i += step )
      if ( (i % w) != w - 1 )
         d.push( Math.abs( a[i + 1] - a[i] ) );
   return fastMedian( d )/0.9539;
}

function boxBlur3( src, tmp, dst, w, h )
{
   var x, y, r;
   for ( y = 0; y < h; ++y )
   {
      r = y*w;
      for ( x = 0; x < w; ++x )
      {
         var xm = (x > 0) ? x - 1 : x, xp = (x < w - 1) ? x + 1 : x;
         tmp[r + x] = (src[r + xm] + src[r + x] + src[r + xp])/3;
      }
   }
   for ( y = 0; y < h; ++y )
   {
      r = y*w;
      var rm = (y > 0) ? r - w : r, rp = (y < h - 1) ? r + w : r;
      for ( x = 0; x < w; ++x )
         dst[r + x] = (tmp[rm + x] + tmp[r + x] + tmp[rp + x])/3;
   }
}

// In-place chamfer distance transform. On entry d = 0 at seed pixels and a
// large value elsewhere; on exit d ~ Euclidean distance in pixels (within ~4%).
function chamferDT( d, w, h )
{
   var a = 0.9619, b = 1.3604, x, y, i, v;
   for ( y = 0; y < h; ++y )
      for ( x = 0, i = y*w; x < w; ++x, ++i )
      {
         v = d[i];
         if ( v == 0 )
            continue;
         if ( x > 0 && d[i - 1] + a < v ) v = d[i - 1] + a;
         if ( y > 0 )
         {
            if ( d[i - w] + a < v ) v = d[i - w] + a;
            if ( x > 0 && d[i - w - 1] + b < v ) v = d[i - w - 1] + b;
            if ( x < w - 1 && d[i - w + 1] + b < v ) v = d[i - w + 1] + b;
         }
         d[i] = v;
      }
   for ( y = h - 1; y >= 0; --y )
      for ( x = w - 1, i = y*w + w - 1; x >= 0; --x, --i )
      {
         v = d[i];
         if ( v == 0 )
            continue;
         if ( x < w - 1 && d[i + 1] + a < v ) v = d[i + 1] + a;
         if ( y < h - 1 )
         {
            if ( d[i + w] + a < v ) v = d[i + w] + a;
            if ( x < w - 1 && d[i + w + 1] + b < v ) v = d[i + w + 1] + b;
            if ( x > 0 && d[i + w - 1] + b < v ) v = d[i + w - 1] + b;
         }
         d[i] = v;
      }
}

// ----------------------------------------------------------------------------
// Detection: connected blobs above threshold
//   D = residual, S = smoothed residual, T = threshold
//   labels: 0 = background, -1 = tiny speck, >0 = index into blobs[]
// ----------------------------------------------------------------------------

function detectBlobs( D, S, labels, w, h, T, minArea )
{
   var N = w*h;
   var blobs = [ null ];
   var stack = new Int32Array( 1 << 16 );
   var members = new Int32Array( Math.max( 1, minArea ) );

   for ( var i = 0; i < N; ++i )
   {
      if ( labels[i] != 0 || !(S[i] > T) )
         continue;

      var id = blobs.length;
      var sp = 0;
      stack[sp++] = i;
      labels[i] = id;

      var ox = i % w, oy = (i - ox)/w;
      var area = 0, flux = 0, peak = 0;
      var sw = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
      var x0 = ox, x1 = ox, y0 = oy, y1 = oy;

      while ( sp > 0 )
      {
         var p = stack[--sp];
         var px = p % w, py = (p - px)/w;
         var d = D[p], wt = S[p];
         if ( d > 0 )
            flux += d;
         if ( d > peak )
            peak = d;
         var dx = px - ox, dy = py - oy;
         sw += wt; sx += wt*dx; sy += wt*dy;
         sxx += wt*dx*dx; syy += wt*dy*dy; sxy += wt*dx*dy;
         if ( px < x0 ) x0 = px;
         if ( px > x1 ) x1 = px;
         if ( py < y0 ) y0 = py;
         if ( py > y1 ) y1 = py;
         if ( area < minArea )
            members[area] = p;
         ++area;

         for ( var ny = py - 1; ny <= py + 1; ++ny )
         {
            if ( ny < 0 || ny >= h )
               continue;
            for ( var nx = px - 1; nx <= px + 1; ++nx )
            {
               if ( nx < 0 || nx >= w )
                  continue;
               var q = ny*w + nx;
               if ( labels[q] == 0 && S[q] > T )
               {
                  labels[q] = id;
                  if ( sp == stack.length )
                  {
                     var ns = new Int32Array( stack.length*2 );
                     ns.set( stack );
                     stack = ns;
                  }
                  stack[sp++] = q;
               }
            }
         }
      }

      if ( area < minArea )
      {
         for ( var j = 0; j < area; ++j )
            labels[members[j]] = -1;
         continue;
      }

      var mx = sx/sw, my = sy/sw;
      var vxx = Math.max( sxx/sw - mx*mx, 0 ) + 1/12;
      var vyy = Math.max( syy/sw - my*my, 0 ) + 1/12;
      var vxy = sxy/sw - mx*my;
      var tr = 0.5*(vxx + vyy);
      var disc = Math.sqrt( Math.max( 0, tr*tr - (vxx*vyy - vxy*vxy) ) );
      var a = Math.sqrt( tr + disc );
      var b = Math.sqrt( Math.max( tr - disc, 1e-6 ) );

      blobs.push( {
         id: id, area: area, flux: Math.max( flux, 1e-12 ), peak: peak,
         cx: ox + mx + 0.5, cy: oy + my + 0.5,
         x0: x0, y0: y0, x1: x1, y1: y1,
         size: Math.sqrt( a*b ), ell: 1 - b/a,
         pa: 0.5*Math.atan2( 2*vxy, vxx - vyy ),
         inst: 0, estG: NaN,
         nG: 0, gFlux: 0, gBright: 99, nPM: 0,
         cls: 0, reason: ""
      } );
   }

   for ( var k = 1; k < blobs.length; ++k )
      blobs[k].inst = -2.5*Math.log( blobs[k].flux )/Math.LN10;

   return blobs;
}

// ----------------------------------------------------------------------------
// Gaia cross-match
//   sources: array of { x, y, G, hasPM } in image pixel coordinates
// ----------------------------------------------------------------------------

function matchGaia( sources, blobs, labels, w, h, radius )
{
   var matched = 0;
   for ( var i = 0; i < sources.length; ++i )
   {
      var s = sources[i];
      var ix = Math.floor( s.x ), iy = Math.floor( s.y );
      s.blob = 0;
      if ( ix < 0 || iy < 0 || ix >= w || iy >= h )
         continue;
      var best = 0, bd = 1e30;
      for ( var y = Math.max( 0, iy - radius ); y <= Math.min( h - 1, iy + radius ); ++y )
         for ( var x = Math.max( 0, ix - radius ); x <= Math.min( w - 1, ix + radius ); ++x )
         {
            var L = labels[y*w + x];
            if ( L > 0 )
            {
               var d = (x - ix)*(x - ix) + (y - iy)*(y - iy);
               if ( d < bd )
               {
                  bd = d;
                  best = L;
               }
            }
         }
      s.blob = best;
      if ( best > 0 )
      {
         var b = blobs[best];
         ++b.nG;
         b.gFlux += Math.pow( 10, -0.4*s.G );
         if ( s.G < b.gBright )
            b.gBright = s.G;
         if ( s.hasPM )
            ++b.nPM;
         ++matched;
      }
   }
   return matched;
}

// ----------------------------------------------------------------------------
// Bright-star guard zones
//   Each bright star gets a polar profile of its own removed footprint
//   (radius reached in every direction), so the zone is round where the star
//   is round and reaches out only along its diffraction spikes.
// ----------------------------------------------------------------------------

#define GUARD_BINS 180

function buildGuards( blobs, labels, w, h, P )
{
   var guards = [];
   for ( var i = 1; i < blobs.length; ++i )
   {
      var b = blobs[i];
      if ( b.nG == 0 || b.gBright >= P.guardMag )
         continue;

      var prof = new Float32Array( GUARD_BINS ), k;
      for ( var y = b.y0; y <= b.y1; ++y )
         for ( var x = b.x0; x <= b.x1; ++x )
            if ( labels[y*w + x] == b.id )
            {
               var dx = x + 0.5 - b.cx, dy = y + 0.5 - b.cy;
               var r = Math.sqrt( dx*dx + dy*dy );
               k = Math.floor( (Math.atan2( dy, dx ) + Math.PI)/(2*Math.PI)*GUARD_BINS ) % GUARD_BINS;
               if ( r > prof[k] )
                  prof[k] = r;
            }

      var tmp = [];
      for ( k = 0; k < GUARD_BINS; ++k )
         tmp.push( prof[k] );
      var core = medianOf( tmp );

      var g = { cx: b.cx, cy: b.cy, rmax: 0, core: 0,
                lim: new Float32Array( GUARD_BINS ), foot: new Float32Array( GUARD_BINS ), arm: [] };
      var mm = [], arm = [];
      for ( k = 0; k < GUARD_BINS; ++k )
      {
         mm.push( Math.max( prof[(k + GUARD_BINS - 1) % GUARD_BINS], prof[k], prof[(k + 1) % GUARD_BINS] ) );
         arm.push( mm[k] > Math.max( 1.5*core, core + 6 ) );
      }
      // A diffraction spike is narrow. A wide run of long bins is nebulosity
      // (or a galaxy) merged into the star's blob, and is not part of the star.
      for ( k = 0; k < GUARD_BINS; ++k )
         if ( arm[k] && !arm[(k + GUARD_BINS - 1) % GUARD_BINS] )
         {
            var len = 0;
            while ( len < GUARD_BINS && arm[(k + len) % GUARD_BINS] )
               ++len;
            if ( len > 12 )
               for ( var q = 0; q < len; ++q )
                  arm[(k + q) % GUARD_BINS] = false;
         }
      g.core = core*(1 + P.guardGrow) + P.guardPad;
      for ( k = 0; k < GUARD_BINS; ++k )
      {
         var m = arm[k] ? mm[k] : Math.min( mm[k], core );
         g.arm.push( arm[k] );
         g.foot[k] = m + P.guardPad;
         g.lim[k] = arm[k] ? m*(1 + P.armGrow) + P.guardPad
                           : core*(1 + P.guardGrow) + P.guardPad;
         if ( g.lim[k] > g.rmax )
            g.rmax = g.lim[k];
      }
      guards.push( g );
   }
   return guards;
}

// Is this unmatched blob a fragment of a bright star (halo piece, spike piece)?
function inGuard( b, guards )
{
   for ( var i = 0; i < guards.length; ++i )
   {
      var g = guards[i];
      var dx = b.cx - g.cx, dy = b.cy - g.cy;
      if ( Math.abs( dx ) > g.rmax || Math.abs( dy ) > g.rmax )
         continue;
      var r = Math.sqrt( dx*dx + dy*dy );
      var phi = Math.atan2( dy, dx );
      var k = Math.floor( (phi + Math.PI)/(2*Math.PI)*GUARD_BINS ) % GUARD_BINS;
      if ( r > g.lim[k] )
         continue;
      if ( !g.arm[k] || r <= g.foot[k] )
         return true;
      // Beyond the end of a spike: only roundish specks and pieces that point
      // along the spike are treated as spike fragments.
      var d = Math.abs( b.pa - phi ) % Math.PI;
      if ( d > 0.5*Math.PI )
         d = Math.PI - d;
      if ( b.ell < 0.3 || d < 0.44 )
         return true;
   }
   return false;
}

// ----------------------------------------------------------------------------
// Classification
// ----------------------------------------------------------------------------

function classifyBlobs( blobs, P, effCompleteMag, log, labels, w, h )
{
   var i, b;

   // --- Photometric zero point from single, well-behaved Gaia stars ---------
   var zpAll = [], zpStat = new BinnedStat( 1.0 ), zp = NaN;
   var photometric = P.linear;
   if ( photometric )
   {
      for ( i = 1; i < blobs.length; ++i )
      {
         b = blobs[i];
         if ( b.nG == 1 && b.nPM == 1 && b.area <= P.maxArea )
         {
            var z = b.gBright - b.inst;
            zpStat.add( b.gBright, z );
            if ( b.gBright > 13 && b.gBright < 18.5 )
               zpAll.push( z );
         }
      }
      if ( zpAll.length < 20 || zpStat.finalize( 8, 0 ) == 0 )
      {
         photometric = false;
         log( "** Too few calibration stars: photometric tests disabled." );
      }
      else
      {
         zp = medianOf( zpAll );
         log( "Zero point: G = " + zp.toFixed( 2 ) + " - 2.5*log10(flux), from " + zpAll.length + " stars." );
      }
   }

   // --- Flux excess over Gaia (also removes blends from the locus) ----------
   function fluxExcess( b )
   {
      if ( !photometric || b.nG == 0 || b.gBright < P.guardMag )
         return false;
      var bin = zpStat.lookup( b.gBright );
      if ( bin == null )
         return false;
      var predG = -2.5*Math.log( b.gFlux )/Math.LN10;
      var blobG = bin.med + b.inst;
      return predG - blobG > Math.max( P.excessMag, 3*bin.sig );
   }

   // --- Stellar locus: size and ellipticity vs instrumental magnitude -------
   var sizeStat = new BinnedStat( 0.5 ), ellStat = new BinnedStat( 0.5 );
   for ( i = 1; i < blobs.length; ++i )
   {
      b = blobs[i];
      b.reason = "";
      b.estG = NaN;
      b.excess = fluxExcess( b );
      if ( photometric )
         b.estG = zp + b.inst;
      if ( b.nG == 1 && b.nPM == 1 && !b.excess && b.area <= P.maxArea )
      {
         sizeStat.add( b.inst, b.size );
         ellStat.add( b.inst, b.ell );
      }
   }
   var haveLocus = sizeStat.finalize( 10, 0.05 ) > 0;
   ellStat.finalize( 10, 0 );
   if ( !haveLocus )
      log( "** Too few Gaia stars for a stellar locus: size test disabled." );

   function extended( b, k )
   {
      var e = ellStat.lookup( b.inst );
      var eThr = (e != null) ? Math.max( P.minEll, e.med + k*Math.max( e.sig, 0.03 ) ) : P.minEll;
      if ( b.ell > eThr )
         return "elongated";
      if ( haveLocus )
      {
         var s = sizeStat.lookup( b.inst );
         if ( (b.size - s.med)/s.sig > k )
            return "larger than stars";
      }
      return "";
   }

   // --- Guard zones around bright stars (spikes, halo fragments) ------------
   var guards = buildGuards( blobs, labels, w, h, P );

   // --- Decide --------------------------------------------------------------
   var count = { star: 0, protect: 0, guard: 0, large: 0, guards: guards };
   for ( i = 1; i < blobs.length; ++i )
   {
      b = blobs[i];
      var why = "";
      b.compact = extended( b, P.sizeSigma ) == "";

      if ( b.area > P.maxArea )
      {
         b.cls = CLS_TOOLARGE;
         ++count.large;
         continue;
      }

      if ( b.nG == 0 )
      {
         if ( inGuard( b, guards ) )
         {
            b.cls = CLS_GUARD;
            ++count.guard;
            continue;
         }

         if ( photometric && b.estG < effCompleteMag )
            why = "bright, not in Gaia";
         else if ( P.faintPolicy == FAINT_PROTECT )
            why = "faint, not in Gaia";
         else if ( P.faintPolicy == FAINT_MORPH )
         {
            why = extended( b, P.sizeSigma );
            if ( why != "" )
               why = "not in Gaia, " + why;
         }
      }
      else
      {
         if ( b.excess )
            why = "brighter than its Gaia source";
         else if ( b.nG == 1 && b.nPM == 0 && b.gBright >= P.guardMag )
         {
            // Gaia source without proper motion/parallax: often a galaxy core.
            why = extended( b, P.sizeSigma + 1 );
            if ( why != "" )
               why = "Gaia source without PM, " + why;
         }
      }

      if ( why != "" )
      {
         b.cls = CLS_PROTECT;
         b.reason = why;
         ++count.protect;
      }
      else
      {
         b.cls = CLS_STAR;
         ++count.star;
      }
   }
   return count;
}

// ----------------------------------------------------------------------------
// Mask: protected blobs grown by a round distance transform with a soft
// cosine edge, never over star footprints
// ----------------------------------------------------------------------------

//   Works on the rectangle (x0,y0,rw,rh) of the image; out has rw*rh elements.
function buildMaskRect( out, x0, y0, rw, rh, labels, blobs, w, dilate, feather )
{
   var x, y, i, j, L, d, v;
   for ( y = 0; y < rh; ++y )
      for ( x = 0, i = (y0 + y)*w + x0, j = y*rw; x < rw; ++x, ++i, ++j )
      {
         L = labels[i];
         out[j] = (L > 0 && blobs[L].cls == CLS_PROTECT) ? 0 : 1e9;
      }
   chamferDT( out, rw, rh );
   for ( y = 0; y < rh; ++y )
      for ( x = 0, i = (y0 + y)*w + x0, j = y*rw; x < rw; ++x, ++i, ++j )
      {
         d = out[j];
         if ( d <= dilate )
            v = 1;
         else if ( feather <= 0 || d >= dilate + feather )
            v = 0;
         else
            v = 0.5*(1 + Math.cos( Math.PI*(d - dilate)/feather ));
         if ( v > 0 )
         {
            L = labels[i];
            if ( L > 0 && (blobs[L].cls == CLS_STAR || blobs[L].cls == CLS_TOOLARGE) )
               v = 0;
         }
         out[j] = v;
      }
}

// ----------------------------------------------------------------------------
// PixInsight glue: astrometry and Gaia
// ----------------------------------------------------------------------------

function hasSolution( win )
{
   var s = win.hasAstrometricSolution;
   if ( typeof s == "function" )
      s = s.call( win );
   return s ? true : false;
}

function pixToSky( win, x, y )
{
   var p = win.imageToCelestial( x, y );
   if ( p == null )
      throw new Error( "imageToCelestial failed at (" + x + "," + y + ")" );
   return p; // x = RA (deg), y = Dec (deg)
}

function skyToPix( win, ra, dec )
{
   return win.celestialToImage( ra, dec ); // Point or null
}

function angSepDeg( ra1, dec1, ra2, dec2 )
{
   var r = Math.PI/180;
   var s1 = Math.sin( 0.5*(dec2 - dec1)*r ), s2 = Math.sin( 0.5*(ra2 - ra1)*r );
   var a = s1*s1 + Math.cos( dec1*r )*Math.cos( dec2*r )*s2*s2;
   return 2*Math.asin( Math.min( 1, Math.sqrt( a ) ) )/r;
}

// Returns { sources: [{x,y,G,hasPM}], release: name, faintMag: p98 of G }
function queryGaia( win, w, h, choice )
{
   var c = pixToSky( win, 0.5*w, 0.5*h );
   var radius = 0;
   var corners = [ [0, 0], [w, 0], [0, h], [w, h] ];
   for ( var i = 0; i < 4; ++i )
   {
      var q = pixToSky( win, corners[i][0], corners[i][1] );
      radius = Math.max( radius, angSepDeg( c.x, c.y, q.x, q.y ) );
   }
   radius *= 1.02;
   console.writeln( format( "Field center: RA %.5f deg, Dec %+.5f deg, radius %.3f deg", c.x, c.y, radius ) );

   // Deepest catalogs first. DR3/SP only holds sources with BP/RP spectra.
   // Gaia process dataRelease values: 1 DR2, 2 EDR3, 3 DR3, 4 DR3/SP.
   // Auto tries the deepest first; DR3/SP only holds sources with BP/RP spectra.
   var releases = (choice == 1) ? [ [3, "DR3"] ] :
                  (choice == 2) ? [ [4, "DR3/SP"] ] :
                  [ [3, "DR3"], [2, "EDR3"], [4, "DR3/SP"], [1, "DR2"] ];
   var S = null, relName = "";
   for ( var r = 0; r < releases.length && S == null; ++r )
   {
      try
      {
         var g = new Gaia;
         g.command = "search";
         g.centerRA = c.x;
         g.centerDec = c.y;
         g.radius = radius;
         g.magnitudeLow = -1.5;
         g.magnitudeHigh = 26;
         g.sourceLimit = 4000000;
         g.requiredFlags = 0;
         g.inclusionFlags = 0;
         g.exclusionFlags = 0;
         g.verbosity = 0;
         g.dataRelease = releases[r][0];
         g.generateTextOutput = false;
         g.generateBinaryOutput = false;
         g.executeGlobal();
         if ( g.sources && g.sources.length > 0 )
         {
            S = g.sources;
            relName = releases[r][1];
         }
      }
      catch ( e )
      {
         // database for this release not configured: try the next one
      }
   }
   if ( S == null )
      throw new Error( (choice == 1 || choice == 2) ?
                       "No Gaia " + releases[0][1] + " sources found. Check that the " + releases[0][1] +
                       " XPSD database files are configured in the Gaia process, or choose another Gaia catalog." :
                       "No Gaia sources found. Check that XPSD database files are configured in the Gaia process." );

   console.writeln( "Gaia " + relName + ": " + S.length + " sources. First row: " + S[0] );

   // Row layout: ra, dec, parx, pmra, pmdec, magG, magBP, magRP, flags
   var years = (new Date).getTime()/3.15576e10 + 1970 - 2016.0;   // since Gaia DR3 epoch
   var out = [], mags = [];
   for ( var i = 0; i < S.length; ++i )
   {
      var row = S[i];
      var ra = row[0], dec = row[1], pmra = row[3], pmdec = row[4], G = row[5], flags = row[8];
      var hasPM = !(flags & 1) && isFinite( pmra ) && isFinite( pmdec ) && (pmra != 0 || pmdec != 0);
      if ( hasPM )
      {
         dec += pmdec*years/3.6e6;
         ra += pmra*years/3.6e6/Math.max( 0.01, Math.cos( dec*Math.PI/180 ) );
      }
      var p = skyToPix( win, ra, dec );
      if ( p == null )
         continue;
      if ( p.x < 0 || p.y < 0 || p.x >= w || p.y >= h )
         continue;
      out.push( { x: p.x, y: p.y, G: G, hasPM: hasPM } );
      mags.push( G );
   }
   mags.sort( function( a, b ) { return a - b; } );
   var faint = mags.length ? mags[Math.floor( 0.98*(mags.length - 1) )] : 0;
   return { sources: out, release: relName, faintMag: faint };
}

// ----------------------------------------------------------------------------
// Engine
//   Cached pipeline. level = how far the cached results are valid:
//   0 nothing, 1 residual + Gaia loaded, 2 objects detected, 3 classified.
//   The full-frame protection mask is cached separately and rebuilt when the
//   classification or the mask settings change.
// ----------------------------------------------------------------------------

function mtf( m, x )
{
   if ( x <= 0 ) return 0;
   if ( x >= 1 ) return 1;
   return (m - 1)*x/((2*m - 1)*x - m);
}

function Engine()
{
   this.level = 0;
   this.hasRun = false;
   this.N = 0;
   this.classStamp = 0;
   this.maskKey = "";

   this.invalidate = function( keep )
   {
      if ( this.level > keep )
         this.level = keep;
   };

   // --- 1: residual removed by the star remover, noise, Gaia -----------------
   this.load = function()
   {
      var SL = P.starlessView, OR = P.originalView, ST = P.starsView;
      var haveOrig = OR != null && !OR.isNull, haveStars = ST != null && !ST.isNull;
      if ( SL == null || SL.isNull )
         throw new Error( "Select the starless image." );
      if ( !haveOrig && !haveStars )
         throw new Error( "Select the original image, the stars image, or both." );

      var img = SL.image, w = img.width, h = img.height, nc = img.numberOfNominalChannels, N = w*h;
      var others = [], i, c;
      if ( haveOrig ) others.push( OR );
      if ( haveStars ) others.push( ST );
      for ( i = 0; i < others.length; ++i )
      {
         var o = others[i].image;
         if ( o.width != w || o.height != h || o.numberOfNominalChannels != nc )
            throw new Error( "All images must have the same dimensions and channels: " + others[i].id );
      }

      var awin = null, cand = [ SL ].concat( others );
      for ( i = 0; i < cand.length && awin == null; ++i )
         if ( hasSolution( cand[i].window ) )
            awin = cand[i].window;
      if ( awin == null )
         throw new Error( "None of the selected images has an astrometric solution. Run ImageSolver first." );

      this.SL = SL; this.OR = OR; this.ST = ST;
      this.haveOrig = haveOrig; this.haveStars = haveStars;
      this.w = w; this.h = h; this.nc = nc;
      if ( N != this.N )
      {
         this.A = new Float32Array( N );
         this.B = new Float32Array( N );
         this.D = new Float32Array( N );
         this.M = new Float32Array( N );
         this.labels = new Int32Array( N );
         this.N = N;
      }
      this.maskKey = "";

      var A = this.A, B = this.B, D = this.D;
      var rect = new Rect( 0, 0, w, h ), SRC = haveOrig ? OR : ST;

      // How does the stars image combine with the starless image? Added
      // (original = starless + stars) or screened (StarX "Unscreen":
      // original = 1 - (1 - starless)(1 - stars)). With added images the sum
      // can never exceed 1, because the original cannot; with screened ones
      // every star core near 1 on a non-black background pushes it past 1.
      // Where neither happens the two formulas differ negligibly.
      this.screen = false;
      this.blendNote = "";
      if ( !haveOrig )
      {
         if ( P.blend == 2 )
         {
            this.screen = true;
            this.blendNote = "screened (set by you)";
         }
         else if ( P.blend == 1 )
            this.blendNote = "added (set by you)";
         else
         {
            var over = 0;
            for ( c = 0; c < nc; ++c )
            {
               SL.image.getSamples( A, rect, c );
               ST.image.getSamples( B, rect, c );
               for ( i = 0; i < N; ++i )
                  if ( A[i] + B[i] > 1.006 )
                     ++over;
            }
            this.screen = over > Math.max( 20, 2e-5*N*nc );
            this.blendNote = (this.screen ? "screened" : "added") + " (detected: " + over +
                             " pixels where starless + stars exceeds 1)";
         }
         console.writeln( "<end><cbr>Stars image: " + this.blendNote );
      }

      for ( i = 0; i < N; ++i )
         D[i] = 0;
      for ( c = 0; c < nc; ++c )
      {
         SL.image.getSamples( A, rect, c );
         SRC.image.getSamples( B, rect, c );
         if ( c == 0 )
            this.sigma = estimateNoise( A, w, h );
         if ( haveOrig )
            for ( i = 0; i < N; ++i ) D[i] += (B[i] - A[i])/nc;
         else if ( this.screen )
            for ( i = 0; i < N; ++i ) D[i] += B[i]*(1 - A[i])/nc;
         else
            for ( i = 0; i < N; ++i ) D[i] += B[i]/nc;
      }
      this.S = B;
      boxBlur3( D, A, this.S, w, h );
      console.writeln( format( "<end><cbr><br>Noise sigma: %.3e", this.sigma ) );

      this.gaia = queryGaia( awin, w, h, P.gaiaRelease );
      this.effComplete = Math.min( P.completeMag, this.gaia.faintMag - 1.0 );
      if ( this.gaia.faintMag - 1.0 < 19 )
         console.warningln( format( "** Gaia %ls only reaches G ~ %.1f here. Install the full DR3 database for deeper matching.",
                                    this.gaia.release, this.gaia.faintMag ) );
      this.level = 1;
   };

   // --- 2: detect removed objects and match them to Gaia ---------------------
   this.detect = function()
   {
      var labels = this.labels;
      for ( var i = 0; i < this.N; ++i )
         labels[i] = 0;
      this.blobs = detectBlobs( this.D, this.S, labels, this.w, this.h, P.kSigma*this.sigma, P.minArea );
      this.matched = matchGaia( this.gaia.sources, this.blobs, labels, this.w, this.h, P.matchRadius );
      this.level = 2;
   };

   // --- 3: classify ----------------------------------------------------------
   this.classify = function()
   {
      this.effComplete = Math.min( P.completeMag, this.gaia.faintMag - 1.0 );
      this.count = classifyBlobs( this.blobs, P, this.effComplete, function( s ) {}, this.labels, this.w, this.h );
      this.guards = this.count.guards;
      ++this.classStamp;
      this.level = 3;
   };

   this.ensure = function( status )
   {
      if ( this.level < 1 ) { status( "Reading images and querying Gaia..." ); this.load(); }
      if ( this.level < 2 ) { status( "Detecting removed objects..." ); this.detect(); }
      if ( this.level < 3 ) { status( "Classifying..." ); this.classify(); }
   };

   this.summary = function()
   {
      var c = this.count;
      return "Gaia " + this.gaia.release + "  |  removed objects: " + (this.blobs.length - 1) + "  |  stars: " + c.star +
             "  |  restored: " + c.protect +
             (this.haveOrig ? "" : "  |  stars image " + (this.screen ? "screened" : "added"));
   };

   // Full-frame protection mask, cached.
   this.fullMask = function()
   {
      var key = this.classStamp + "|" + P.dilate + "|" + P.feather;
      if ( key != this.maskKey )
      {
         buildMaskRect( this.M, 0, 0, this.w, this.h, this.labels, this.blobs, this.w, P.dilate, P.feather );
         this.maskKey = key;
      }
      return this.M;
   };

   // Preview bitmap. With P.view == null: the whole image, each OVERVIEW_STEP x
   // OVERVIEW_STEP block averaged into one pixel. Otherwise the P.view region
   // at full resolution. Returns { bitmap, x0, y0, step }.
   this.render = function()
   {
      var w = this.w, h = this.h, nc = this.nc;
      var x0 = 0, y0 = 0, x1 = w, y1 = h, s = OVERVIEW_STEP;
      if ( P.view != null )
      {
         x0 = P.view.x0; y0 = P.view.y0; x1 = P.view.x1; y1 = P.view.y1;
         s = 1;
      }
      var cw = x1 - x0, ch = y1 - y0, n = cw*ch;
      var ow = Math.ceil( cw/s ), oh = Math.ceil( ch/s ), on = ow*oh;
      var rect = new Rect( x0, y0, x1, y1 ), i, c, x, y;
      var mode = P.previewMode, color = nc >= 3 && mode < 2;
      var M = this.fullMask();

      // block sizes (edge blocks can be smaller)
      var bx = new Int32Array( cw ), by = new Int32Array( ch ), inv = new Float32Array( on );
      for ( x = 0; x < cw; ++x ) bx[x] = Math.floor( x/s );
      for ( y = 0; y < ch; ++y ) by[y] = Math.floor( y/s );
      for ( y = 0; y < ch; ++y )
         for ( x = 0; x < cw; ++x )
            inv[by[y]*ow + bx[x]] += 1;
      for ( i = 0; i < on; ++i )
         inv[i] = 1/inv[i];

      if ( this.rbA == null || this.rbA.length < n )
      {
         this.rbA = new Float32Array( n );
         this.rbB = new Float32Array( n );
      }
      var a = this.rbA.subarray( 0, n ), b = this.rbB.subarray( 0, n );
      var out = new Float32Array( on );
      var img = new Image( ow, oh, color ? 3 : 1, color ? ColorSpace_RGB : ColorSpace_Gray, 32, SampleType_Real );
      var orect = new Rect( 0, 0, ow, oh );

      // screen stretch from the starless image in this region, shared by the
      // result, starless and removed-light views
      var c0 = 0, mb = 0.5, med = 0;
      this.SL.image.getSamples( a, rect, 0 );
      if ( P.linear )
      {
         var smp = [], step = Math.max( 1, Math.floor( n/40000 ) );
         for ( i = 0; i < n; i += step )
            smp.push( a[i] );
         med = medianOf( smp );
         var mad = sigmaMAD( smp, med );
         c0 = Math.max( 0, med - 2.8*mad );
         mb = mtf( 0.25, Math.max( 1e-6, med - c0 ) );
      }
      function show( v )
      {
         if ( !P.linear )
            return (v < 0) ? 0 : (v > 1) ? 1 : v;
         return mtf( mb, (v - c0)/(1 - c0) );
      }
      // block average of a region array (indexed like a and b) into out
      function accumulate( arr )
      {
         for ( i = 0; i < on; ++i )
            out[i] = 0;
         for ( y = 0; y < ch; ++y )
         {
            var row = by[y]*ow, jj = y*cw;
            for ( x = 0; x < cw; ++x, ++jj )
               out[row + bx[x]] += arr[jj];
         }
         for ( i = 0; i < on; ++i )
            out[i] *= inv[i];
      }
      // copy a full-frame array's region into b
      function region( full )
      {
         for ( y = 0; y < ch; ++y )
         {
            var p = (y0 + y)*w + x0, jj = y*cw;
            for ( x = 0; x < cw; ++x )
               b[jj + x] = full[p + x];
         }
      }

      if ( mode == 2 )
      {
         region( M );
         accumulate( b );
         img.setSamples( out, orect, 0 );
      }
      else if ( mode == 3 )
      {
         region( this.D );
         accumulate( b );
         for ( i = 0; i < on; ++i )
            out[i] = show( med + out[i] );
         img.setSamples( out, orect, 0 );
      }
      else
      {
         var SRC = this.haveOrig ? this.OR : this.ST;
         for ( c = 0; c < (color ? 3 : 1); ++c )
         {
            if ( c > 0 )
               this.SL.image.getSamples( a, rect, c );
            if ( mode == 0 )
            {
               SRC.image.getSamples( b, rect, c );
               for ( y = 0; y < ch; ++y )
               {
                  var pm = (y0 + y)*w + x0, jm = y*cw;
                  for ( x = 0; x < cw; ++x )
                  {
                     var m = M[pm + x];
                     if ( m > 0 )
                     {
                        var j2 = jm + x, v = a[j2];
                        a[j2] = v + m*(this.haveOrig ? b[j2] - v : this.screen ? b[j2]*(1 - v) : b[j2]);
                     }
                  }
               }
            }
            accumulate( a );
            for ( i = 0; i < on; ++i )
               out[i] = show( out[i] );
            img.setSamples( out, orect, c );
         }
      }
      return { bitmap: img.render(), x0: x0, y0: y0, step: s };
   };

   // Console report; returns true if there is something to restore.
   this.report = function()
   {
      var c = this.count, i;
      console.writeln( "<end><cbr><br>Removed objects detected: " + (this.blobs.length - 1) +
                       ", Gaia sources in frame: " + this.gaia.sources.length + ", matched: " + this.matched );
      console.writeln( "Stars: " + c.star + "   Restored: " + c.protect +
                       "   Near bright stars (ignored): " + c.guard + "   Too large (ignored): " + c.large );
      var prot = [];
      for ( i = 1; i < this.blobs.length; ++i )
         if ( this.blobs[i].cls == CLS_PROTECT )
            prot.push( this.blobs[i] );
      prot.sort( function( a, b ) { return b.flux - a.flux; } );
      console.writeln( "<br>      x       y   area   est.G   size    ell  reason" );
      for ( i = 0; i < prot.length && i < 300; ++i )
      {
         var b = prot[i];
         console.writeln( format( "%7.1f %7.1f %6d %7.2f %6.2f %6.2f  %ls",
                                  b.cx, b.cy, b.area, isFinite( b.estG ) ? b.estG : 0, b.size, b.ell, b.reason ) );
      }
      if ( prot.length > 300 )
         console.writeln( "... and " + (prot.length - 300) + " more." );
      return prot.length > 0;
   };

   // Full-frame mask image and restoration. Invalidates the cache.
   this.apply = function()
   {
      var w = this.w, h = this.h, nc = this.nc, N = this.N, i, c;
      var SL = this.SL, OR = this.OR, ST = this.ST;
      if ( !this.report() )
      {
         console.writeln( "Nothing to restore." );
         return;
      }

      console.writeln( "<br>Building protection mask..." );
      console.flush();
      var m = this.fullMask(), A = this.A, B = this.B, rect = new Rect( 0, 0, w, h );
      this.level = 0;
      this.hasRun = false;
      this.maskKey = "";

      var mw = new ImageWindow( w, h, 1, 32, true, false, "StellarShield_mask" );
      mw.mainView.beginProcess( UndoFlag_NoSwapFile );
      mw.mainView.image.setSamples( m, rect, 0 );
      mw.mainView.endProcess();
      mw.show();

      if ( P.maskOnly )
      {
         console.writeln( "Mask only: images were not modified." );
         return;
      }

      console.writeln( "Restoring..." );
      console.flush();
      SL.beginProcess();
      if ( this.haveStars )
         ST.beginProcess();
      for ( c = 0; c < nc; ++c )
      {
         SL.image.getSamples( A, rect, c );
         if ( this.haveOrig )
         {
            OR.image.getSamples( B, rect, c );
            for ( i = 0; i < N; ++i )
               if ( m[i] > 0 )
                  A[i] += m[i]*(B[i] - A[i]);
         }
         else
         {
            ST.image.getSamples( B, rect, c );
            for ( i = 0; i < N; ++i )
               if ( m[i] > 0 )
               {
                  var add = this.screen ? B[i]*(1 - A[i]) : B[i];
                  A[i] = Math.min( 1, A[i] + m[i]*add );
               }
         }
         SL.image.setSamples( A, rect, c );

         if ( this.haveStars )
         {
            ST.image.getSamples( B, rect, c );
            for ( i = 0; i < N; ++i )
               if ( m[i] > 0 )
                  B[i] *= 1 - m[i];
            ST.image.setSamples( B, rect, c );
         }
      }
      SL.endProcess();
      if ( this.haveStars )
         ST.endProcess();

      console.writeln( "Done: restored to " + SL.id +
                       (this.haveStars ? " and removed from " + ST.id : "") + ". Use Undo to revert." );
   };
}

// ----------------------------------------------------------------------------
// Remembered settings and image selections
//   Settings persist between PixInsight sessions. Image selections are stored
//   by view id and only reused while that image is still open in the same
//   PixInsight session.
// ----------------------------------------------------------------------------

function sessionToken()
{
   try
   {
      if ( typeof CoreApplication != "undefined" && typeof CoreApplication.pid != "undefined" )
         return "" + CoreApplication.pid;
   }
   catch ( e )
   {
   }
   return "";
}

var VIEW_KEYS = [ "starlessView", "originalView", "starsView" ];

// Factory values of every remembered setting, captured before anything is loaded.
var DEFAULTS = {};
for ( var defKey in PERSIST )
   DEFAULTS[defKey] = P[defKey];

function saveSettings()
{
   for ( var k in PERSIST )
      Settings.write( SETTINGS + k, (PERSIST[k] == 'b') ? DataType_Boolean : DataType_Double, P[k] );
   for ( var i = 0; i < VIEW_KEYS.length; ++i )
   {
      var v = P[VIEW_KEYS[i]];
      Settings.write( SETTINGS + VIEW_KEYS[i], DataType_String, (v != null && !v.isNull) ? v.id : "" );
   }
   Settings.write( SETTINGS + "session", DataType_String, sessionToken() );
}

function loadSettings()
{
   for ( var k in PERSIST )
   {
      var type = (PERSIST[k] == 'b') ? DataType_Boolean : DataType_Double;
      var val = Settings.read( SETTINGS + k, type );
      if ( !Settings.lastReadOK )
         val = Settings.read( "GalaxyGuard/" + k, type );   // settings saved under the old name
      if ( Settings.lastReadOK && val != null )
         P[k] = val;
   }
   P.previewMode = Math.min( 3, Math.max( 0, Math.round( P.previewMode ) ) );
   P.blend = Math.min( 2, Math.max( 0, Math.round( P.blend ) ) );
   P.gaiaRelease = Math.min( 2, Math.max( 0, Math.round( P.gaiaRelease ) ) );
   P.feather = Math.min( 50, Math.max( 0, P.feather ) );

   var session = Settings.read( SETTINGS + "session", DataType_String );
   if ( !Settings.lastReadOK || session != sessionToken() )
      return;
   for ( var i = 0; i < VIEW_KEYS.length; ++i )
   {
      var id = Settings.read( SETTINGS + VIEW_KEYS[i], DataType_String );
      if ( !Settings.lastReadOK || id == null || id == "" )
         continue;
      var v = View.viewById( id );
      if ( v != null && !v.isNull )
         P[VIEW_KEYS[i]] = v;
   }
}

// ----------------------------------------------------------------------------
// Dialog
// ----------------------------------------------------------------------------

function GGDialog( engine )
{
   this.__base__ = Dialog;
   this.__base__();

   var dlg = this;
   var labelW = this.font.width( "Faint sources not in Gaia:" ) + 8;
   this.busy = false;
   this.pending = false;
   this.setters = {};                // key -> function that shows P[key] in its control
   this.cache = {};                  // rendered previews by display mode

   // --- live update plumbing -------------------------------------------------

   // keep = highest engine level still valid after the change:
   //   0 reload images, 1 re-detect, 2 re-classify, 3 re-render only
   this.changed = function( keep, keepCache )
   {
      if ( !keepCache )
         dlg.cache = {};
      engine.invalidate( keep );
      if ( keep == 0 )
      {
         engine.hasRun = false;
         P.view = null;
         dlg.pv.shown = null;
         dlg.pv.message = "Press Preview to analyze the selected images.";
         dlg.pv.update();
         dlg.status.text = "";
         dlg.updateViewButton();
         return;
      }
      if ( engine.hasRun )
         dlg.schedule();
   };

   this.schedule = function()
   {
      try
      {
         dlg.timer.stop();
         dlg.timer.start();
      }
      catch ( e )
      {
         dlg.refresh();
      }
   };

   this.setStatus = function( s )
   {
      dlg.status.text = s;
      processEvents();
   };

   this.refresh = function()
   {
      if ( dlg.busy )
      {
         dlg.pending = true;
         return;
      }
      dlg.busy = true;
      try
      {
         dlg.setStatus( "Working..." );
         engine.ensure( dlg.setStatus );
         if ( dlg.cache[P.previewMode] == null )
         {
            dlg.setStatus( "Rendering..." );
            dlg.cache[P.previewMode] = engine.render();
         }
         dlg.pv.shown = dlg.cache[P.previewMode];
         dlg.pv.scaled = null;
         dlg.setStatus( engine.summary() + "  |  " +
                        ((P.view == null) ? "whole image at 1/" + OVERVIEW_STEP + " resolution: drag a box to zoom in"
                                          : "full resolution: double-click or Full image to go back") );
      }
      catch ( e )
      {
         engine.hasRun = false;
         dlg.pv.shown = null;
         dlg.pv.message = "" + e.message;
         dlg.status.text = "";
         console.criticalln( "*** " + TITLE + ": " + e.message );
      }
      dlg.pv.update();
      dlg.updateViewButton();
      dlg.busy = false;
      if ( dlg.pending )
      {
         dlg.pending = false;
         dlg.schedule();
      }
   };

   this.timer = new Timer;
   this.timer.interval = 0.35;
   this.timer.periodic = false;
   this.timer.onTimeout = function()
   {
      dlg.timer.stop();
      dlg.refresh();
   };

   // --- control helpers ------------------------------------------------------

   function label( text )
   {
      var l = new Label( dlg );
      l.text = text;
      l.textAlignment = TextAlign_Right | TextAlign_VertCenter;
      l.minWidth = labelW;
      return l;
   }

   function viewRow( text, tip, key, preset )
   {
      var vl = new ViewList( dlg );
      vl.getMainViews();
      vl.toolTip = tip;
      if ( P[key] == null || P[key].isNull )
         P[key] = (preset != null && !preset.isNull) ? preset : null;
      if ( P[key] != null )
         vl.currentView = P[key];
      vl.onViewSelected = function( v ) { P[key] = v; dlg.changed( 0 ); };
      var s = new HorizontalSizer;
      s.spacing = 4;
      s.add( label( text ) );
      s.add( vl, 100 );
      return s;
   }

   function numeric( text, tip, key, lo, hi, prec, keep )
   {
      var n = new NumericControl( dlg );
      n.label.text = text;
      n.label.minWidth = labelW;
      n.setRange( lo, hi );
      n.slider.setRange( 0, 200 );
      n.slider.minWidth = 160;
      n.setPrecision( prec );
      n.setValue( P[key] );
      n.toolTip = tip;
      n.onValueUpdated = function( v )
      {
         P[key] = (prec == 0) ? Math.round( v ) : v;
         dlg.changed( keep );
      };
      dlg.setters[key] = function() { n.setValue( P[key] ); };
      return n;
   }

   function check( text, tip, key, keep )
   {
      var cb = new CheckBox( dlg );
      cb.text = text;
      cb.toolTip = tip;
      cb.checked = P[key];
      cb.onCheck = function( on ) { P[key] = on; if ( keep >= 0 ) dlg.changed( keep ); };
      dlg.setters[key] = function() { cb.checked = P[key]; };
      var s = new HorizontalSizer;
      s.addSpacing( labelW + 4 );
      s.add( cb );
      s.addStretch();
      return s;
   }

   function combo( text, tip, key, items, keep, keepCache )
   {
      var cb = new ComboBox( dlg );
      for ( var i = 0; i < items.length; ++i )
         cb.addItem( items[i] );
      cb.currentItem = Math.min( items.length - 1, Math.max( 0, Math.round( P[key] ) ) );
      cb.toolTip = tip;
      cb.onItemSelected = function( i )
      {
         P[key] = i;
         dlg.changed( keep, keepCache );
      };
      dlg.setters[key] = function() { cb.currentItem = Math.round( P[key] ); };
      var s = new HorizontalSizer;
      s.spacing = 4;
      s.add( label( text ) );
      s.add( cb, 100 );
      return s;
   }

   // --- preview pane ---------------------------------------------------------
   //   Shows the whole image at 1/3 resolution. Drag a box to render that
   //   area at full resolution (dragging again inside it zooms further);
   //   double-click or "Full image" returns to the whole image.

   var pv = new Control( this );
   this.pv = pv;
   pv.setMinSize( 640, 520 );
   pv.shown = null;                  // { bitmap, x0, y0, step }
   pv.scaled = null;
   pv.scale = 1;
   pv.ox = 0;
   pv.oy = 0;
   pv.dragging = false;
   pv.message = "Press Preview to analyze the selected images.";
   pv.toolTip = "Drag a box to view that area at full resolution. Double-click to return to the whole image.";

   pv.onPaint = function()
   {
      var g = new Graphics( this );
      g.fillRect( 0, 0, this.width, this.height, new Brush( 0xff202020 ) );
      if ( this.shown != null )
      {
         var bmp = this.shown.bitmap;
         var sc = Math.min( this.width/bmp.width, this.height/bmp.height );
         if ( this.scaled == null || sc != this.scale )
         {
            this.scale = sc;
            this.scaled = (sc == 1) ? bmp : bmp.scaled( sc );
         }
         this.ox = (this.width - this.scaled.width) >> 1;
         this.oy = (this.height - this.scaled.height) >> 1;
         g.drawBitmap( this.ox, this.oy, this.scaled );
         if ( this.dragging )
         {
            g.pen = new Pen( 0xffffd040, 1 );
            g.brush = new Brush( 0x30ffd040 );
            g.drawRect( Math.min( this.dragX, this.curX ), Math.min( this.dragY, this.curY ),
                        Math.max( this.dragX, this.curX ), Math.max( this.dragY, this.curY ) );
         }
      }
      else
      {
         g.pen = new Pen( 0xffc0c0c0 );
         g.drawText( 12, 24, this.message );
      }
      g.end();
   };

   pv.onResize = function()
   {
      this.update();
   };

   // display coordinates -> image coordinates
   pv.toImage = function( x, y )
   {
      var s = this.shown;
      return { x: s.x0 + (x - this.ox)/this.scale*s.step,
               y: s.y0 + (y - this.oy)/this.scale*s.step };
   };

   pv.onMousePress = function( x, y )
   {
      if ( this.shown != null && !dlg.busy )
      {
         this.dragging = true;
         this.dragX = this.curX = x;
         this.dragY = this.curY = y;
      }
   };

   pv.onMouseMove = function( x, y )
   {
      if ( this.dragging )
      {
         this.curX = x;
         this.curY = y;
         this.update();
      }
   };

   pv.onMouseRelease = function( x, y )
   {
      if ( !this.dragging )
         return;
      this.dragging = false;
      this.update();
      if ( Math.abs( x - this.dragX ) < 6 || Math.abs( y - this.dragY ) < 6 )
         return;
      var p0 = this.toImage( Math.min( x, this.dragX ), Math.min( y, this.dragY ) );
      var p1 = this.toImage( Math.max( x, this.dragX ), Math.max( y, this.dragY ) );
      var x0 = Math.max( 0, Math.floor( p0.x ) ), y0 = Math.max( 0, Math.floor( p0.y ) );
      var x1 = Math.min( engine.w, Math.ceil( p1.x ) ), y1 = Math.min( engine.h, Math.ceil( p1.y ) );
      if ( x1 - x0 < 16 || y1 - y0 < 16 )
         return;
      P.view = { x0: x0, y0: y0, x1: x1, y1: y1 };
      dlg.changed( 3 );
   };

   pv.onMouseDoubleClick = function()
   {
      dlg.fullImage();
   };

   this.fullImage = function()
   {
      if ( P.view == null )
         return;
      P.view = null;
      dlg.changed( 3 );
   };

   this.status = new Label( this );
   this.status.text = "";

   var previewBtn = new PushButton( this );
   previewBtn.text = "Preview";
   previewBtn.toolTip = "Analyze the selected images and show the result. Afterwards the preview follows the controls.";
   previewBtn.onClick = function()
   {
      engine.hasRun = true;
      dlg.refresh();
   };

   var showRow = combo( "Show:", "What the preview displays.", "previewMode",
                        [ "Result (restored starless)", "Starless before", "Protection mask", "Removed by the star tool" ], 3, true );
   var flip = new PushButton( this );
   flip.text = "Before / After";
   flip.toolTip = "Switch the preview between the starless image before StellarShield and the restored result.";
   flip.onClick = function()
   {
      P.previewMode = (P.previewMode == 0) ? 1 : 0;
      dlg.setters.previewMode();
      if ( engine.hasRun )
         dlg.refresh();
   };
   showRow.add( flip );

   var fullBtn = new PushButton( this );
   fullBtn.text = "Full image";
   fullBtn.toolTip = "Return to the whole image.";
   fullBtn.onClick = function() { dlg.fullImage(); };
   this.updateViewButton = function()
   {
      fullBtn.enabled = P.view != null;
   };
   showRow.add( fullBtn );

   var right = new VerticalSizer;
   right.spacing = 6;
   right.add( pv, 100 );
   right.add( showRow );
   right.add( this.status );

   // --- left column ----------------------------------------------------------

   var info = new Label( this );
   info.useRichText = true;
   info.wordWrapping = true;
   info.margin = 4;
   info.text = "<p><b>" + TITLE + " v" + VERSION + "</b> | Zachary Eichelberger</p>" +
               "<p>Restores non-stellar objects that star-removal tools take away, specifically galaxies, " +
               "back into the targeted starless image.</p>" +
               "<p>Press <i>Preview</i> to see results before applying.</p>";

   var active = ImageWindow.activeWindow;
   var preset = (active != null && !active.isNull) ? active.mainView : null;

   var reset = new PushButton( this );
   reset.text = "Reset";
   reset.toolTip = "Return every setting to its default. The selected images are kept.";
   reset.onClick = function()
   {
      var reload = P.blend != DEFAULTS.blend || P.gaiaRelease != DEFAULTS.gaiaRelease;
      for ( var k in PERSIST )
      {
         P[k] = DEFAULTS[k];
         if ( dlg.setters[k] )
            dlg.setters[k]();
      }
      dlg.changed( reload ? 0 : 1 );
   };
   var apply = new PushButton( this );
   apply.text = "Apply";
   apply.toolTip = "Restore on the full image and close.";
   apply.onClick = function() { dlg.ok(); };
   var cancel = new PushButton( this );
   cancel.text = "Close";
   cancel.onClick = function() { dlg.cancel(); };
   var buttons = new HorizontalSizer;
   buttons.spacing = 6;
   buttons.add( previewBtn );
   buttons.add( reset );
   buttons.addStretch();
   buttons.add( apply );
   buttons.add( cancel );

   var left = new VerticalSizer;
   left.spacing = 6;
   left.add( info );
   left.add( viewRow( "Starless image:", "Output of the star remover. Modified in place (undoable).", "starlessView", preset ) );
   left.add( viewRow( "Original image:", "Image before star removal (preferred).", "originalView", null ) );
   left.add( viewRow( "Stars image:", "Optional. Restored objects are removed from it. Required if no original is given.", "starsView", null ) );
   left.add( combo( "Gaia catalog:",
                    "Which Gaia XPSD database files to match against. Auto uses the deepest set configured " +
                    "(DR3, then EDR3, DR3/SP, DR2). DR3 is much deeper than DR3/SP, which only holds sources " +
                    "with BP/RP spectra (the set SPCC uses).",
                    "gaiaRelease", [ "Auto (deepest available)", "DR3", "DR3/SP" ], 0 ) );
   left.add( combo( "Stars image blend:",
                    "Only used when no original image is selected: how the stars image combines with the starless image. " +
                    "Auto detects it. Added: original = starless + stars (StarX default, StarNet). " +
                    "Screened: StarXTerminator with Unscreen. A wrong setting over-brightens restored galaxy cores.",
                    "blend", [ "Auto", "Added", "Screened (StarX Unscreen)" ], 0 ) );
   left.add( check( "Linear data (enable photometric tests)", "Uncheck for stretched images: only Gaia matching and morphology are used.", "linear", 2 ) );
   left.addSpacing( 4 );
   left.add( numeric( "Detection threshold:", "Residual threshold in units of background noise sigma. Lower finds fainter removed objects.", "kSigma", 0.5, 10, 2, 1 ) );
   left.add( numeric( "Minimum area (px):", "Removed objects smaller than this are ignored.", "minArea", 2, 50, 0, 1 ) );
   left.add( numeric( "Maximum area (px):", "Removed objects larger than this are left alone.", "maxArea", 100, 50000, 0, 2 ) );
   left.add( numeric( "Gaia limit (G mag):", "Gaia is assumed complete brighter than this: a brighter removed object with no Gaia match is restored.", "completeMag", 15, 21, 1, 2 ) );
   left.add( numeric( "Flux excess (mag):", "A removed object this much brighter than its matched Gaia sources is restored (galaxy core in Gaia, or star on a galaxy).", "excessMag", 0.2, 3, 2, 2 ) );
   left.add( combo( "Faint sources not in Gaia:",
                    "What to do with removed objects that are fainter than the Gaia limit and have no Gaia match. " +
                    "At high galactic latitude most of these are galaxies, so 'Protect all' is worth trying.",
                    "faintPolicy", [ "Protect if elongated or larger than stars", "Protect all", "Treat as stars" ], 2 ) );
   left.add( numeric( "Shape sensitivity:", "Sigmas above the stellar locus in size or ellipticity. Lower protects more.", "sizeSigma", 1, 8, 1, 2 ) );
   left.add( numeric( "Minimum ellipticity:", "Objects are never called elongated below this.", "minEll", 0.1, 0.8, 2, 2 ) );
   left.add( numeric( "Bright star guard (G mag):", "Unmatched fragments inside the halo or along the diffraction spikes of stars brighter than this are ignored. The zone follows each star's own shape.", "guardMag", 4, 16, 1, 2 ) );
   left.addSpacing( 4 );
   left.add( numeric( "Mask growth (px):", "How far the protection mask extends beyond each restored object (round).", "dilate", 0, 15, 0, 3 ) );
   left.add( numeric( "Mask feather (px):", "Width of the soft edge of the mask.", "feather", 0, 50, 0, 3 ) );
   left.add( check( "Mask only (do not modify images)", "Apply creates StellarShield_mask and the console report, without changing the images.", "maskOnly", -1 ) );
   left.addStretch();
   left.add( buttons );

   this.sizer = new HorizontalSizer;
   this.sizer.margin = 8;
   this.sizer.spacing = 8;
   this.sizer.add( left );
   this.sizer.add( right, 100 );

   this.updateViewButton();
   this.windowTitle = TITLE + " " + VERSION;
   this.adjustToContents();
}

GGDialog.prototype = new Dialog;

function main()
{
   console.show();
   loadSettings();
   P.view = null;
   var engine = new Engine;
   var dialog = new GGDialog( engine );
   var ok = dialog.execute();
   saveSettings();
   if ( !ok )
      return;
   try
   {
      var t0 = (new Date).getTime();
      engine.ensure( function( s ) { console.writeln( s ); console.flush(); } );
      engine.apply();
      console.writeln( format( "Elapsed: %.1f s", ((new Date).getTime() - t0)/1000 ) );
   }
   catch ( e )
   {
      console.criticalln( "*** " + TITLE + ": " + e.message );
      (new MessageBox( e.message, TITLE, StdIcon_Error, StdButton_Ok )).execute();
   }
}

main();
