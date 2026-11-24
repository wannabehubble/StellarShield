# StellarShield

A PixInsight script that restores non-stellar objects that star-removal tools (StarXTerminator, StarNet) take away, specifically galaxies, back into the starless image.

## How it works

1. The light the star remover took is the original image minus the starless image (or the stars image itself).
2. That light is split into separate objects.
3. Each object is matched against the Gaia catalogue and measured (brightness, size, elongation).
4. Objects judged non-stellar are copied back into the starless image through a soft-edged mask, and the stars image is corrected to match.

## Requirements

- PixInsight with the Gaia process configured with XPSD database files. The full DR3 set is strongly preferred; DR3/SP is much shallower.
- A plate-solved image (any one of the selected images).

## Installation

Download `StellarShield.js`, then in PixInsight choose **Script > Feature Scripts... > Add**, point it at the folder containing the file and click **Done**. The script appears under **Script > Utilities > StellarShield**. You can also run it once with **Script > Execute Script File...**.

## Usage

1. Run the star remover. Either work on a clone and keep the original, or generate a stars image.
2. Open StellarShield and select the starless image plus the original (preferred) and/or the stars image.
3. Press **Preview**. The preview shows the whole image at 1/3 resolution; drag a box to view an area at full resolution, and double-click or press **Full image** to go back. **Before / After** toggles between the starless image and the result.
4. Adjust the settings and press **Apply**. Changes are undoable (Ctrl+Z). A `StellarShield_mask` image shows what was restored.

## Settings

- **Gaia catalog** – which Gaia database files to match against: Auto (deepest available), DR3 or DR3/SP.
- **Stars image blend** – how the stars image combines with the starless image when no original is selected. Auto detects it; choose Screened for StarXTerminator's Unscreen Stars option.
- **Linear data** – enables the brightness tests; uncheck for stretched images.
- **Detection threshold / Minimum area / Maximum area** – which removed objects are considered.
- **Gaia limit / Flux excess / Faint sources not in Gaia / Shape sensitivity / Minimum ellipticity** – how objects are judged non-stellar.
- **Bright star guard** – fragments in the halo or spikes of stars brighter than this are ignored.
- **Mask growth / Mask feather** – size and soft edge of the restoration mask.
- **Mask only** – create the mask and console report without changing the images.

Settings are remembered between runs; **Reset** returns them to their defaults.

## Author

Zachary Eichelberger
