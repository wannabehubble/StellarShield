# StellarShield

By **Zachary Eichelberger**

Restores non-stellar objects that star-removal tools take away, specifically galaxies, back into the starless image.

Star removers such as StarXTerminator and StarNet often take small galaxies along with the stars. StellarShield finds what was removed, splits it into separate objects and checks each one against the Gaia catalog. It then measures each object's brightness, size and elongation. Objects judged non-stellar are copied back into the starless image through a soft-edged mask, and the stars image is corrected to match.

<img width="1583" height="1200" alt="Starless image with galaxies restored by StellarShield" src="images/surf_galax_resto.jpg" />

## Download

StellarShield is a PixInsight script.

| You use | Download | Version |
|---|---|---|
| **PixInsight** | [`StellarShield.js`](StellarShield.js) | 1.3.2 |

## Features

- **Gaia cross-matching** – removed objects are matched against Gaia DR3 or DR3/SP. Objects with no Gaia star, or with too much light or too large a shape for their Gaia star, are judged non-stellar.
- **Round, feathered masks** – restoration masks follow each object's shape, with an adjustable soft edge.
- **Bright star guard** – halos and diffraction spikes of bright stars are kept out of the restoration.
- **Works with original or stars image** – use the original image (preferred) or the stars-only image. Screened stars images (StarXTerminator's Unscreen Stars) are detected automatically.
- **Live preview** – the whole image at 1/3 resolution. Drag a box to view an area at full resolution, with a before/after toggle.
- **Remembered settings** – settings and image choices are kept between runs, with a **Reset** button.

Works on **linear and stretched (non-linear)** images.

## Requirements

- PixInsight with the Gaia process configured with XPSD database files. The full DR3 set is strongly preferred; DR3/SP is much shallower.
- A plate-solved image (any one of the selected images).

## Installation

Download `StellarShield.js`, then in PixInsight:

1. Choose **Script > Feature Scripts... > Add**.
2. Point it at the folder containing the file and click **Done**.

The script appears under **Script > Utilities > StellarShield**. You can also run it once with **Script > Execute Script File...**.

## Usage

1. Run the star remover. Either work on a clone and keep the original, or generate a stars image.
2. Open StellarShield and select the starless image plus the original (preferred) and/or the stars image.
3. Press **Preview**.
   - The preview shows the whole image at 1/3 resolution.
   - Drag a box to view an area at full resolution; double-click or press **Full image** to go back.
   - **Before / After** toggles between the starless image and the result.
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

## Contributing

Bug reports and suggestions are welcome. Please open an issue with:
- your PixInsight and StellarShield version numbers
- the settings you used
- a screenshot, if possible

## License

Released under the [MIT License](LICENSE). © 2026 Zachary Eichelberger.
