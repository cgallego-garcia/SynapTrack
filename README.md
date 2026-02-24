# SynapTrack

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17795689.svg)](https://doi.org/10.5281/zenodo.17795689)

SynapTrack is an automated Fiji/ImageJ pipeline for quantifying synaptic contacts in neuronal cultures and tissue sections. 
Built on the SynQuant segmentation engine, it performs fully unsupervised detection and quantification of synaptic puncta from fluorescence image sets that include nuclei (DAPI), dendrites (MAP2), and pre- and postsynaptic markers for excitatory or inhibitory synapses.

<img width="800" height="800" alt="imagen" src="https://github.com/user-attachments/assets/626d0e91-a11e-4d69-984e-ecf7f3c9dcb3" />

---

## Main Features

* Fully automated segmentation and synapse detection
* Compatible with both cultured neurons and tissue sections
* Based on SynQuant with multiple channels
* Export of synaptic density per cell and per dendrite length and dendrite measurements

---

## Getting Started
### Requirements
**Before using SynapTrack, make sure you have the following installed:**
- **Fiji / ImageJ**  
  - Recommended versions: Fiji ≥ 2.16 / ImageJ ≥ 1.53f  
  - [Download Fiji here]([https://imagej.net/software/fiji/downloads](https://imagej.net/software/fiji/downloads))
    
- **Java**  
  - Version ≥ 8

- **Required Plugins** (install most plugins via `Fiji's Help > Update...` by activating the corresponding update site):

| Plugin | Version |
|--------|---------|
| Bio-Formats | ≥ 8.4.0 |
| StarDist 2D | ≥ 0.7.0 |
| CSBDeep | ≥ 2.0 |
| TensorFlow | ≥ 0.6.0 |
| IJPB plugins (MorphoLibJ) | ≥ 7.0.0 |
| Neuronanatomy | ≥ 1.6.0 |
| SynQuant | 1.2.8 |
> **Note:** SynQuant requires manual installation as follows:  
> 1. Download the plugin from [SynQuant GitHub Repository](https://github.com/yu-lab-vt/SynQuant).
> 2. Copy it into: `Fiji.app/plugins`

---

### **Basic Use**
1. **Download SynapTrack** and the associated macros from this repository.
2. **Copy the `SynapTrack` folder** into the Fiji plugin directory: `Fiji.app/plugins`
3. Restart Fiji to ensure the plugin is recognized.

4. Execute SynapTrack. It could be done in two-ways:

    (A) Drag and drop the `SynapTrack.ijm` file and press `Run`

    (B) SynapTrack will appear as an option under `Plugins > SynapTrack`

5. Specify the parameters for your analysis (e.g., dendrite channel, synaptic channels, thresholds) and run the analysis.  
   - *Adjust parameters (*see Interactive Parameters panel below*) as needed.*
6. Run the analysis to generate synaptic density and associated metrics.

---

### Image Preparation

SynapTrack expects files to follow this naming convention: `<ExpPrefix>_<Channel>_<Index>.tif`

Where:

   - `<ExpPrefix>` - Experiment identifier (e.g., Exp01)
   - `<Channel>` - Specifies the marker one of the following (preEx, preInh, postEx, postInh, MAP2, or DAPI)
   - `<Index>` - two-digit replicate identifier (e.g., 01, 02)

If your data are in proprietary formats (.czi, .vsi, .nd2, .lif, etc.), you can run `SynapTrack_FileConversion.ijm` to generate SynapTrack-compatible TIFF files.

---

### Interactive Parameters

<img width="612" height="626" alt="imagen" src="https://github.com/user-attachments/assets/06ac4adb-8c6c-4348-92cc-a7c7ffef3d9c" />

#### Input Data

-   Input images folder
-   Experiment prefix
-   Image index range (First-Last image sets to be analysed)
-   Synapse type (excitatory / inhibitory)

#### Image Calibration *(If images are not calibrated)*
-   Pixel size (µm/pixel)
-   Camera binning (if metadata is missing)

#### Image Preprocessing

-   Dendrite-specific CLAHE and Channel-specific background subtraction
-   Background subtraction rolling ball radii defined in µm
-   Z-stack processing: Max Intensity or Sum of Slices

#### Nuclei Segmentation

-   StarDist-based detection in nuclei (DAPI) images
-   Size (in µm) and circularity filtering
-   Automated cell counting

#### SynQuant

-   z-score threshold
-   Size (µm²) and roundness filters for puncta

#### Dendrite Mask

-   Expansion of MAP2 mask by a user-defined distance
-   Ensures detection of puncta near dendrites

---

### Outputs

SynapTrack generates a Results folder containing:

Top-level

-   Results.xls — global summary of all image sets

-   Per-image subfolder (Image_/)
      - Dendrite_enlarged.tif
      - Dendrite_skeleton.tif
      - Dendrite_stats.csv
      - ROIset_nuclei.zip
      - ROIset_Synapse.zip
      - SynapseDetection.tif

Example:

    Results/
        Results.xls
        1/
            Dendrite_enlarged.tif
            Dendrite_skeleton.tif
            Dendrite_stats.csv
            ROIset_nuclei.zip
            ROIset_Synapse.zip
            SynapseDetection.tif


## Changelog
### v1.0.1 (2026-01-20)
- Added defensive checks for missing Results and Summary tables
- Improved error messages and safe termination


