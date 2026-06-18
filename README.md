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
  - Download Fiji [here](https://imagej.net/software/fiji/downloads)

- **Required Plugins** (install most plugins via `Fiji's Help > Update...` by activating the corresponding update site):

| Plugin | Version |
|--------|---------|
| Bio-Formats | ≥ 8.4.0 |
| StarDist 2D | ≥ 0.7.0 |
| CSBDeep | ≥ 2.0 |
| TensorFlow | ≥ 0.6.0 |
| IJPB plugins (MorphoLibJ) | ≥ 7.0.0 |
| Neuronanatomy | ≥ 1.6.0 |
| [SynQuant](https://github.com/yu-lab-vt/SynQuant#getting-started) | 1.2.8 |
> **Note:** SynQuant requires manual installation as follows:  
> 1. Download the plugin from [SynQuant GitHub Repository](https://github.com/yu-lab-vt/SynQuant#getting-started).
> 2. Copy it into: `Fiji.app/plugins`

---

### **Basic Use**
1. Download **SynapTrack** and the associated macros from this repository using the **<> Code → Download ZIP** option at the top of this repository.

<img width="380" height="307" alt="imagen" src="https://github.com/user-attachments/assets/446bbf35-9fd5-434b-8347-9dd6b102beca" />

2. **Extract and Copy the `SynapTrack` folder (SynapTrack-main)** into the Fiji plugin directory: `Fiji.app/plugins`

3. Restart Fiji to ensure the plugin is recognized.

4. Execute SynapTrack. It could be done in two-ways:

    (A) Drag and drop the `SynapTrack.ijm` file and press `Run`

    (B) SynapTrack will appear as an option under `Plugins > SynapTrack`

5. SynapTrack can be tested using the sample images provided in `SynapTrack-Main/test_images`.

> **IMPORTANT:** Before using the macro, please ensure that your images are named correctly. Refer to → [Image Preparation](https://github.com/cgallego-garcia/SynapTrack#Image-Preparation)

6. Specify the parameters for your analysis (e.g., dendrite channel, synaptic channels, thresholds) and run the analysis.  
   - *Adjust parameters (*see [Interactive Parameters panel](https://github.com/cgallego-garcia/SynapTrack#Interactive-Parameters)*) as needed.*

7. Run the analysis to generate synaptic density and associated metrics.

8. Once the analysis is complete, a `Results` folder will be created within the image directory, following the structure described in → [Outputs](https://github.com/cgallego-garcia/SynapTrack#Output-structure).

A [workflow example](https://github.com/cgallego-garcia/SynapTrack#Example-Workflow) could be found at the bottom of this repository

---

### Image Preparation

SynapTrack expects files to follow this naming convention: `<ExpPrefix>_<Channel>_<Index>.tif`

Where:

   - `<ExpPrefix>` - Experiment identifier (e.g., ExpInhib_01 | ExpExcit_01)
   - `<Channel>` - Specifies the marker one of the following (preEx, preInh, postEx, postInh, MAP2, or DAPI)
   - `<Index>` - two-digit replicate identifier (e.g., 01, 02)

#### - Example of Inhibitory images:

`ExpInhib_01_DAPI_01.tif | ExpInhib_01_MAP2_01.tif | ExpInhib_01_preInh_01.tif | ExpInhib_01_postInh_01.tif`

#### - Example of Excitatory images:

`ExpExcit_01_DAPI_01.tif | ExpExcit_01_MAP2_01.tif | ExpExcit_01_preEx_01.tif | ExpExcit_01_postEx_01.tif`

If your data are in proprietary formats (.czi, .vsi, .nd2, .lif, etc.), you can run `SynapTrack_FileConversion.ijm` to generate SynapTrack-compatible TIFF files.

---

## Interactive Parameters

The following panel is used to configure SynapTrack before running the analysis.

<img width="707" height="800" alt="imagen" src="https://github.com/user-attachments/assets/50f9b064-cc68-4d5e-a3af-bb578a70a01f" />

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

## Example Workflow

This example demonstrates the analysis of inhibitory synapses in cultured neurons using the sample dataset located in:

`SynapTrack-main/test_images/neuronal_cultures/inhibitory`

### Input Images

For each field of view, SynapTrack requires four images:

| Channel | Marker | Function |
|----------|----------|----------|
| DAPI | Cell nuclei | Cell counting |
| MAP2 | Dendrites | Dendrite segmentation and length measurements |
| preInh | Presynaptic marker (e.g., VGAT) | Presynaptic puncta detection |
| postInh | Postsynaptic marker (e.g., Gephyrin) | Postsynaptic puncta detection |

Example image set:

`ExpInhib_01_DAPI_01.tif | ExpInhib_01_MAP2_01.tif | ExpInhib_01_preInh_01.tif | ExpInhib_01_postInh_01.tif`

### Running the Analysis

1. Open SynapTrack in Fiji.
2. Select the folder containing the example images.
3. Set:
   - Experiment prefix: `ExpInhib_01`
   - First image index: `01`
   - Last image index: `03` (or last available image)
   - Synapse type: `Inhibitory`
4. Use default parameters for the first test run.
5. Click **Run**.

### Output structure

Upon completion of the analysis, SynapTrack generates a `Results` folder inside the analyzed image folder containing:

### Results File

The final output file `Results.xls` is a summary table for all analyzed image sets containing:

- Number of detected cells
- Number of synaptic contacts
- Synapses per cell
- Total dendrite length (µm)
- Synapses per 10 µm dendrite length

Each analyzed image set generates an `Image_/` subfolder containing:

- Dendrite mask - `Dendrite_enlarged.tif`
- Dendrite skeleton - `Dendrite_skeleton.tif`
- Dendrite measurements - `Dendrite_stats.csv`
- Nuclei segmentation (ROIs) -`ROIset_nuclei.zip`
- Synapse ROIs - `ROIset_Synapse.zip`
- Synapse detection overlay images - `SynapseDetection.tif`

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

---

## Changelog
### v1.0.1 (2026-01-20)
- Added defensive checks for missing Results and Summary tables
- Improved error messages and safe termination

### v1.1 (2026-04-05)
- Added Plugin availability check function in SynapTrack macros

