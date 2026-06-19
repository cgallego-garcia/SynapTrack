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
* Export of synapse number per cell and per dendrite length and dendrite measurements

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
> 1. Download the file `SynQuantVid_-1.2.8.jar` from [SynQuant GitHub Repository](https://github.com/yu-lab-vt/SynQuant/releases).
> 2. Copy it into: `Fiji.app/plugins`

---

### **Basic Use**
1. Download **SynapTrack-main.zip** by clicking the **<> Code → Download ZIP** option at the top of this repository.

<img width="380" height="307" alt="imagen" src="https://github.com/user-attachments/assets/446bbf35-9fd5-434b-8347-9dd6b102beca" />

2. **Extract and Copy `SynapTrack-main` folder** into the Fiji plugin directory: `Fiji.app/plugins`

3. Restart Fiji to ensure the plugin is loaded.

4. Execute SynapTrack. It could be done in two-ways:

    (A) Drag and drop the `SynapTrack.ijm` file within main FIJI window and press `Run`

    (B) Click SynapTrack option under `Plugins > SynapTrack`

5. SynapTrack example analysis can be done using the sample images provided in `SynapTrack-Main/test_images`. 

    A [workflow example](https://github.com/cgallego-garcia/SynapTrack#Example-Workflow) could be found at the bottom of this repository

> **IMPORTANT:** Before running the macro, please ensure that the images file names have the correct syntax. Refer to → [Image Preparation](https://github.com/cgallego-garcia/SynapTrack#Image-Preparation)

6. Run SynapTrack and specify the parameters for your analysis in the [StartUp panel](https://github.com/cgallego-garcia/SynapTrack#StartUp-Panel) (e.g., dendrite channel, synaptic channels, thresholds).

7. Run SynapTrack analysis.

8. Once the analysis is complete, a `Results` folder will be created within the image folder, following the [Output structure](https://github.com/cgallego-garcia/SynapTrack#Output-structure).

---

### Image Preparation

SynapTrack requires image filenames with the following syntax: `<ExpPrefix>_<Channel>_<Index>.tif`

Where:

   - `<ExpPrefix>` - User-defined experiment identifier (e.g., ExpInhib_01 | ExpExcit_01).
   - `<Channel>` - Must be one of the following markers: preEx, preInh, postEx, postInh, MAP2, or DAPI.
   - `<Index>` - Must be a two-digit replicate identifier starting at 01.

#### - Example of Inhibitory set of images:

    `ExpInhib_01_DAPI_01.tif | ExpInhib_01_MAP2_01.tif | ExpInhib_01_preInh_01.tif | ExpInhib_01_postInh_01.tif`

#### - Example of Excitatory set of images:

    `ExpExcit_01_DAPI_01.tif | ExpExcit_01_MAP2_01.tif | ExpExcit_01_preEx_01.tif | ExpExcit_01_postEx_01.tif`

Use `SynapTrack_FileConversion.ijm` to generate SynapTrack-compatible TIFF files, if your data are in proprietary formats (.czi, .vsi, .nd2, .lif, etc.).

---

## StartUp Panel

<img width="806" height="925" alt="image" src="https://github.com/user-attachments/assets/e8e4c4bc-06d3-4191-9079-ed684a45ec0f" />


#### Input Data

-   Input images folder: Enter the folder of your image files.
-   Experiment prefix: User-defined experiment identifier
-   Image index range (First-Last image sets to be analysed)
-   Synapse type: Either Excitatory or Inhibitory

#### Image Calibration *(Only if images are not calibrated)*
-   Pixel size (µm/pixel)
-   Camera binning (if metadata is missing)

#### Image Preprocessing

-   Background subtraction: Rolling ball radii (defined in µm for Nuclei, Synapses and Dendrites)
-   Enhanced Local Contrast Dendrites: Block size, Histogram bins and Maximum Slope
-   Dendrites Gaussian filter: Sigma filter (µm)
-   Z-stack processing. Z-handling: Max Intensity (MaxIP) or Sum of Slices (SumSlices)

#### Nuclei Segmentation

-   StarDist-based detection in nuclei (DAPI) images
-   Nuceli Area (in µm²) and Nuclei circularity

#### SynQuant

-   z-score threshold and zScore adjustment for probabilistic puncta selection
-   Min and Max Synapse size (µm²) and Max roundness (WH) ratio for puncta

#### Dendrite Mask

-   Enlarge by (µm): Expands the MAP2 mask to ensure detection of puncta close to dendrites

---

## Example Workflow

This example shows the analysis of inhibitory synapses in cultured neurons from `ExpInhib_01` dataset located in:

`SynapTrack-main/test_images/neuronal_cultures/inhibitory`

### Input Images

For each replicate, SynapTrack requires four images:

| Channel | Marker | Function |
|----------|----------|----------|
| DAPI | Cell nuclei | Cell counting |
| MAP2 | Dendrites | Dendrite segmentation and length measurement |
| preInh | Presynaptic marker (e.g., VGAT) | Presynaptic puncta detection |
| postInh | Postsynaptic marker (e.g., Gephyrin) | Postsynaptic puncta detection |

-Example image replicate:

    `ExpInhib_01_DAPI_01.tif | ExpInhib_01_MAP2_01.tif | ExpInhib_01_preInh_01.tif | ExpInhib_01_postInh_01.tif`

### Running the Analysis

1. Run SynapTrack in Fiji.
2. Input the folder containing the example images.
3. Set:
   - Experiment prefix: `ExpInhib_01`
   - First image index: `01`
   - Last image index: `05`
   - Synapse type: `Inhibitory`
4. Use the default parameters for the first test run.
5. Click **Run**.

### Output structure

Upon analysis completion, SynapTrack generates a `Results` folder inside the image folder containing:

    #### Results File
    
    The final output file `Results.xls` is a summary table for all image replicates containing:
    
        - Number of detected cells
        - Number of synaptic contacts
        - Synapses per cell
        - Total dendrite length (µm)
        - Synapses per 10 µm dendrite length
    
    Each analyzed image relpicate generates a subfolder containing:
    
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

