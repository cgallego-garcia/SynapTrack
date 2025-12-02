# SynapTrack

SynapTrack is an automated Fiji/ImageJ pipeline for quantifying synaptic contacts in neuronal cultures and tissue sections. 
Built on the SynQuant segmentation engine, it performs fully unsupervised detection and quantification of synaptic puncta from fluorescence image sets that include nuclei (DAPI), dendrites (MAP2), and pre- and postsynaptic markers for excitatory or inhibitory synapses.

---

## Main Features

* Fully automated segmentation and synapse detection
* Compatible with both cultured neurons and tissue sections
* Based on SynQuant with multiple channels
* Export of synaptic density per cell and per dendrite length and dendrite measurements

---

## Requirements

* Fiji/ImageJ (recommended ≥ 2.16/1.53f) (Download [here](https://imagej.net/software/fiji/downloads))
* Java ≥ 8
* Required plugins:
	* Bio-Formats Importer
	* StarDist 2D
	* CSBDeep / TensorFlow
	* IJPB plugins (MorphoLibJ)
	* Neuronanatomy
	* [SynQuant 1.2.8](https://github.com/yu-lab-vt/SynQuant)

---

## Installation and Basic Use

1. Download SynapTrack and the associated macros from this repository
2. Copy the downloaded SynapTrack folder in `Fiji.app/plugins`
3. Restart FIJI

4. Run SynapTrack. It could be done in two-ways:

    (A) Drag and drop the `SynapTrack.ijm` file and press `Run`

    (B) SynapTrack will appear as an option under `Plugins > SynapTrack`

5. Adjust parameters.
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

#### Image Calibration

(If images are not calibrated)
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
