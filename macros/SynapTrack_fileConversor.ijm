/* ============================================================================
 *  IMAGE FORMAT PREPARATION MACRO
 *  ---------------------------------------------------------------------------
 *  This macro is intended for users whose original data come in *non–TIFF*
 *  proprietary formats (Zeiss .czi / .lsm, Nikon .nd2, Leica .lif).
 *  It performs:
 *  1) Workspace cleanup
 *  2) User input dialog (folder, microscope system, channel assignment)
 *  3) Batch import using Bio-Formats
 *  4) Channel splitting
 *  5) Systematic renaming into a standardized 4-channel TIFF structure
 *  
 *  Output:
 *  	/splitTiffs/
 *  		Experiment_DAPI_01.tif
 *  		Experiment_MAP2_01.tif
 *  		Experiment_preEx_01.tif   or preInh
 *  		Experiment_postEx_01.tif  or postInh
 *  		
 *  This macro is *required* to normalize datasets so they can be processed by
 *  SynapTrack v1 regardless of the microscope manufacturer.
============================================================================ */
	
	initWorkspace(); // Clean Fiji environment

/* ============================ PARAMETERS DIALOGUE ======================================
   GUI for user to supply experiment paths, system and channel mapping.
   ======================================================================================= */
	channelArray = newArray("1","2","3","4");
	excitatory = newArray("preEx","postEx");
	inhibitory = newArray("preInh","postInh");

	Dialog.create("Images_Folder");
	Dialog.addDirectory("Input Folder", File.getDefaultDir);
	Dialog.addChoice("Microscope system", newArray("Zeiss-Nikon-Evident","Leica"));
	Dialog.addString("Experiment Prefix", "NB_03"); Dialog.addToSameRow();
	Dialog.addChoice("Synapse Type", newArray("Excitatory", "Inhibitory"), "Excitatory");
	// Channel mapping
	Dialog.addMessage("Channel order selection");
	Dialog.addChoice("Nuclei", channelArray, channelArray[0]); Dialog.addToSameRow();
	Dialog.addChoice("Dendrites", channelArray, channelArray[1]); Dialog.addToSameRow();
	Dialog.addChoice("PreSyn", channelArray, channelArray[2]); Dialog.addToSameRow();
	Dialog.addChoice("PostSyn", channelArray, channelArray[3]);
    Dialog.addNumber("If z-stack, How many Zs would you like per image?:", 3);
	Dialog.addCheckbox("Batch Mode?", true);
	Dialog.show();
	
	inputFolder 		= Dialog.getString();
	outputFolder 		= inputFolder+File.separator+"splitTiffs"+File.separator; File.makeDirectory(outputFolder);
	imgFormat 			= Dialog.getChoice();
	experimentPrefix 	= Dialog.getString();
	synapseType 		= Dialog.getChoice(); 
	nucleiCh 			= Dialog.getChoice(); 
	dendritesCh 		= Dialog.getChoice();
	preSynCh 			= Dialog.getChoice();
	postSynCh 			= Dialog.getChoice();
	slicesPerFragment	= Dialog.getNumber();
	batchMode 			= Dialog.getCheckbox();
	digitN				= 1;
	setBatchMode(batchMode);
	
/* ============================ PROCESSING LOOP ======================================
   Main loop: iterate over image sets.
   ================================================================================ */

	imgList = getFileList(inputFolder);
	
	for (i = 0; i < imgList.length; i++){
		// Prevent the user from selecting duplicated channels
		if (nucleiCh != dendritesCh && nucleiCh != preSynCh && nucleiCh != postSynCh && dendritesCh != preSynCh && dendritesCh != postSynCh && preSynCh != postSynCh) {
			// Select correct import function based on system
			if (imgFormat == "Zeiss-Nikon-Evident") {
				if (endsWith(imgList[i], ".czi") || endsWith(imgList[i], ".lsm") || endsWith(imgList[i], ".nd2")) { imgName = zeissNikonOpenning(inputFolder, imgList); }
			}
			
			else if (imgFormat == "Leica") {
				if (endsWith(imgList[i], ".lif")) { imgName = leicaOpenning(inputFolder, imgList); }
			}
			
			if (nImages > 0 && nImages < 2) {
				getDimensions(width, height, channels, slices, frames);
				// Basic safety check
				if (channels < 4) {
					// Skip this file
					print("[WARN] Dataset "+imgList[i]+ " does not contains the minimum channels");
					run("Close All");
				}
				else {
					digitN = imgProcessing(imgList, imgName, experimentPrefix, digitN, excitatory, inhibitory, synapseType, slicesPerFragment, nucleiCh, dendritesCh, preSynCh, postSynCh, outputFolder); 
				}
			}
		}
		
		else { exit("There was an error selecting the channels and some of them are repeated"); }
	}
	
/* ============================ END MESSAGE ======================================
   ================================================================================ */
	Dialog.create("Processing Ended");
	Dialog.addMessage("Image Preparation has finished correctly",14);
	Dialog.show();
	
/* ============================ FUNCTIONS ===========================================
   Utility functions used across the macro.
   ================================================================================ */
	function initWorkspace(){
	// Close images and helper windows and force garbage collection.
		if (nImages > 0) { run("Close All"); } // Close all open images
		close("ROI Manager"); // Clean & close ROIManager
		close("Log"); close("Results"); // close log, results & Data  windows
		run("Collect Garbage");	// Free up unused ram memory
		run("Bio-Formats Macro Extensions"); // Activate Bio-Formats Macro Extensions
	}

	function zeissNikonOpenning(inputFolder, imgList) {
	// Generic importer for Zeiss, Nikon and Olympus formats
		run("Bio-Formats Importer", "open=["+inputFolder+imgList[i]+"] autoscale color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		imgName = File.nameWithoutExtension;
		rename(imgName); // Rename window for easy handling
		return imgName;
	}
	
	function leicaOpenning(inputFolder, imgList) {
	// LIF files may contain multiple series; this iterates through series.
	// Returns the series name of the first series imported.
		Ext.setId(inputPath+imgList[i]);
		Ext.getSeriesCount(seriesCount);
		
		for (s = 1; s < seriesCount+1; s++) {
			Ext.setSeries(s-1);
			Ext.getSeriesName(seriesName);
			
			run("Bio-Formats Importer", "open=["+inputFolder+imgList[i]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s);
			rename(seriesName);
			return seriesName;
		}
	}
	
	function imgProcessing(imgList, imgName, experimentPrefix, digitN, excitatory, inhibitory, synapseType, slicesPerFragment, nucleiCh, dendritesCh, preSynCh, postSynCh, outputFolder) {
	// Splits channels and saves them as standardized TIFFs
		getDimensions(width, height, channels, slices, frames);
		
		if (slices <= 1 || slicesPerFragment == slices) {
			splitAndSave(imgList, imgName, experimentPrefix, digitN, excitatory, inhibitory, synapseType, slicesPerFragment, nucleiCh, dendritesCh, preSynCh, postSynCh, outputFolder);
			digitN++; run("Close All");
		}
		
		else {
			rawID = getImageID();
			nFragments = floor(slices / slicesPerFragment);
		 	slicesValid = slicesPerFragment * nFragments;   // Slices que usaremos realmente
		 	slicesDiscarded = slices - slicesValid;         // Los que se descartan del final
		 	
		 	start = 1;
		 	
		 	for (f=1; f <= nFragments; f++) {
		 		end = start + slicesPerFragment - 1;
		 		tmpName = "tmpFrag";
		 		run("Duplicate...", "title="+tmpName+" duplicate slices=" + start + "-" + end);
		 		splitAndSave(imgList, tmpName, experimentPrefix, digitN, excitatory, inhibitory, synapseType, slicesPerFragment, nucleiCh, dendritesCh, preSynCh, postSynCh, outputFolder);
		 		selectImage(rawID); close("\\Others");
		 		digitN++;
		 		start = end + 1;
		 	}
		}
		run("Close All");
		return digitN;
	}
	
	function splitAndSave(imgList, imgName, experimentPrefix, digitN, excitatory, inhibitory, synapseType, slicesPerFragment, nucleiCh, dendritesCh, preSynCh, postSynCh, outputFolder) {
		run("Split Channels"); // Produces C1-, C2-, C3-, C4- windows
		d = d2(digitN); // Convert index → 2-digit string
	// Save Nuclei (DAPI)
		selectImage("C"+nucleiCh+"-"+imgName);
		saveAs("Tiff", outputFolder + experimentPrefix + "_DAPI_"+d+".tif");
	// Save Dendrites (MAP2)
		selectImage("C"+dendritesCh+"-"+imgName);
		saveAs("Tiff", outputFolder + experimentPrefix + "_MAP2_"+d+".tif");
	// Save Pre-Synaptic
		selectImage("C"+preSynCh+"-"+imgName);
		if (synapseType == "Excitatory") { saveAs("Tiff", outputFolder + experimentPrefix + "_"+excitatory[0]+"_"+d+".tif"); }
		else if (synapseType == "Inhibitory") { saveAs("Tiff", outputFolder + experimentPrefix + "_"+inhibitory[0]+"_"+d+".tif"); }
	// Save Post-Synaptic
		selectImage("C"+postSynCh+"-"+imgName);
		if (synapseType == "Excitatory") { saveAs("Tiff", outputFolder + experimentPrefix + "_"+excitatory[1]+"_"+d+".tif"); }
		else if (synapseType == "Inhibitory") { saveAs("Tiff", outputFolder + experimentPrefix + "_"+inhibitory[1]+"_"+d+".tif"); }
	}
	
	function d2(n){
	// Utility: converts integer → 2-digit string (1 → "01")
		if (n < 10) { return "0"+n; }
		return ""+n;
	}

	
	
	
	
	
	


