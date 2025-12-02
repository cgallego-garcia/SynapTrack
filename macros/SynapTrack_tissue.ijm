
/*========================================================================================
* SynapTrack - Tissue version 1.0
* 
* Description  : Robust, faster synapse quantification for 2D and z-stacks
*
* Authors      : Carlos Gallego-Garcia, Elena Martinez-Blanco & Fco. Javier Diez-Guerra
* Facility     : Advanced Light Microscopy Facility (SMOA), CBM
* Created      : 2025-06-28
*
* Requirements : Fiji (with Bio-Formats), StarDist 2D, CSBDeep, TensorFlow,
*                IJPB plugins (MorphoLibJ), Neuronanatomy, SynQuant 1.2.8
*========================================================================================*/


/* ============================ GLOBALS ======================================
 *  Minimal initialization and global constants. This section prepares Fiji
 *  and sets channel name patterns used across the macro.
 =========================================================================== */
	requires("1.53f");		// Ensure the required ImageJ/Fiji version
	row = initWorkspace();	// Clean/initialize the workspace (defined below) and set Row index for summary table accumulation
	
// Channel name patterns used to build file names and logical mapping
	excitatory = newArray("preEx","postEx");
	inhibitory = newArray("preInh","postInh");
	logical   = newArray("preSyn","postSyn");
	
/* ============================ PARAMETERS DIALOGUE ======================================
 *  GUI for user to supply experiment paths and algorithm parameters.
 ======================================================================================= */
	Dialog.create("SynapTrack Tissue: Analysis Parameters");

// ====== EXPERIMENT DATA ======
	Dialog.addMessage("Input Data", 16);
	Dialog.addDirectory("Input images folder:", File.getDefaultDir);
	Dialog.addString("Experiment Prefix:", "TissuEX_"); Dialog.addToSameRow();
	Dialog.addNumber("First ImageSet to analyze:", 1); Dialog.addToSameRow();
	Dialog.addNumber("Last ImageSet to analyze:", "");
	Dialog.addChoice("Synapse type:", newArray("Excitatory","Inhibitory"), "Excitatory");
	
// ====== IMAGE CALIBRATION ======
	Dialog.addMessage("Image Calibration", 16);
	Dialog.addChoice("Are the images Calibrated?", newArray("Yes","No"), "Yes"); Dialog.addToSameRow();
	Dialog.addMessage("If not, please introduce the following:");
	Dialog.addNumber("Image Calibration (µm/px):", 0.16); Dialog.addToSameRow();
	Dialog.addNumber("Image Binning:", 1);
	
// ====== IMAGE PREPROCESSING ======
	Dialog.addMessage("Image Preprocessing", 16);
	// Background subtraction
	Dialog.addMessage("Background subtraction: Rolling Ball radius (µm)"); 
	Dialog.addNumber("Synapses:", 2);
	// Z-stack processing
	Dialog.addCheckbox("Input are z-stacks?", false); Dialog.addToSameRow();
	Dialog.addChoice("Z handling:", newArray("MaxIP","SumIP"), "MaxIP");
	
// ====== SYNQUANT PARAMETERS ======
	Dialog.addMessage("SynQuant", 16);
	// Size thresholds
	Dialog.addNumber("z-score threshold:", 2); Dialog.addToSameRow();
	Dialog.addNumber("Min synapse size (µm²):", 0.2); Dialog.addToSameRow();
	Dialog.addNumber("Max synapse size (µm²):", 1.2);
	// Intensity thresholds
	Dialog.addNumber("Min fill (0-1):", 0.5); Dialog.addToSameRow();
	Dialog.addNumber("Max WH ratio:", 4); Dialog.addToSameRow();
	// Advanced
	Dialog.addNumber("zScore adjustment:", 2);
	Dialog.show();

/* ============================ VARIABLES ======================================
 *  Read the dialog results in the exact order they were added above.
 =========================================================================== */
// Experiment Parameters
	inDir     = Dialog.getString(); 
	outpuDir = inDir+"Results"+File.separator; File.makeDirectory(outpuDir);
	expPrefix = Dialog.getString();
	firstIdx  = Dialog.getNumber();
	lastIdx   = Dialog.getNumber();
	sType     = Dialog.getChoice();
// Image calibration
	calImg	  = Dialog.getChoice();
	px_um     = Dialog.getNumber() * Dialog.getNumber(); // returns -> calibration * binning
// Background
	rollSyn   = Dialog.getNumber();
// Z-stacks handling
	isZ       = Dialog.getCheckbox();
	zMode     = Dialog.getChoice();
// SynQuant parameters
	zScore    = Dialog.getNumber();
	minSize   = Dialog.getNumber();
	maxSize   = Dialog.getNumber();
	min_0     = Dialog.getNumber();
	max_0     = Dialog.getNumber();
	zScore2   = Dialog.getNumber();
	
// Parameters validation
	if (px_um == 0) { exit("Invalid calibration values"); }


/* ============================ PROCESSING LOOP ======================================
 *  Main Loop: iterate over image sets form firstidx to lastidx
 *  Per set:
 *  	Check file existance
 *  	Open channels, and scale/project if needed
 *  	Preprocess (background substraction
 *  	SynQuant detection and ROI filtering
 *  	Append results to summary table
 =================================================================================== */

	for (i = firstIdx; i <= lastIdx; i++) {
		idStr = d2(i); // 2-digit index (01, 02, ...)
		if (sType == "Excitatory") { chList = excitatory; }
		else if (sType == "Inhibitory") { chList = inhibitory; }
		
	// Image checking: build expected paths and ensure files exist
		ok = true;
		pathList = newArray(chList.length);
		for (c=0; c<chList.length; c++) {
			p = inDir + expPrefix + chList[c] + "_" + idStr + ".tif";
			pathList[c] = p;
			if (!File.exists(p)) {
				// Warn and exit: user probably selected wrong synapse type or prefix
				print("[WARN] Missing: "+p+" (skipping set "+idStr+")");
				ok = false;
				exit("Input data does not match with your parameters selection");
			}
		}
		if (!ok) { continue; }
		
	// Results directory for each image set
		outDir = outpuDir+i+File.separator;
		if(!File.exists(outDir)){ File.makeDirectory(outDir); }
		
	// Opening and preprocessing: open channels and optionally project z-stacks
		imgs = newArray(chList.length);
		for (c = 0; c < pathList.length; c++) {
			open(pathList[c]);
			title = getTitle();
			rawID = getImageID();
			
			
		// Scale images: if images are not already calibrated, set scale using px_um
			if (calImg != "Yes") { imgArea = scale_current(px_um, rawID);
 }
			else { getStatistics(imgArea, mean, min, max, std, histogram); } 
			
		// z-Handling (optional): create a projection for stacks, then select the new image
			if (isZ) {
				if (zMode == "MaxIP") { run("Z Project...", "projection=[Max Intensity]"); selectImage(rawID); close(); }
				else if (zMode == "MeanIP") { run("Z Project...", "projection=[Sum Slices]"); selectImage(rawID); }
				// After projection, pick the last opened image (project result)
				titles = getList("image.titles");
				selectImage(titles[lengthOf(titles)-1]);
			}
			
			else if (!isZ) {
				getDimensions(width, height, channels, slices, frames);
				// If user didn't mark Z, but file has multiple slices, abort to avoid wrong processing
				if (slices > 1) { cleanWorkspace(); exit("Input data are Z-stacks and Z-stack option was not checked"); }
			}
			
			rename(logical[c]);
		// rename to a predictable logical name
			imgs[c] = getTitle();
	// store the new title in the imgs array
		}

		tileImages(); // arrange open images on screen for convenience
	
		// Background subtraction
		getPixelSize(unit, pxSize, pixelHeight);
		selectImage("postSyn");
		rb_sub(rollSyn, pxSize);	// background for post-synaptic channel
		selectImage("preSyn");
		rb_sub(rollSyn, pxSize);	// background for pre-synaptic channel
		
	// SynQuant analysis
	// Check whether pre/post synapse images are still open after preprocessing
		if (!isOpen("preSyn") || !isOpen("postSyn")) {
			warnAndCleanup("Missing required images after preprocessing.", outDir);
			continue;
		}
		
		// Convert synapse size thresholds from µm² to pixels (value expected by SynQuant)
		getPixelSize(unit, pixelWidth, pixelHeight);
		minSizeN = umPixConvert(minSize, pixelWidth);
		maxSizeN = umPixConvert(maxSize, pixelWidth);
		
		// Run SynQuant; parameters are composed into a long argument string
		run("SynQuantVid ", "z-score="+zScore+" min="+minSizeN+" max="+maxSizeN+" min_0="+min_0+" max_0="+max_0+" post-synapse=postSyn pre-synapse=preSyn way=intersect dendrite=Null extended=1 z=1 zscore="+zScore2);
		tileImages();
		
	// SynQuant Data saving: collect and store SynQuant outputs (ROIs, images, tables)
		if (isOpen("Synapse Quantification Coefficents Table")) { selectWindow("Synapse Quantification Coefficents Table"); run("Close"); } // Not used downstream

		if (isOpen("Synapse detection results")) {
			selectWindow("Synapse detection results");
			roiManager("save", outDir+"ROIset_Synapse.zip");
			saveAs("Tiff", outDir+"Synapse_detection_results.tif");
			totRes = roiManager("count");	// total detected synapses remaining
			
			close("Synapse detection results");
		}


	// SynQuant derived Data: export summary and features if available
		if (isOpen("Synapse Quantification Summary Table")) {
			selectWindow("Synapse Quantification Summary Table");
			saveAs("Results", outDir+"Synapse_Quantification_Summary.csv");
			run("Close");
		}

		if (isOpen("Synapse Quantification Feature Table")) {
			selectWindow("Synapse Quantification Feature Table");
			totDen = Table.getColumn("Label");
			saveAs("Results", outDir+"Synapse Quantification Feature Table.xls");
			run("Close");
		}
		
		getStatistics(imgArea, mean, min, max, std, histogram);
		
	// Append results into the global "Summary" table (one row per image set)
		selectWindow("Summary");
		Table.set("Experiment", row, expPrefix+i);
		Table.set("Total # synapses", row, totRes);
		Table.set("Synapse density (/10 µm²)", row, (totRes/imgArea)*10);
		Table.update;
		row++;
	// Clean workspace to free memory before next iteration
		cleanWorkspace();
	}

/* ============================ SUMMARY STATISTICS ======================================
 *  After looping all image sets, compute mean / STD / SEM for each metric and save
 *  results.xls containing per-set rows + aggregated statistics.
 ====================================================================================== */
	selectWindow("Summary");
	totalSyn = Table.getColumn("Total # synapses");
	SynperArea = Table.getColumn("Synapse density (/10 µm²)");
	
	Table.set("Experiment", row, "Avg");
	Table.set("Experiment", row+1, "STD");
	Table.set("Experiment", row+2, "SEM");
	
	// For each vector, compute statistics and place in the table (AVG / STD / SEM)
	
	Array.getStatistics(totalSyn, min, max, mean, stdDev);
	Table.set("Total # synapses", row, mean);
	Table.set("Total # synapses", row+1, stdDev);
	Table.set("Total # synapses", row+2, stdDev/sqrt(totalSyn.length));
	
	Array.getStatistics(SynperArea, min, max, mean, stdDev);
	Table.set("Synapse density (/10 µm²)", row, mean);
	Table.set("Synapse density (/10 µm²)", row+1, stdDev);
	Table.set("Synapse density (/10 µm²)", row+2, stdDev/sqrt(SynperArea.length));
	
	Table.update;
	
	saveAs("Results", outpuDir+"Results.xls");
	run("Close");
	
	
/* ============================ FUNCTIONS ===========================================
 *  Utility functions used across the macro.
 ================================================================================ */
	function initWorkspace(){
	// Close images and helper windows, force garbage collection, and set default measurement settings and colors.
		if (nImages > 0) { run("Close All"); } // Close all open images
		close("ROI Manager"); // Clean & close ROIManager
		close("Log"); close("Results"); // close log, results & Data  windows
		closeOpenData();
		call("java.lang.System.gc");
		run("Collect Garbage");	// Free up unused ram memory
		run("Bio-Formats Macro Extensions"); // Activate Bio-Formats Macro Extensions
		run("Set Measurements...", "area mean min perimeter shape area_fraction limit display redirect=None decimal=3");
		setForegroundColor(255, 255, 255);
		setBackgroundColor(0, 0, 0);
		row = 0;
		return row;
	}

	function cleanWorkspace(){
	// Similar to init but lighter: used between image-set iterations
		if (nImages > 0) { run("Close All"); } // Close all open images
		close("ROI Manager"); // Clean & close ROIManager
		close("Log"); close("Results"); // close log, results & Data  windows
		run("Collect Garbage");	// Free up unused ram memory
		call("java.lang.System.gc");
	}
	
	function closeOpenData(){
	// Close known SynQuant windows if open and ensure a clean "Summary" table exists
		if (isOpen("Synapse Quantification Summary Table")) { close("Synapse Quantification Summary Table"); }
		if (isOpen("Synapse Quantification Coefficents Table")) { close("Synapse Quantification Coefficents Table"); }
		if (isOpen("Synapse Quantification Feature Table")) { close("Synapse Quantification Feature Table"); }
		if (isOpen("Synapse detection results")) { selectWindow("Synapse detection results"); run("Close"); }
		if (isOpen("Summary")) {
			selectWindow("Summary");
			// If table has rows, delete them (reset)
			if (Table.size > 0) { for (t = Table.size; t >= 0; t--) { Table.deleteRows(0, t, "Summary"); } }
		}
		else { Table.create("Summary"); }
 // create if missing
		tablePosition();
	}
	
	function d2(n){
	// Return a 2-digit string for indices: 1 -> "01", 10 -> "10"
		if (n < 10) { return "0"+n; }
		return ""+n;
	}

	function contains(s, sub) {
	// Contains function (case-insensitive)
		s = toLowerCase(s); sub = toLowerCase(sub);
		if (indexOf(s, sub) < 0) { return false }
		return true;
	}
	
	function containsAny(s, arr) {
	// Return true if any element in arr is contained (case-insensitive) in s
		for (k = 0; k < arr.length; k++) {
	        if (contains(s, arr[k])) { return true; }
		}
		return false;
	}

	function rb_sub(rolling, pxSize){
	// Substract background using rolling ball
	// rolling: rolling ball radius in µm, pxSize: µm per pixel
	// round(rolling/pxSize) converts µm -> pixels, required by the Subtract Background command
		run("Subtract Background...",  "rolling="+round(rolling/pxSize));
	}

	function scale_current(f,rawID){
	// This sets the image scale so measurements are in µm and not pixels
	// f: pixel size in µm, rawID: numeric ID of the opened image
		run("Set Scale...", "distance=1 known="+f+" unit=um");
		getStatistics(imgArea, mean, min, max, std, histogram);
		return imgArea;
	}
	
	function tileImages() {
	// Organize open images across the screen to improve visualization during processing.
	// This function computes a layout based on screen resolution and number of images.
		n = nImages();
	// Scren resolution
		scrW = screenWidth;
		scrH = screenHeight;
	// Screen margins
		marginX = scrW * 0.05;
		marginY = scrH * 0.05;
	// Usable screen space
		usableWidth = scrW - 2*marginX;
		usableHeight = scrH - 2*marginY;
	// Image size (all images share same height; width is divided)
		imgWidth = usableWidth / n;
		imgHeight = usableHeight;
		
	// Image positioning loop: compute X/Y and set each window location and size
		for (img = 0; img < n; img++) {
			selectImage(img+1);
			x = marginX + img * imgWidth;
			y = marginY;                 
			setLocation(x, y, imgWidth, imgHeight);
		}
	}
	
	function tablePosition() {
	// Position the Summary table in the lower-left quadrant of the screen
		selectWindow("Summary");
	// Scren resolution
		scrW = screenWidth;
		scrH = screenHeight;
	// Screen margins
		marginX = scrW * 0.05;
		marginY = scrH * 0.4;
	// Table positioning
		bottomLeftX = marginX;
		bottomLeftY = scrH - marginY;
		setLocation(bottomLeftX, bottomLeftY);
	}
	
	function warnAndCleanup(msg, outDir){
	// Print a warning and attempt to save the current log for debugging, then close everything.
		print("[WARN] "+msg);
	// Try to save whatever is open for debugging
		if (isOpen("Log")) { File.saveString(getInfo("log"), outDir+"debug_log.txt"); }
		run("Close All");
	}
	
	function umPixConvert(synSize, pixSize) {
	/* Convert an area in µm² (synSize) to pixels^2 given pixSize (µm per pixel).
	 *  The formula: number_of_pixels = synSize / pixel_area
	 *  pixel_area = pixSize^2
	 *  Use round() because SynQuant expects integer pixel sizes.
	*/
		pxArea = round(synSize/pow(pixSize, 2));
		return pxArea
	}
	
	
	