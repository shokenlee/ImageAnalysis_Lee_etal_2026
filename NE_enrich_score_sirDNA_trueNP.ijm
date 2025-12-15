// For quantifying "NE enrichment score"

Clear_Results_ROImanager();
print("\\Clear");
run("Clear Results");

// ----- Parameters -----
// Formats
file_format = "nd2";
save_format = "tif";

// Channel numbers
query_channel = 1; // e.g. AH-mNG-n|NLS
nuclear_marker_channel = 2; // e.g. sirDNA

// Background value based on manual measurement
background = 112;

// Gaussian blur sigma
gaussian_sigma = 1;

// Auto threshold
method = "Li";
min_area_of_roi = 200;
max_area_of_roi = 500;

// Fraction of area that is considered nucleoplasm out of "Intranucleus"
Nucleoplasm_Area_Fraction = 0.3;

// How many times of erosion to be done to make a "Intra-nucleus" mask
erosion_repeats = 3;

// Measured items in "Set measurement"
measured_values = "area mean";

// Get the directory
dataFolder = getDirectory("Choose the folder you want to process");


// ---- Measurement ---
// Make a file list and repeat the analysis for all images in the directory
filelist = getFileList(dataFolder);
for (i = 0; i < lengthOf(filelist); i++) {
    file = filelist[i];
    if (endsWith(file, file_format)) { 
        doAnalysis(file);
        print("Done: " + file);
        run("Close All");
    } 
}

// --- Functions ---
// --- Primary functions ---

function doAnalysis(file) { 
	// Open a file and perform analysis for all nuclei
	open(dataFolder + file);
	rename(file);
	
	// Count and define how many nuclei you have to analyze
	waitForUser("Count how many nuclei you have to analyze");
	n_of_nuclei =  0;
	while (n_of_nuclei ==  0) {
		n_of_nuclei = getNumber("How many nuclei do you want to analyze?", 0);
	};
	
	// Make arrays, which measured values are going to be stored in
	Areas = newArray(n_of_nuclei);
	Means_Whole = newArray(n_of_nuclei);
	Means_Nucleoplasm = newArray(n_of_nuclei);
	Means_Periphery = newArray(n_of_nuclei);
	NEESs = newArray(n_of_nuclei);
	
	for (i = 1; i <= n_of_nuclei; i++) {
		// Show counter
		print(i + "out of " + n_of_nuclei + " nuclei being in analysis...");

		// Clean up results and ROI manager
		Clear_Results_ROImanager();
		
		// Crop the image to mNG and sirDNA images by user-defined rectangular area and slice
		waitForUser("Crop the area of interest");
		Stack.getPosition(channel, slice, frame);
		run("Duplicate...", "duplicate channels=" + query_channel + " slices=" + slice); rename("Cropped");
		selectWindow(file);
		run("Duplicate...", "duplicate channels=" + nuclear_marker_channel + " slices=" + slice); rename("sirDNA");
		
		// Auto-threshold sirDNA image
		run("Gaussian Blur...", "sigma=" + gaussian_sigma);
		setOption("BlackBackground", true);
		run("Auto Threshold", "method=" + method + " white");

		// Define the 1st ROI from sirDNA image
		selectWindow("sirDNA");
		run("Set Measurements...", "area mean min centroid redirect=" + "sirDNA" +" decimal=3");
		run("Analyze Particles...", "size=" + min_area_of_roi + "-" + max_area_of_roi + "exclude include add");

		// Judges if a ROI wrapping the nucleus was created
		// case 1: ROI was not found
		if (roiManager("count") == 0) {
			
			// No measurement is performed, just save the image and go to another nucleus
			selectWindow("Cropped");
			saveAs(save_format, dataFolder + file + "_" + i + "_Failed");
		}
		
		// case 2: ROI was found. Perform measurement
		else {	
			// Display ROIs on original image and save it
			selectWindow("Cropped");
			roiManager("Show None");
			roiManager("Show All with labels");
			run("Flatten");
			saveAs(save_format, dataFolder + file + "_" + i);
			
			// Make a mask for making a 2nd ROI that encircles intra-nucleus
			roi_of_interest =  0;
			while (roi_of_interest == 0) {
				roi_of_interest = getNumber("Give the ROI number", 0);
			};
			roiManager("select", roi_of_interest-1);
			run("Create Mask");
			selectWindow("Mask");
			for (j = 0; j < erosion_repeats; j++) {
				run("Erode");
			}
			run("Divide...", "value=255"); run("16-bit");
			
			// Make a second ROI based on the eroded binary image
			selectWindow("Mask"); setAutoThreshold("Default dark no-reset");
			run("Analyze Particles...", "size=50-Infinity include add");
			
			// Make measurments for the two ROIs
			measure_in_ROI("Cropped", measured_values);
			
			// Store the values of Area, Mean and Total intensity to new variables A1, M1, T1
			Area_Whole = getResult("Area", roi_of_interest-1);
			Mean_Whole = getResult("Mean", roi_of_interest-1);
			Total_Whole = Area_Whole * Mean_Whole;
			print(Area_Whole, Mean_Whole, Total_Whole);
			
			// Measure "cropped" image using the 2nd ROI
			// Obtain area and Mean
			Area_IntraNucleus = getResult("Area", roi_of_interest);
			Mean_IntraNucleus = getResult("Mean", roi_of_interest);
			Total_IntraNucleus = Area_IntraNucleus * Mean_IntraNucleus;
			
			// Obtain values of nuclear rim mean intensity
			Total_Periphery = Total_Whole - Total_IntraNucleus;
			Area_Periphery = Area_Whole - Area_IntraNucleus;
			Mean_Periphery = Total_Periphery / Area_Periphery - background;
			
			// Obtain the approximate value of "Nucleoplasm" from the histogram
			// Get the histgram of "Cropped" within the 2nd ROI
			selectWindow("Cropped");
			roiManager("Show None"); roiManager("select", 1);
			getHistogram(values, counts, 256);
			// Identify the index where the pixel count is maximal
			// - intentionally omit j = 255 because it can be maximal when nucleoli intensity is saturated
			max_index = 0;
			temporary_max = 0;
			for (j = 0; j < 255; j++) {
				if (counts[j] > temporary_max) {
					temporary_max = counts[j];
					max_index = j;
				};
			};
			print(max_index);
			// Get the average intensity of the pixels in the range of +/- defined steps around the peak
			// - which should be close to the true mean intensity of the nucleoplasm
			sum_value_around_peak = 0;
			sum_count_around_peak = 0;
			steps = 3;
			for (k = max_index - steps; k < max_index + steps; k++) {
				sum_value_around_peak += values[k] * counts[k];
				sum_count_around_peak += counts[k];
			};
			Mean_Nucleoplasm = sum_value_around_peak / sum_count_around_peak - background;
			print(sum_value_around_peak, sum_count_around_peak);
			
			// Obtain the ratio, i.e., NE enrichment score
			NE_enrichment_score = Mean_Periphery / Mean_Nucleoplasm;
		
		    // Store the values to the arrays
		    Areas[i-1] = Area_Whole;
		    Means_Whole[i-1] = Mean_Whole;
		    Means_Nucleoplasm[i-1] = Mean_Nucleoplasm;
		    Means_Periphery[i-1] = Mean_Periphery;
		    NEESs[i-1] = NE_enrichment_score;
		}

	    // Close all images other than the original image
		selectWindow(file);
		close("\\Others");
	}
	
	
	// Now all the desired nuclei are measured from the given image, store the values to CSV
	run("Clear Results");
	for (i = 0; i < n_of_nuclei; i++) {
	    setResult("FileName", i, file);
	    setResult("Mean_whole", i, Means_Whole[i]);
	    setResult("Mean_Nucleoplasm", i, Means_Nucleoplasm[i]);
	    setResult("Mean_periphery", i, Means_Periphery[i]);
	    setResult("Ratio", i, NEESs[i]);
	}
	
	selectWindow("Results");
	saveAs("Measurements", dataFolder + "_" + file + "_Results.csv");

}


// --- Secondary functions ---

function Clear_Results_ROImanager() { 
	// Clean up the result table and ROI manager
	run("Clear Results");
	if (roiManager("count") != 0) {
		roiManager("deselect");
		roiManager("delete");
	};
}

function measure_in_ROI(imageName, measured_values) { 
	// Perform
	run("Set Measurements...", measured_values + " redirect=" + imageName +" decimal=3");
	roiManager("deselect");
	roiManager("measure");
}




