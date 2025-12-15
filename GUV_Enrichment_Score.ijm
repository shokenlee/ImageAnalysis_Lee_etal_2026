Clear_Results_ROImanager();

// Formats
file_format = "nd2";
save_format = "tif";

// Channel numbers
protein_channel = 1;
lipid_channel = 2;

// Gaussian blur sigma
gaussian_sigma = 1;

// Autothreshold method
method = "Li";

// Minimum size of GUV to be measured
min_guv_area = 20;

// Get the directory
dataFolder = getDirectory("Choose the folder you want to process");
//dataFolder = "/Users/ShokenLEE/Desktop/DATA/1-2023 10-12/10-26-23 P-290 GUV AH-AcGFP/tifs/";


// Make a file list and repeat the analysis for all images in the directory
filelist = getFileList(dataFolder);

for (i = 0; i < lengthOf(filelist); i++) {
    file = filelist[i];
   
    if (endsWith(file, file_format)) {
   
    // Open a file and perform analysis for all nuclei
open(dataFolder + file); rename(file);

// Determine how many guvs you need to analyze
waitForUser("Count how many guvs you have to analyze");
n_of_guvs =  0;
while (n_of_guvs ==  0) {
n_of_guvs = getNumber("How many guvs do you want to analyze?", 0);
};

// Make arrays, which measured values are going to be stored in
Rim_values = newArray(n_of_guvs);
Outside_values = newArray(n_of_guvs);
Ratios = newArray(n_of_guvs);
Areas = newArray(n_of_guvs);

for (j = 1; j <= n_of_guvs; j++) {
// Show counter
print(j + "out of " + n_of_guvs + " guvs being in analysis...");
   
    // Crop and split the image
waitForUser("Crop the area of interest");
Stack.getPosition(channel, slice, frame);
run("Duplicate...", "duplicate channels=" + protein_channel + " slices=" + slice); rename("Protein");
selectWindow(file);
run("Duplicate...", "duplicate channels=" + lipid_channel + " slices=" + slice); rename("Lipid");

//------
// Make a binary image from lipid channel that fits the whole vesicle
selectWindow("Lipid"); run("Gaussian Blur...", "sigma=" + gaussian_sigma);
setOption("BlackBackground", true);
run("Auto Threshold", "method=" + method + " white");
run("Fill Holes");
for (k = 0; k < 1; k++) { run("Erode"); }

//------ A1
// Define and measure in ROI
define_ROI("Lipid", "Protein"); measure_in_ROI("Protein", "area mean");
A1 = getResult("Area", 0); M1 = getResult("Mean", 0); T1 = A1* M1;
print(A1);
selectWindow("Protein"); run("Flatten"); saveAs("tif", dataFolder + file + "_ROI-" + j);
Clear_Results_ROImanager();


//------ A2
// Measure in a smaller region
selectWindow("Lipid"); for (k = 0; k < 2; k++) { run("Erode"); }
// Define and measure in ROI
define_ROI("Lipid", "Protein"); measure_in_ROI("Protein", "area mean");
A2 = getResult("Area", 0); M2 = getResult("Mean", 0); T2 = A2* M2;
Clear_Results_ROImanager();

//------ A3
// Measure in a smaller region
selectWindow("Lipid"); for (k = 0; k < 3; k++) { run("Erode"); }
// Define and measure in ROI
define_ROI("Lipid", "Protein"); measure_in_ROI("Protein", "area mean");
A3 = getResult("Area", 0); M3 = getResult("Mean", 0); T3 = A3* M3;
Clear_Results_ROImanager();

//------ A4
// Alter
selectWindow("Lipid"); for (k = 0; k < 10; k++) { run("Dilate");; }
// Define and measure in ROI
define_ROI("Lipid", "Protein"); measure_in_ROI("Protein", "area mean");
A4 = getResult("Area", 0); M4 = getResult("Mean", 0); T4 = A4* M4;
Clear_Results_ROImanager();


//------ A5
// Alter
selectWindow("Lipid"); for (k = 0; k < 5; k++) { run("Dilate");; }
// Define and measure in ROI
define_ROI("Lipid", "Protein"); measure_in_ROI("Protein", "area mean");
A5 = getResult("Area", 0); M5 = getResult("Mean", 0); T5 = A5* M5;
Clear_Results_ROImanager();

// Enrichment score caluculation
Background_mean = M3;
Rim_values[j-1] = ((T1 - T2) / (A1 - A2)) - Background_mean;
Outside_values[j-1] = ((T5 - T4) / (A5 - A4)) - Background_mean;
Ratios[j-1] = Rim_values[j-1] / Outside_values[j-1];
Areas[j-1] = A1;

// Close all images other than the original image
selectWindow(file); close("\\Others");
    }
   
    // Now all the desired nuclei are measured from the given image, store the values to CSV
run("Clear Results");
for (j = 0; j < n_of_guvs; j++) {
   setResult("FileName", j, file);
   setResult("Outside_value", j, Outside_values[j]);
   setResult("Rim_value", j, Rim_values[j]);
   setResult("Ratio", j, Ratios[j]);
   setResult("Area", j, Areas[j]);
}
selectWindow("Results"); saveAs("Measurements", dataFolder + "_" + file + "_Results.csv");
run("Close All"); Clear_Results_ROImanager();
    }
}

// --- Functions ---
function Clear_Results_ROImanager() {
// Clean up the result table and ROI manager
run("Clear Results");
if (roiManager("count") != 0) {
roiManager("deselect");
roiManager("delete");
};
}

function define_ROI(roi_image, target_image) {
run("Set Measurements...", "area mean redirect=" + target_image + " decimal=3");
selectWindow(roi_image); run("Analyze Particles...", "size=" + min_guv_area + "-Infinity include add");
selectWindow(target_image); roiManager("Show None"); roiManager("Show All");
// waitForUser("Check");
}

function measure_in_ROI(imageName, measured_values) {
// Perform
run("Set Measurements...", measured_values + " redirect=" + imageName +" decimal=3");
roiManager("deselect");
roiManager("measure");
}