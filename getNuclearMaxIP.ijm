// Set batch mode to suppress all GUI prompts
setBatchMode(true);

// Define the input and output directories
inputDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/7wk/raw/";
outputDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/7wk/nucleiUnsegmented/";

// Get a list of .nd2 files in the input directory
fileList = getFileList(inputDir);

function replaceExtension(filename, newExtension) {
    dotIndex = lastIndexOf(filename, ".");
    if (dotIndex != -1) {
        return substring(filename, 0, dotIndex) + newExtension;
    } else {
        return filename + newExtension;
    }
}

// Iterate through each file
for (i = 0; i < fileList.length; i++) {
    if (endsWith(fileList[i], ".nd2")) {
        // Open the file using Bio-Formats
        open(inputDir + fileList[i]);
        
        // Split channels (assuming the .nd2 has multiple channels)
        // run("Split Channels");
        
        // Select the second channel (assuming the channels are named C1, C2, etc.)
        selectWindow(fileList[i] + " - C=2");
        
        // Perform Maximum Intensity Projection
        run("Z Project...", "projection=[Max Intensity]");
        
        // Save the MaxIP result
        outputPath = outputDir + "MaxIP_" + replaceExtension(fileList[i], ".tif");
		saveAs("Tiff", outputPath);

        // Close the opened images to free up memory
        close("*");
    }
}
