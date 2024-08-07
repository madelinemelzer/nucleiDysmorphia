// Set batch mode to suppress all GUI prompts
setBatchMode(true);

// Define the base input and output directories
baseDir = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/newImages/";

subfolders = newArray("12wk", "24wk", "52wk"); 

for (subfolderIndex = 0; subfolderIndex < subfolders.length; subfolderIndex++) {
    inputDir = baseDir + subfolders[subfolderIndex] + "/raw/";
    outputDir = baseDir + subfolders[subfolderIndex] + "/unsegmented/";

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
            
            // Initialize variables for dimensions
	        width = 0; height = 0; channels = 0; slices = 0; frames = 0;
	        Stack.getDimensions(width, height, channels, slices, frames);
	
	        // Check if the image is multi-channel
	        if (channels > 1) {
	            // Split channels
	            run("Split Channels");
	        }
	        
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
}
