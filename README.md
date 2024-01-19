# Pipeline for multi-feature analysis of aortic endothelial nuclear shape 
## Salvador, J., … Melzer, M.E.,... Goyal, Y., … Iruela-Arispe, M.L., 2024
Created by MEM on 20231021, last update by MEM 20240118

![pipelineSteps](https://github.com/madelinemelzer/nucleiDysmorphia/assets/121903053/141ec955-0366-437e-a881-fec88755fe2f)

Source image: MaxIP_2023_0702_52week_Male1_568LaminA_Arch_60x_3_cp_masks

## Input files required
Nikon .nd2 file with endothelial nuclei specific staining

## Steps Overview
1. Raw .nd2 file to nuclei-only, maxIP, .tif image (Fiji)
2. Segmentation and quality control of nuclei (Python/Jupyter, using Cellpose)
3. Feature-extraction (Python/Jupyter)
4. Creating UMAP embedding (R)

## Step 1: Raw .nd2 file to nuclei-only, maxIP, .tif image (Fiji)
### Scripts used:
"getNuclearMaxIP.ijm" 
### Directories specified:
Input folder of raw images

Output folder of maxIP images of nuclei-only
### Instructions
1. Open Fiji
2. Navigate to Plugins > Macros > Edit
3. Click Edit and select the macro file (getNuclearMaxIP.ijm). This should open up the macro editing window.
4. Modify the directories to match where the .nd2 files are located (inputDir) and where you would like the output .tif files to be located (outputDir). If you have multiple subfolders, specify those as well. For example:

			baseDir = “.../data/age/7wk/raw/”

			outputDir = “.../data/age/7wk/unsegmented/”
5. Click “Run” in the macros editor to begin the segmentation. Your outputDir will populate with .tif images after a minute or so. This step is the longest step to run in the whole pipeline.

## Step 2: Nuclear segmentation (Python/Jupyter, using Cellpose)
### Scripts used:
“segmentation.ipynb"
### Directories specified:
Input folder of maxIP .tif images of nuclei only

Output folder of segmented, binary .tif images of nuclei only

QC output folder of segmented nuclei but filtered by removing peripheral nuclei cutoff by the image field of view

### Instructions
1. Open segmentation Jupyter notebook, “segmentation.ipynb”
2. Modify the input and output directories to match where the unsegmented .tif files are located (dataPath) and where you would like to the output segmented .tif images to be located (savePath). For example:

			dataPath = ".../data/age/7wk/unsegmented/"

			segPath = ".../data/age/7wk/segmented/"
3. Execute the code cell that imports all necessary libraries.
4. Execute the code cell where the image file names and samples are read in from the dataPath.
5. Running the Cellpose segmentation model

	* Modify the model directory to match where the segmentation model you would like to use is. You can also use a default CellPose model. For Salvador et al., a custom model was trained on the segmented nuclei. For example:
pretrained_model="dysmorphicNucleiCellposeModel"). Cellpose documentation regarding model choice can be found here: https://cellpose.readthedocs.io/en/latest/models.html 

	* Specify the channels, flowThreshold, and cellprobThreshold according to the guidelines outlined in the code cell. These guidelines were adapted from CellPose, and related documentation can be found here: https://cellpose.readthedocs.io/en/latest/settings.html 
6. Execute the code cell where the segmented output masks are saved. Here, you have the option to save the masks, flows, outlines, and *_seg.npy files in different formats. For the purposes of this pipeline, you will only need segmented .tif files of the masks (only leave the “tif=True, # save masks as TIFFs” line uncommented.
7. Execute the code cells that define and run the remove_peripherals() function, which removes nuclei that are cutoff by the field of view of the image (peripheral nuclei). For example:

		qcPath = “.../data/age/7wk/filtered/”

## Step 3: Feature-extraction with metrics (Python/Jupyter)
### Scripts used:
“metrics.ipynb”
### Directories specified
Input folder of segmented, quality-controlled .tif images

Output directory for holding .csv table of all metric tables

### Instructions:
1. Open metric calculation Jupyter notebook, “metrics.ipynb”
2. Execute the code cell that imports all necessary libraries.
3. Execute every code cell that defines a function.
4. Modify input and output directories to reflect the location of your segmented, quality-controlled images and where the results will be located.
5. Execute the code cell which calls the main() function and saves the output to dfResult, df_distance, and then to .csv files.
6. Execute the last cell in the notebook which provides a histogram of the nuclear area. Take note of the upper and lower bounds, which you will use as area cutoffs when making the UMAP. 

## Step 4: Creating UMAP embedding (R)
### Scripts used:
metricPlots.R
### Directories specified:
Input feature table (.csv) from Step 3

Output UMAP embedding with associated metadata (.csv)

Plot directory for all resultant UMAP .svg and .png files
### Instructions:
1. Change your file paths to reflect the locations of the above directories.
2. Edit the lower and upper bounds of the area threshold to reflect the results from Step 3.
3. Execute all the code up until line 129, where the scree plot for the PCAs is displayed. Identify the number of PCAs to use based on where the plot tapers off and no additional variance is captured by adding another PCA. For example, in the following plot, the number of PCAs to choose was 11. 

![pcaScree](https://github.com/madelinemelzer/nucleiDysmorphia/assets/121903053/bc0f186e-0d02-4722-86c1-d1d0c9ca21cd)

4. Edit the number of PCAs that make the UMAP based on the results from the above scree plot.
5. Execute lines 130-138 to create the UMAP dataframe, which you can save with associated metadata.
6. Execute lines 142-160 to plot the UMAP and optionally save results. 
