# BME3503C - Brain Tumor Segmentation Using Image Processing Techniques
Final Project for Computer Applications for BME

For options 1 and 2 of this code, we used data from the Medical Segmentation Decathlon for brain tumors.
It can be downloaded from this link: http://medicaldecathlon.com/dataaws/

The first folder from this data_set (Task01_BrainTumour) should be in the same directory as your Main.m file. The SliceDecider and nii2png_BT functions should also be in the sime directory as Main.m so that they can properly be called. 

The nii2png_BT function is a slight modification from the nii2png.m code created by Alexander Laurence.
The original code can be found at this link: https://www.mathworks.com/matlabcentral/fileexchange/71567-nifti-image-converter

The SliceDecider often does not choose the best slice for the MRI scan, so there may be issues seeing a tumor. To see a slice where the tumor is visible, I recommend looking at scans 235, 484, 135, and 435 for correct visualizations.

When doing option 3, it is important to properly type in the filepath and filename for the code to run, so please try to type it in the same format as the example. There may be some issues when attempting to input a 3D nifti scan as well, such as with the image registration or conversion to png. 

COLOR.png and Reference_Scan.png are required for the last 2 steps of our analysis of the image, which is the image registration and tumor location detection.
