#PDE Based Image Processing for Matlab

##1 INTRODUCTION
Matlab functions related to PDE based image processing including optical flow, disparity and segmentation.
These functions were done as part of my phd thesis which explains, in details, how the models work, how these can be
discretized and solved efficiently. The thesis is available from http://www.jarnoralli.fi, so if you want to
understand how the models work, I suggest that you take a look there. I have also included tutorials and other
material related to PDE based image processing at my site.

##2 ABOUT THE FUNCTIONS

You can execute all the examples by running
>runme.m

###2.1 OPTICAL FLOW
There are several different optical flow functions based on different mathematical formulas, for example:

  1. Horn&Schunck
  2. Late linearization with image warping
  2. Early linearization (without image warping)

###2.2 STEREO DISPARITY
There are several different functions for calculating stereo disparity maps with different constraints, for example:

  1. Late linearization with image warping
  2. Symmetric flow

###2.3 SEGMENTATION BASED ON DYNAMIC IMPLICIT SURFACES
I have included two different functions for segmenting stereo disparity maps (sparse and dense maps) based on dynamic
implicit surfaces, also known as level sets.

###2.4 GEOMETRIC ACTIVE CONTOUR CODE
There are two different versions of geometric active contour:

  1. Geometric active contour
  2. Geodesic active contour

###2.5 TOTAL VARIATION BASED IMAGE DENOISING
There are two variants of image denoising based on total variation:

  1. Non-linear (4 neareast neighbours)
  2. Anisotropic (8 nearest neighbours)

##3 DIRECTORIES
Here is a brief explanation of the directory structure and what they contain. 

###3.1 images-directory
Includes directories containing the test images. One of the directories is called 'middlebury'. In this directory you will find images available from:

-http://vision.middlebury.edu/flow/
-http://vision.middlebury.edu/stereo/data/

This is an excerp from their site:
"How to cite our datasets:
We grant permission to use and publish all images and disparity maps on this website. However, if you use our datasets, we request that you cite the appropriate paper(s): [1] for the 2001 datasets, [2] for the 2003 datasets, and [3] or [4] for the 2005 and 2006 datasets.

References:
[1]	D. Scharstein and R. Szeliski. A taxonomy and evaluation of dense two-frame stereo correspondence algorithms.
International Journal of Computer Vision, 47(1/2/3):7-42, April-June 2002.
[2]	D. Scharstein and R. Szeliski. High-accuracy stereo depth maps using structured light.
In IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 2003), volume 1, pages 195-202, Madison, WI, June 2003.
[3]	D. Scharstein and C. Pal. Learning conditional random fields for stereo.
In IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.
[4]	H. Hirschmüller and D. Scharstein. Evaluation of cost functions for stereo matching.
In IEEE Computer Society Conference on Computer Vision and Pattern Recognition (CVPR 2007), Minneapolis, MN, June 2007.
".

###3.2 matlab-directory
This directory contains all the matlab functions.

###3.3 mex-directory
This directory contains all the mex-functions called by some of the Matlab-functions defined in the 'matlab'-directory. You can compile all the mex-functions by running:
>buildAll.m

####3.3.1 build-directory
Directory where all the compiled mex-functions are placed.

####3.3.2 source-directory
Source code for the mex-functions.

##4 KNOWN ISSUES
MEX-codes using OpenMP, when compiled with certain versions of gcc, fail to execute properly (segmentation fault). I think the problem is related to stack size, but I have not confirmed this yet. If you have problems using the functions, please do let me know.

Cheers,
Jarno

