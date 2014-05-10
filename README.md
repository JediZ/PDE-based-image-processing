Matlab functions related to PDE based image processing including optical flow, disparity and segmentation.

1.0 DIRECTORIES

You can execute the examples by running:
>runme.m

1.1 images
Contains directories containing the test images. One of the directories is called 'middlebury'. In this directory you will find images available from:
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

1.2 matlab
This directory contains all the matlab functions.

1.3 mex
This directory contains all the mex-functions called by some of the Matlab-functions defined in the 'matlab'-directory. You can compile all the mex-functions by running:
>buildAll.m

1.3.1 build
Directory where all the compiled mex-functions are placed.

1.3.2 source
Source code for the mex-functions.

