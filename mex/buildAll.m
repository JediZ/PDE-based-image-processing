%Compilation is only tested on Linux and will most probably not work in Windows until the file separator etc. issues are fixed

%Compile optical flow related functions
mex -I./source/library ./source/Oflow_lhs_elin4_2d.c ./source/library/opticalflowSolvers.c -outdir ./build
mex -I./source/library ./source/Oflow_lhs_llin4_2d.c ./source/library/opticalflowSolvers.c -outdir ./build
mex -I./source/library ./source/Oflow_sor_elin4_2d.c ./source/library/opticalflowSolvers.c -outdir ./build
mex -I./source/library ./source/Oflow_sor_llin4_2d.c ./source/library/opticalflowSolvers.c -outdir ./build

%Compile disparity related functions
mex -I./source/library ./source/Disp_sor_llin4_2d.c ./source/library/disparitySolvers.c -outdir ./build
mex -I./source/library ./source/Disp_sor_llin_sym4_2d.c ./source/library/disparitySolvers.c -outdir ./build
mex -I./source/library ./source/DdiffWeights.c ./source/library/imageDiffusionWeights.c -outdir ./build

%Compile image derivative related functions
mex -I./source/library ./source/FstDerivatives5.c ./source/library/imageDerivatives.c -outdir ./build
mex -I./source/library ./source/SndDerivatives5.c ./source/library/imageDerivatives.c -outdir ./build

%Compile image interpolation related functions
mex -I./source/library ./source/BilinInterp_2d.c ./source/library/imageInterpolation.c -outdir ./build

%Compile level-set related functions
%Ransac is linked against blas and lapack (comes with Matlab)
mex -largeArrayDims -I./source/library ./source/SurfaceEquation.c ./source/library/ransac.c -outdir ./build -lmwblas -lmwlapack

if isunix
	mex CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -I./source/library ./source/CV_solver_2d.c ./source/library/levelsetSolvers.c -outdir ./build
	mex CFLAGS="\$CFLAGS -fopenmp -msse" LDFLAGS="\$LDFLAGS -fopenmp -msse" -I./source/library ./source/AC_solver_2d.c ./source/library/levelsetSolvers.c -outdir ./build
	mex CFLAGS="\$CFLAGS -fopenmp -msse" LDFLAGS="\$LDFLAGS -fopenmp -msse" -I./source/library ./source/Reinit.c ./source/library/levelsetSolvers.c -outdir ./build
else
	mex CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" -I./source/library ./source/CV_solver_2d.c ./source/library/levelsetSolvers.c -outdir ./build
	mex CFLAGS="$CFLAGS -fopenmp -msse" LDFLAGS="$LDFLAGS -fopenmp -msse" --I./source/library ./source/AC_solver_2d.c ./source/library/levelsetSolvers.c -outdir ./build
	mex CFLAGS="$CFLAGS -fopenmp -msse" LDFLAGS="$LDFLAGS -fopenmp -msse" -I./source/library ./source/Reinit.c ./level_sets_c.c ./source/library/levelsetSolvers.c -outdir ./build
end


