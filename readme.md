# Sub-second Volumetric Additive Manufacturing by Digital Incoherent Synthesis of Holographic Light Fields

This source code has been tested on the Windows 10 operating system with `MATLAB R2023a`.

Run `tbatch_1_coarse.m` and `tbatch_2_fine.m` in sequence to perform pattern optimization. 

- The input is a 3D model in TIF format, and its file name is specified by the variable `sample_filename`.
- The output is an optimized projection pattern sequence in MAT file format in the folder named `result_samples/<samplename>/`, where the image sequence is stored as a cell of two-dimensional matrices.
	- The image sequence obtained in the coarse step (corresponding to the beam patterns inside the medium, without holographic optimization) is stored in `<samplename>_coarse_*.mat`.
	- The image sequence obtained in the fine step (corresponding to the DMD projection patterns, with holographic optimization) is stored in `<samplename>_fine_*.mat`.

Demo models are stored in `model/`. It takes approximately 8 minutes to run `model/PointsTest/points p-3.tif` with default parameters, and the expected output is stored in `result_samples/points p-3/points p-3.zip`. 

