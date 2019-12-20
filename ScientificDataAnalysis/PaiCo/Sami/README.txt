README
-----------

Suuporting online material for the article: "Pairwise Comparisons to Reconstruct Mean Temperature in the Arctic Atlantic Region Over the Last 2000 Years". The package contains the necessary files to run the experiments for comparing the reconstruction methods, as well as reconstructing the temperature signal of the Arctic Atlantic. 

1. Directory structure

barcast		Matlab implementation of BASCAST by Tingley et al. (2010)
common		Common Matlab-utility functions for the reconstruction methods.
data		Data for all reconstrucitons. Contains the proxy data for Arctic Atlantic and 
			area-weighted mean annual instrumental temperature for the Arctic Atlantic area calculated from HadCRUT3. 
lna			Matlab implementation of Li et al. (2010) method.
manuscript	Matlab scripts for running the experiments and reporting the results in the manuscript.
other		Matlab implementations of Method-of-moments (MoM), Ordinary Least Squares (OLS), and Principal Component Regression (PCReg).
paico		Matlab implementation of PaiCo.
regem		Optimized Matlab implementation of RegEM with TTLS. Other regularizations are currently not supported.

Notice that some of the methods use native code (MEX-files) to speed up the processes. Precompiled versions for Windows 32bit and 64bit are available in the directories, but for other operating systems the files need to be compiled by hand.

2. Initialization

The file init.m in the root of the package initializes the code by adding to path all the subdirectories. This is required before using any of the other methods.

If the Brownian motion method is used, you need to download the package Generalized Linear models package for Matlab from:
http://www-stat.stanford.edu/~tibs/glmnet-matlab/
Furthermore, the extracted directory has to be included in the search path. Notice that the Brownian motion method is used to compare the actual reconstruction methods to LASSO-fitted noise form Brownian motion. Therefore, it is not an actual reconstruction method.

3. Qualitative comparison

The methods need to be initialized (see above) before running the quantitative comparisons. Note that all of the methods require the Parallel Computing Toolbox for Matlab. Remember to initialize matlabpool to take use of all of the CPU cores you have. Depending on your computer, the execution of each of experiments may take anywhere from a half a day to several days.

The qualitative comparisons in Section 3.2.1 are run by executing the command:
>> runSimpleExperiment;
When the calculations are finished, the analysis is carried out and corresponding figures are created by running:
>> createSimpleFigures;
This also saves the figures to the folder manuscript/figures

The quantitative comparisons in Section 3.2.2 are run by first initializing the data:
>> initPseudoProxy(datafile);
where datafile is the name of the file to where the data is to be stored. The datafile will take more than 1Gb of space.
The actual calculations are carried by running
>> runPseudoProxy(datafile);
where datafile is the name of the file where the data is stored. When the calculations are finished, the results are analyzed and figures plotted  by running:
>> createPseudoFigures;
This also saves the figures to the folder manuscript/figures

The uncertainty experiments of Section 3.2.2 need separate calculations and are run:
>> runUncertAnalysis;
When the calculations are finished, the analysis is carried out and corresponding figures are created by running:
>> createUncertFigures;
This also save the figures to the folder manuscript/figures

The power spectrum analysis is carried out and figures plotted by running:
>> createErrorFigures;
Again, this saves the figures to the folder manuscript/figures

4. Reconstructing the Arctic Atlantic

The rerun the reconstruction for the Arctic Atlantic, simply run:
>> runArcticAtlantic;
The figures are created by running:
>> createAtlanticFigures;

The uncertainty quantifications of the statements in the manuscript are calculated by running:
>> calcAtlanticUncertainties;

5. Using the reconstruction methods for own purposes

The reconstruction methods are written in a way that they accept the same kind of input. To understand the structure of the input, type:
help rtstruct

And to gain more information about individual methods, run, for example,
help paico

All of the methods provide information about the parameters and their default values. To get that information, type
[options, info] = paico('info');

6. Linux and Mac-users

The CPP source files in directories /common and /paico/private should be compiled to MEX-files to use the methods. If mex-compiler has not been setup, type
mex -setup
and select an appropriate C++ -compiler. Notice that older Matlab versions were shipped with a compiler that does not support C++.

After that, all of the MEX-files can be compiled by running in the root folder the command:
setup
