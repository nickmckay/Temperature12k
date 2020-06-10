# gridding-temp12k
Package to grid the Temperature12K database

This code has been developed by Philipp S. Sommer ([@Chilipp](https://github.com/Chilipp)) and is further maintained and developed in the [Chilipp/gridding-temp12k](https://github.com/Chilipp/gridding-temp12k) repository.

The folder contains the scripts to compute a gridded representation of the
Temperature12k LiPD files. The important scripts are

- [scripts/compute-annual-gmst.ipynb](scripts/compute-annual-gmst.ipynb): This
  notebook computes the mean surface temperature for the 6 latitudinal bands
	using grid-cell specific GAMs for
	Kaufman et al., *Holocene global mean surface temperature: A multi-method
	reconstruction approach*, Scientific Data, submitted. Please see the
	corresponding publication and/or the notebook for a more detailed description
- [scripts/grid-temp12k.ipynb](scripts/grid-temp12k.ipynb): This notebook
  (which is not yet finalized) computes the anomaly per grid cell with reference
	to PI.

## Installation
To install the necessary dependencies, please

1. Clone (or [download](https://github.com/Chilipp/gridding-temp12k/archive/master.zip)) this repository
2. install [anaconda](https://conda.io/en/latest/miniconda.html) for your
   platform and run
3. Open the terminal (or Anaconda Prompt on Windows) and navigate to the
   downloaded folder from step 1
4. Install the dependencies via `conda env create`
5. Activate the new environment via `conda activate temperature12k`
6. install the pyleogrid dependency via
   `pip install git+https://github.com/Chilipp/pyleogrid.git@master`
7. run the jupyter notebook server on your machine via `jupyter notebook` (an
	 internet browser will open)
8. and in your internet browser, navigate to the desired jupyter notebook in
   the [scripts](scripts) directory

If you are not familiar with jupyter notebooks, you can find some tutorials
at: [jupyter.org/try](https://jupyter.org/try)
