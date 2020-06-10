# Kaufman-etal-2020-ScientificData
This code has been developed by Michael Erb([@michealerb](https://github.com/michaelerb)) and is maintained and developed in the [michaelerb/Kaufman-etal-2020-ScientificData](https://github.com/michaelerb/Kaufman-etal-2020-ScientificData) repository.


This directory contains scripts for processing the climate reconstructions described in Kaufman et al., 2020, Scientific Data, "Holocene global mean surface temperature -- A multi-method reconstruction approach", and for producing Figs. 2-6 in that paper.  The code is written in python3 and is offered as-is.

### Description of scripts:

make_final_data_files.py - This script reads in the climate reconstructions produced using five different methods, restandardizes the data, and produces new files.  The code used to create the climate reconstructions can be found at https://github.com/nickmckay/Temperature12k/tree/master/ScientificDataAnalysis.

GMST_figures.py - This script does some analysis and creates Figs. 3, 4, and 6 in the paper.
 
latitude_band_figures.py - This script does some analysis and creates Figs. 2 and 5 in the paper as well as an additional figure.
