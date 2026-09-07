### Installation instructions
The pipeline was written in Matlab. Full installation of Matlab takes less than one hour. Matlab version 23.2.0.2859533 (R2023b) Update 10 was used. Other versions of Matlab have not been tested.

You must also download the bio-formats plugin code for Matlab, which can be accessed from https://www.openmicroscopy.org/bio-formats/

This pipeline was created using BioFormats version 6.10.1. Updated versions of BioFormats may change functionality, and have not been tested. Installation of BioFormats takes minutes.

The software has been run on both Windows and Mac.

This script is distributed on the CC BY-SA license, created by Greg Reeves 2026, Chemical Engineering Department, Texas A&M University

### Demos for RICS image analysis

Several demos have been provided. To run these demos, first install Matlab, BioFormats, and the code in this repository and the RICS_pipeline linked above.

If the demos require the demo czi files or Mat files, please download the contents of those respective folders from the data repository (link below), and place them in the corresponding (empty) folders from the github tree.
https://doi.org/10.18738/T8/LLXUJG



# Demo to analyze timecourse of one embryo file at fixed location (as in Figs. 1,2):
To run the demo for timecourse data, run the script sc_demo_timecourse.m. This demo shows the opening of the microscopy (czi) file (primary data), calculation of the time course ACFs and CCF (secondary data), and fitting of models to these ACFs and CCF for the estimation of parameter values (tertiary data) as functions of time.

The demo should take less than 5 minutes to run and the expected output is a series of graphs for one part of nc14 for one embryo, including:
- nuclear concentration vs time
- free concentration vs time
- correlated concentration vs time
- uncorrelated concentration vs time
- dose/response curve of correlated concentration vs free concentration
- histone concentration vs time
- free fraction vs time
- correlated fraction vs time
- uncorrelated fraction vs time
- dose/response curve of uncorrelated concentration vs free concentration



# Demo to analyze one embryo file at multiple locations (with tiles; as in Fig. S1f-h):
To run the demo for the rectangular data, run the script sc_demo_rectangle.m. This demo shows the opening of the microscopy (czi) file (primary data), calculation of the rectangle ACFs and CCF (secondary data), and fitting of models to these ACFs and CCF for the estimation of parameter values (tertiary data) as functions of location.

The demo should take less than 5 minutes to run and the expected output is a series of graphs for one part of nc14 for one embryo, including:
- nuclear concentration vs location
- free concentration vs location
- correlated concentration vs location
- uncorrelated concentration vs location
- dose/response curve of correlated concentration vs free concentration
- histone concentration vs location
- free fraction vs location
- correlated fraction vs location
- uncorrelated fraction vs location
- dose/response curve of uncorrelated concentration vs free concentration



# Demo for plotting the entire spatiotemporal data sets (Figure 1,2):
A set of Matlab "mat" files containing the tertiary data (concentrations and fractions as functions of location) for the entire AP-axis data set (at nc14) is provided. To run the demo to plot the outputs of this data set, run the script sc_demo_AP_dataset.m

The demo should run almost instantly and the expected output is a series of graphs for all embryos including:
- nuclear concentration vs time
- free concentration vs time
- correlated concentration vs time
- uncorrelated concentration vs time
- dose/response curve of correlated concentration vs free concentration
- histone concentration vs time
- free fraction vs time
- correlated fraction vs time
- uncorrelated fraction vs time
- dose/response curve of uncorrelated concentration vs free concentration


...as well as another plot, including:
- nuclear concentration vs AP coordinate
- free concentration vs AP coordinate
- correlated concentration vs AP coordinate
- uncorrelated concentration vs AP coordinate
- dose/response curve of correlated concentration vs free concentration
- histone concentration vs AP coordinate
- free fraction vs AP coordinate
- correlated fraction vs AP coordinate
- uncorrelated fraction vs AP coordinate
- dose/response curve of uncorrelated concentration vs free concentration
