# Plugins to Phy2
These plugins add additional features to Phy2

The repository has been created from phy1-plugins and made compatible with Phy2 by Thomas Hainmueller. 

## Features
* Reclustering. Reclustering with KlustaKwik 2.0 - dependent on a local version of KlustaKwik, which is provided in the zip file for Windows 10) and python package: pandas. To install write “pip install pandas” in the terminal in your phy environment.
* Outlier removal using the Mahalanobis distance. Standard threshold is 16 standard deviations (adjustable threshold).
* K-means clustering. Standard separation into clusters (adjustable number).
* visualization for spikes violating the refractory period and a filter for the raw traces.

All new features are accessible from the top menu labeled clustering.

## ControllerSettings - Extra columns in ClusterView (not implemented yet)
ControllerSettings also allows you to adjust the number of spike displayed i FeatureView (increased to 15,000) and WaveformView (standard: 300). I recommend to delete the local .phy folder in your data folder, when adjusting these parameters.

## Installation 
* place the content in your plugins directory (~/.phy/), replacing the existing files and plugins folder.
* Installl the python package panda: write “pip install pandas” in the terminal from your phy environment.
* Copy 'tempdir.py' in "*YourPhyDirectory*/phy/utils". The exact path can be found in the error message from phy you will receive when this step has not been completed.

For the time being, you will also have to find the path to your phy installation and copy 'tempdir.py' to *~YourPhyPath*/phy/utils.

### Known issue
If you receive the following error after installation: 'GUI' object has no attribute 'edit_actions' replace SplitShortISI.py with SplitShortISI_v2.py (both files located in the plugins folder).

## How to cite
Please use below DOI for citing these plugins.

<a href="https://zenodo.org/badge/latestdoi/126424002"><img src="https://zenodo.org/badge/126424002.svg" alt="DOI"></a>
