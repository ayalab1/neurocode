# KilosortWrapper
Allows you to load Kilosort from a .xml and a .dat file compatible with Neurosuite. 

## Installation
Download and add the KilosortWrapper directory to your Matlab path.

## Settings
Most settings are defined in the KilosortConfiguration.m file. Some general settings are defined in the KilosortWrapper file, including: 

* Path to SSD
* CreateSubdirectory: Allows you to save the output files from Kilosort to a sub directory (labeled by data and time).

You can supply a config version input to use another config file (configuration files should be stored in the ConfigurationFiles folder).
 
## Features
Skip channels: To skip dead channels, synchronize the anatomical groups and the spike groups in Neuroscope and remove the dead channels in the spike groups. The synchronization is necessary for maintaining the correct waveform layout in Phy.
Define probe layouts: The wrapper now supports probes with staggered, poly3 and poly5 probe layouts. Open your xml file and define your probe layout in the Notes field (General information). Kilosort assumes a staggered probe layout without any input.

CreateSubdirectory: Allows you to save the output files from Kilosort to a sub directory (labeled by data and time).

## Outputs
The Kilosort wrapper allows you to save the output in Neurosuite and Phy compatible files. 

### Phy (rezToPhy_KSW)


### Neurosuite (Kilosort2Neurosuite)
Creates all classical files used in the Neurosuite format. For this the dat file is filtered, waveforms are extracted and global PCA features are calculated. 
