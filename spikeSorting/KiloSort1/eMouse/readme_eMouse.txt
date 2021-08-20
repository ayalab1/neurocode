We made some scripts to generate artificial data with similar statistics to real recordings. 
Run master_eMouse to verify that Kilosort has been installed correctly, as well as to understand
 how the various options are being passed in and used, and what you should be seeing in Phy 
after running Kilosort on your own data. This example has been made intentionally hard, so that 
the are still some clustering errors which you can visualize (and correct for!) in Phy. 

The scripts also measure the accuracy of the algorithm at two different stages (before and 
after the AUTO merges). The accuracy is measured both in terms of the clusters found by Kilosort, 
and in terms of the best achievable accuracy after merging those clusters optimally (with knowledge 
of the ground truth). The last set of results is supposed to mimic the kind of results you would 
get if you were only doing merges in Phy, and you knew exactly what merges are best. 