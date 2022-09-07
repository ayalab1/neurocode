// clean results
run("Clear Results");

//data_title = getString("file name?", "R");
data_title = File.nameWithoutExtension;
//first use polygon to select the region of interest for cell detection
save_path="D:/微云同步盘/George's Cloud/Abroad/NYU/工作/Buzsaki Lab/Projects/Roman/Development/Histology analysis/ImageJ results/"
run("ROI Manager...");
roiManager("Add");
saveAs("ROI",save_path+"/"+data_title+".roi");


run("Find Maxima...", "prominence=50 strict exclude output=List");
saveAs("Results",save_path+"/"+data_title+".csv");

open(save_path+"/"+data_title+".roi");
run("Find Maxima...", "prominence=80 strict exclude output=[Point Selection]");