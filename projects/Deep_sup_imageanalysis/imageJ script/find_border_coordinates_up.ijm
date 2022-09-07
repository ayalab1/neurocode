//draw a segmented line and store the coordinates in a csv file

data_title = getInfo("image.title");
//data_title = File.nameWithoutExtension;
data_title=data_title+"_up";
data_title=replace(data_title,".png","");
run("ROI Manager...");
roiManager("Add");
save_path="D:/微云同步盘/George's Cloud/Abroad/NYU/工作/Buzsaki Lab/Projects/Roman/Development/Histology analysis/ImageJ results/"
getSelectionCoordinates( x, y );
saveAs("ROI",save_path+data_title+".roi");
f = File.open(save_path+data_title+".txt"); 
for (i=0; i<x.length; i++) print(f,x[i]);
for (i=0; i<y.length; i++) print(f,y[i]);