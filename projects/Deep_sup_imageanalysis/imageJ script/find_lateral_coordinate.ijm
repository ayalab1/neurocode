//get the coordnate of a point at the lateral end of CA1 pyramidal layer
getSelectionCoordinates( x, y );
data_title = getTitle();
data_title = replace(data_title,".png","");
data_title=data_title+"_lateral";
f = File.open("D:/微云同步盘/George's Cloud/Abroad/NYU/工作/Buzsaki Lab/Projects/Roman/Development/Histology analysis/ImageJ results/"+data_title+".txt"); 
for (i=0; i<x.length; i++) print(f,x[i]);
for (i=0; i<y.length; i++) print(f,y[i]);