name = getInfo("image.filename");
selectWindow(name);
run("Subtract Background...", "rolling=50");
run("Split Channels");
selectWindow(name+" (green)");
close();
selectWindow(name+" (blue)");
close();
selectWindow(name+" (red)");

