
clear all 
protiendata = readtable("YOUR_FILE_HERE.txt");
protiendata = table2array(protiendata);
pointsize = 5; 
figure
scatter(protiendata(:,1),protiendata(:,2), pointsize, 'filled')