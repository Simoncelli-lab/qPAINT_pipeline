
clear all 
PALMsieverdata = readtable("YOUR_FILE_HERE.csv");
PALMsieverdata = table2array(PALMsieverdata);
pointsize = 3; 
figure
scatter(PALMsieverdata(:,2),PALMsieverdata(:,3), pointsize,PALMsieverdata(:,7), 'filled')
