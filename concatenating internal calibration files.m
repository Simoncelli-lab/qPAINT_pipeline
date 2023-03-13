%% Concatenating files

clear all

files = dir('*.txt');
for i=1:length(files)
    load(files(i).name, '-ascii');
end

All_ROIs = [ROI1' ROI2'];All_ROIs = All_ROIs'; %add here the amount of ROIs you have e.g. the next one is ROI3'

dlmwrite('All_ROIS_internal_calibration_200nm.txt',All_ROIs);


