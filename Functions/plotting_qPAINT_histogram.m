
clear all

%% import data

All_ROIs = importdata('All_ROIS_pCD3_internal_calibration_100nm.txt');

%% get rid of empty rows and filter out large clusters

rows = any(isnan(All_ROIs),2);
All_ROIs(rows,:) = [];

data = num2cell(All_ROIs); 
my_Table = cell2table(data); 
new_table = my_Table(~(my_Table.data7 > 10 | my_Table.data2 > 1000),:); % filter and remove qPAINT index above 10 and locs per cluster over 1000

%%
qPAINT_index_filtered = new_table.data7;
nbis = 200;
%rowsToDelete = BSinternalcalibrationdarktimes > 5;
%BSinternalcalibrationdarktimes(rowsToDelete) = [];

figure(1), hist(qPAINT_index_filtered, nbis);
%figure(2), hist(darktimes, nbis);
%figure(3), hist(loc_percluster,nbis)

histObject = histogram(qPAINT_index_filtered,nbis);
counts = histObject.Values;
bins = histObject.BinEdges + histObject.BinWidth/2; n=length(bins);bins =bins(1:1:n-1);  
