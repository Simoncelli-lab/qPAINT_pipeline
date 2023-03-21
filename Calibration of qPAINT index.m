%% Calculation of qPAINT index for 1 docking site via Internal Calibration

% Step 1: 
% Calculates qPAINT index (qi) (i.e., inverse of dark times) for all 
% clusters with a maximum point distance of 150 nm in each ROI using the  
% information obtained after running DBSCAN cluster analysis on the 
% DNA-PAINT single-molecule localisation microscopy data. 

% In short, time stamps (frame number) of localisations within the same cluster
% are used to reconstruct the sequence of dark times per cluster as continuous 
% frame times that did not contain an event. All the dark times per cluster 
% are pooled and used to obtain a normalised cumulative histogram of the 
% dark times which is then fitted to 1 - exp(t/tau_d). This allows to estimate
% the dark time, tau_d, per cluster. 

% This step will generate a new .csv file in each of the subfolders where the
% ROIs are stored with the "_Internal Calibration_150" add-on the file
% name. These files are structured as:

% Column 1: clusterID
% Column 2: locs per cluster
% Column 3: mean frames	
% Column 4: darktimecounts	
% Column 5: darktime_mean 	
% Column 6: muMLE = darktime via cumulative distribution function fit
% Column 7: qPAINT_index for each cluster. The value is provided as 100x
%the true value in 1/s.

% It will also generate a merged file containg all the qPAINT indexes
% (Column 7) of all the files generated in each Cell or FOV subfolders.

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the directory where your the ROI files for each Cell or FOV are
% stored

Directory = '/Users/sabri/Documents/qPAINT analysis'; 

%Define if you have 1 or two colour data sets
NChannels = 1; %Number 1 is the default if having 1 color data

% Set integration time per frame used in image acquisition, units: seconds
integration_time        = 0.1;

% Set the column number that contains the cluster ID information for each
% single molecule localisation. 
N_clusterID = 7; % Number 7 is the default if you have used ROIs and DBSCAN clustering step.

% Define the name of the file that will contain all the information of the
% merged qPAINT values. % Change Channel ID in the file name.
Merged_File_Name = 'Merged_allROIs_Ch1_Internal_Calibration_150_qPAINTindex.txt'
Merged_File_Name_2 = 'Merged_allROIs_Ch2_Internal_Calibration_150_qPAINTindex.txt'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run function Internal Calibration qPAINT index:

Internal_Calibration_qPAINT_index(Directory,integration_time,N_clusterID, Merged_File_Name,Merged_File_Name_2,NChannels);

%% Step 2:

% Go to the directory path where your file 'Merged_allROIs_Ch1_Internal_Calibration_150_qPAINTindex.txt'
% is stored.

% Plot qPAINT histogram for Internal Calibration and fit the data using
% Curve Fitting Matlab APP. Use as input: x = bins; y = counts
% Use the drop down menu to choose different types of Guassian functions to
% see which fits your data the best (i.e., two, three (or more) multi-peak Gauss function.

% Look at the parameters of the fit, b1 will be the first qPAINT index, 
% corresponding to a cluster of points for one DNA docking site (i.e., 1 protein). 
% Also check the numbers make sense by seeing if the other "b" parameters are multiples 
% of the first (i.e. b2 is double b1). Ideally, you want to fit a custom 
% equation using the gaussian function but replacing b2 with 2*b1, b3 with
% 3*b1 (and so on). Save the graph in sfit format with file, save session and make a note of
% the qPAINT index for one protein ? you will need this in the next
% analysis code.

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qPAINT_index_all = importdata('Merged_allROIs_Ch1_Internal_Calibration_150_qPAINTindex.txt'); % Name of the Merged file for the Internal Calibration
nbis = 100;

%Use the following lines to load the data of Localisations per Cluster if required.
%LocPerCluster_all = importdata('Merged_allROIs_Ch1_Internal_Calibration_150_LocPerCluster.txt'); % Name of the Merged file for the Internal Calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qPAINT_index_all

rowsToDelete = qPAINT_index_all > 10; % filter and remove qPAINT index above 10
qPAINT_index_all(rowsToDelete) = [];

figure(1), hist(qPAINT_index_all, nbis);

histObject = histogram(qPAINT_index_all,nbis);
counts = histObject.Values;
bins = histObject.BinEdges + histObject.BinWidth/2; n=length(bins);bins =bins(1:1:n-1);  

%Use the following lines to plot the histograms of the Localisations per Cluster if required.
%figure(3), hist(LocPerCluster_all,nbis)
