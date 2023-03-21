%% Reconstruction of proteins maps and quantification

% Step 1: Convert single molecule coordinates to x,y molecular coordinates

% Generates x,y molecular coordinates from the DNA-PAINT source data of x,y,sd single
% molecule (SM) localizations for each region of interest (ROI) using the
% information obtained after running DBSCAN cluster analysis. To convert SM
% localizations into molecular coordinate the function uses k-means clustering 
% and the qPAINT index value corresponding to one single docking site
% (i.e., 1 protein). This value can be obtained running the "Calibration of
% qPAINT index.m" code.

% This step will generate two new .csv files in each of the subfolders where the
% ROIs are stored with the "_Protein_Maps" and "_Protein_Info" add-on the file
% names. These files are structured as:

% (1)"_Protein_Maps": Column 1 and 2: x, y coordinates of proteins in nm. 
% (2)"_Protein_Info":

% Column 1: clusterID
% Column 2: locs per cluster
% Column 3: mean frames	
% Column 4: darktimecounts	
% Column 5: darktime_mean 	
% Column 6: muMLE = darktime via cumulative distribution function fit
% Column 7: qPAINT_index for each cluster. The value is provided as 100x
%the true value in 1/s.
% Column 8: maxClusterDistance. This information is usefull to estimate the
% size of each cluster
% Column 9: number of proteins per cluster

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the directory where your the ROI files for each Cell or FOV are
% stored
Directory = '/Users/sabri/Documents/qPAINT analysis'; 

%Define if you have 1 or two colour data sets
NChannels = 1; %Number 1 is the default if having 1 color data

% Set integration time per frame used in image acquisition, units: seconds,
% and the total number of frames recorded:
integration_time        = 0.1;
Frame_number            = 15000;  % number of frames adquired

% Set the column number that contains the cluster ID information for each
% single molecule localisation. 
N_clusterID = 7; % Number 7 is the default if you have used ROIs and DBSCAN clustering step.

% Set the qPAINT index calibration value corresponding to one docking site
qPAINT_index_1BS = 1; % this value corresponds to x100*qPAINT index in units 1/s for Ch1
qPAINT_index_1BS_Ch2 = 1; % this value corresponds to x100*qPAINT index in units 1/s for Ch2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run function to convert single-molecule localisations into proteins maps:
SMLocalizations_to_Protein_Coordinates(Directory,integration_time,N_clusterID, Frame_number,qPAINT_index_1BS, qPAINT_index_1BS_Ch2,NChannels)

%% (OPTIONAL) Step 2: Plot SM localizations and molecular coordinates in the same graph

% Visualize the x,y SM localizations and x,y protein coordinates for a specific ROI 
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Go to the directory path where the DBSCAN and Protein Maps outputs are located.
% Set the file names of the specific ROI that you would like to visualise.
Data_SingleMolecules = load('Cell1_ROI1_Ch1_DBSCAN_E15_P10.txt');
Data_Proteins_Coordinates = load('Cell1_ROI1_Ch1_DBSCAN_E15_P10.txt_Protein_Maps.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
sz = 30; tz =2; 
figure(1)
scatter (Data_SingleMolecules(:,2),Data_SingleMolecules(:,3),tz,'b','filled');pbaspect([1 1 1])
hold on
scatter(Data_Proteins_Coordinates(:,1),Data_Proteins_Coordinates(:,2),sz,'r','filled')
