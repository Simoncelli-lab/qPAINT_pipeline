%% Step 1: Protein Clustering Quantification

% Calculate different parameters associated to protein cluster
% quantification to compare conditions (i.e. KO vs WT; activated vs non-activated, etc). 

% The code calculates the following parameters: 
%(1) Total Number of Proteins : 'TotalProteins.txt'
%(2) Number of Clusters (cluster defined as > 3 proteins): 'NumberofClusters.txt'
%(3) Percentage of proteins in clusters: 'PercetangeClusteredProteins.txt'
%(4) Cluster Size (max distance btw two points in the cluster): 'ClusterSize.txt'
%(5) Cluster Size in terms of number of proteins (number of proteins per
%cluster): 'NumProteinsperCluster.txt'
% EXTRA:(6): Number of different Cluster Types (i.e. Small, Medium, Large):'ClusterTypes.txt'

% In each file you will have in the 1st Column the quantity
% per ROI, and if applicable, in the 2nd Column the quantity per
% micrometer^2. 

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the directory where your the ROI files for each Cell or FOV are
% stored
Directory = '/Users/sabri/Documents/qPAINT analysis';

%Define if you have 1 or two colour data sets
NChannels = 2; %Number 1 is the default if having 1 color data

% Set the size of the ROI in micrometers
size_ROI = 3; 

% Set the cut-off values for the max number of proteins to define a cluster
% as "Small", "Medium" or "Large". Small clusters are defined as having 3
% to N1 proteins; medium clusters are defined as having N1 to N2 proteins, 
% and large clusters are defined as having more than N2 proteins.
N1 = 6 % this is an example
N2 = 12 % this is an example

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run function to generate files with quantitative data:
Protein_Clustering_Quantification(Directory,size_ROI,N1,N2,NChannels);

 %% Step 2: Nearest neighbour analysis
 
 % Nearest neighbour distances (NND) for proteins of the same identify 
 % (i.e, Ch1-Ch1; Ch2-Ch2; ..) or different types of proteins (i.e, Ch1-Ch2) 
 % are calculated using the proteins maps obtained via qPAINT analysis 
 % and k-means clustering (using Protein Map Converter.m code, and files 
 % with the _Protein_Maps.csv extension).
 
% To evaluate the significance of the NND distributions (i.e., Ch1-Ch1, 
% Ch2-Ch2, Ch1-Ch2) this code also generates random positions for each 
% protein type to calculate the random NND distribution for each case. 

% The resulting histogram of the nearest neighbour distances for both
% experimental data sets and the randomly distributed data is then normalized 
% using the total number of NND calculated per ROI to quantify the percentage 
% of the populate with distances smaller than a set threshold value (T1 and T2).

% The code generates the following output files by merging together the 
% information of all the ROIs in all the Cells/FOVs analysed:

% Distribution of Nearest Neighbour Distances (NND) for 1st neighbour:
%(1) NND between Ch1 to Ch1 for data and Rnd:'NND_Ch1.txt' & 'NND_Ch1_Rnd.txt'
%(2) NND between Ch2 to Ch2 for data and Rnd:'NND_Ch2.txt' & 'NND_Ch2_Rnd.txt'
%(3) NND between Ch1 to Ch2 for data and Rnd:'NND_Ch1_Ch2.txt' & 'NND_Ch1_Ch2_Rnd.txt'
%(4) NND between Ch2 to Ch1 for data and Rnd:'NND_Ch2_Ch1.txt' & 'NND_Ch2_Ch1_Rnd.txt'

% Summary descriptors of distributions for data and random distributions:

%(5) Summary of NND between Ch1 to Ch1 for data and Rnd:'Summary_NND_Ch1'
%(6) Summary of NND between Ch2 to Ch2 for data and Rnd:'Summary_NND_Ch1' 
%(7) Summary of NND between Ch1 to Ch2 for data and Rnd:'Summary_NND_Ch1_Ch2.txt' 
%(8) Summary of NND between Ch2 to Ch1 for data and Rnd:'Summary_NND_Ch1_Ch2.txt'

% The summary files are organised as follows, each row corresponds to one ROI:

% Column 1 & 2: Average value of the NND for data and random, respectively
% Column 3 & 4: Median value of the NND for data and random, respectively
% Column 5 & 6: % of NND below thershold T1 (i.e., 25 nm) for data and random, respectively
% Column 7 & 8: % of NND below thershold T2 (i.e., 50 nm) for data and random, respectivel    
    
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the directory where your the ROI files for each Cell or FOV are
% stored. The data for Ch1 and Ch2 should be located in the same folder,
% for example: Cell2_ROI1_Ch1_DBSCAN_E15_P10.txt and
% Cell2_ROI1_Ch2_DBSCAN_E15_P10.txt. For the code run properly, it is 
% important that the ChID is located after the ROI number!!!

Directory = '/Users/sabri/Documents/qPAINT analysis'; 

% Set the size of the ROI in micrometers
size_ROI = 3; 

% Set the cut-off values (T1, T2) for calcualting the % of proteins that are 
% located at a max NND with respect to the other population. 

T1 = 25 % this is an example 
T2 = 50 % this is an example

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run function to generate files with NND information:
NND(Directory,size_ROI,T1,T2);

 %% (OPTIONAL) Step 3: Plotting Histrograms of NND
 
% Loading NND distributions, go to the Path Directory where your data is
% stored:

% Load NND Data:
load 'NND_Ch1_Ch2.txt';
load 'NND_Ch2_Ch1.txt';
load 'NND_Ch1.txt';
load 'NND_Ch2.txt';

% Load NND Random:
load 'NND_Ch1_Ch2_Rnd.txt';
load 'NND_Ch2_Ch1_Rnd.txt';
load 'NND_Ch1_Rnd.txt';
load 'NND_Ch2_Rnd.txt';
 
nbis = 100; % change number of bins 
maxDistance = 500; % change value, this corresponds to cut-off distribution in nm

rowsToDelete1 = NND_Ch1_Ch2 > maxDistance;
NND_Ch1_Ch2(rowsToDelete1) = [];
 
rowsToDelete2 = NND_Ch1_Ch2_Rnd > maxDistance;
NND_Ch1_Ch2_Rnd(rowsToDelete1) = [];
 
rowsToDelete3 = NND_Ch2_Ch1 > maxDistance;
NND_Ch2_Ch1(rowsToDelete3) = [];
 
rowsToDelete4 = NND_Ch2_Ch1_Rnd > maxDistance;
NND_Ch2_Ch1_Rnd(rowsToDelete4) = [];
 
figure(1)
h= histogram(NND_Ch1_Ch2,nbis,'FaceColor','r', 'FaceAlpha', 0.2, 'EdgeColor','none');
hold on
h2= histogram(NND_Ch1_Ch2_Rnd,nbis,'FaceColor','k', 'FaceAlpha', 0.2, 'EdgeColor','none');
title('NND Ch1-Ch2'); 

figure(2)
h3= histogram(NND_Ch2_Ch1,nbis,'FaceColor','g', 'FaceAlpha', 0.2, 'EdgeColor','none');
hold on
h4= histogram(NND_Ch2_Ch1_Rnd,nbis,'FaceColor','k', 'FaceAlpha', 0.2, 'EdgeColor','none');
title('NND Ch2-Ch1'); 

figure(3)
h= histogram(NND_Ch1,nbis,'FaceColor','#7E2F8E', 'FaceAlpha', 0.2, 'EdgeColor','none');
hold on
h2= histogram(NND_Ch1_Rnd,nbis,'FaceColor','k', 'FaceAlpha', 0.2, 'EdgeColor','none');
title('NND Ch1-Ch1'); 

figure(4)
h3= histogram(NND_Ch2,nbis,'FaceColor','#77AC30', 'FaceAlpha', 0.2, 'EdgeColor','none');
hold on
h4= histogram(NND_Ch2_Rnd,nbis,'FaceColor','k', 'FaceAlpha', 0.2, 'EdgeColor','none');
title('NND Ch2-Ch2'); 


