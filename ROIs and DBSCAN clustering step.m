%% Step 1: Convert .hdf5 Picasso localisation file to .txt format for posterior analysis.

% Convert single molecule coordinates stored in .hdf5 files from Picasso 
% to x,y molecular coordinates in .csv format for posterior analysis. The
% final format of the .csv file is:
%
% Column 1 and 2: xloc(nm) and yloc(nm);
% Column 3: frame number
% Column 4: Photons
% Column 5 and 6: Localisation error in x (nm) and Localisation error in y (nm)

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
filename = 'Cell1_ROI1_Ch1.hdf5'; %Enter the name of your .hdf5 localisation file
pixelsize = 130; % Enter the pixel size of your optical configuration in nm (i.e. 1 pixel = 130 nm) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = h5info(filename);
data = h5read(filename,'/locs');
frame = data.frame;
xloc = data.x; xloc = xloc.*pixelsize;
yloc = data.y; yloc = yloc.*pixelsize;
sx = data.lpx; sx = sx.*pixelsize;
sy = data.lpy; sy = sy.*pixelsize;
photons = data.photons;

listLocalizations = [xloc yloc frame photons sx sy];


dlmwrite('Cell1_Ch1.csv',listLocalizations); % CHANGE FILE NAME HERE!

%% Step 2: Select Region of Interest (ROIs) in each cell (typicall 3 ROI's per cell of 3x3 micrometers)

%DO STEPS "A-D" ONCE; STEPS "E" ONWARDS PER CHANNEL IF HAVING MULTI-COLOUR
%DATA (i.e, ChID = Ch1, Ch2, Ch3...).

%(A) Go to APPS (top) and select RegionFinder (https://github.com/quokka79/RegionFinder)
%(B) Change: 1. Region Size = 3000 nm 
%            2. Colummns: x->1, y->2, ChannelID = 0, FrameID-> 3
%            3. Image Size: 60 um @1nm/pix 
%(C) Choose 3 x ROIs per cell and save snapshot of the selected regions for future reference.
%(D) Save coord file with the ROIs coordinates of the regions
%(E) Remember to move the Data_ChID.txt file to a new Analysis forder for that
%ChID and also copy there the coords.txt file.
%(F) Change in the coord.txt file (ROI's coord) the number of the first
%column to the ChID and the numbers of the second column to 1,2,3,4....
%(G) Run the following code and remember to change the number of Ch ID!!!
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
foldern=['ROIs']; %Enter the subfolder where the ROIs map will be stored in subfolders
Size_ROI_size=3000; %Enter the size of your ROI (in nm), e.g Size_ROI_size = 3000 for a 3000nm x 3000nm ROI
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SkipToRegionNumber=1;
[FileNameC,PathNameC] = uigetfile({'*.txt*';'*.*'},'Choose your Coordinates file');
coordinates=dlmread(FileNameC, '\t', 0,0);
NumberOfRegions=size(coordinates,1); 
Size_ROI=Size_ROI_size/2;
if NumberOfRegions < SkipToRegionNumber 
    error(horzcat('You are trying to skip ahead to Region ', num2str(SkipToRegionNumber),' but there are only ', num2str(NumberOfRegions(1)),' region(s) in the coordinates file. Check your settings...'));
end
TotalNumberofDataTables=max(coordinates(:,1));
StartingTable=(coordinates(SkipToRegionNumber,1));
StartingTableCoordsBeginAtLine = find(coordinates(:,1)==StartingTable,1);
Headers=[{'x,y,frame, photons, sdx, sdy'}];
mkdir(foldern)
F=[foldern,'/'];
i=0;
 

for CurrentTableID=StartingTable:TotalNumberofDataTables;
    CurrentRegionFileName=strcat(num2str(CurrentTableID), strcat('.csv'));
    FileID = fopen(CurrentRegionFileName);
    datac=dlmread(CurrentRegionFileName,',', 1, 0); %read data from row 1 col 0
    fclose(FileID);
    coordinates2=[];
    coordinates2=coordinates;
    sheet=find(coordinates2(:,1)~=CurrentTableID);
    coordinates2(sheet,:)=[];
    TotalNumberOfReg=size(coordinates2,1);
        
 
 for CurrentRegionID=1:TotalNumberOfReg
             
            data2=[];
            data2=datac;
            xcrops=find(data2(:,1)<coordinates2(CurrentRegionID,3)-Size_ROI);
            data2(xcrops,:)=[]; 
            xcropl=find(data2(:,1)>coordinates2(CurrentRegionID,3)+Size_ROI);
            data2(xcropl,:)=[]; 
            ycrops=find(data2(:,2)<coordinates2(CurrentRegionID,4)-Size_ROI);
            data2(ycrops,:)=[];
            ycropl=find(data2(:,2)>coordinates2(CurrentRegionID,4)+Size_ROI);
            data2(ycropl,:)=[]; 
            data2(:,1)=data2(:,1)-(coordinates2(CurrentRegionID,3)-Size_ROI);
            data2(:,2)=data2(:,2)-(coordinates2(CurrentRegionID,4)-Size_ROI);
            i=i+1;
            foldername1=[foldern,'/',num2str(i)];
            mkdir(foldername1)
            foldername=[foldername1,'/Cell1_Ch1.csv'] % change here number of channel depeding on the Img
            dlmwrite(foldername, Headers, 'delimiter','');   
            dlmwrite(foldername, data2, '-append');
                 
         end
end

%% Step 3: Run DBSCAN using PALMsiever

%(A) Open PALMsiever in Matlab by navigating to the folder it is stored. (https://github.com/PALMsiever/palm-siever/wiki)
%(B)   1. Load your .txt/.csv file of the ROI (Import, Generic text file)
%      2. Go to Plugins -> Cluster_DBSCAN -> Input "minPts" and "Epsilon"
%      and run.
%(C)   Save DBSCAN output to a .txt file by running the following code,
%      which will save the data in the following format:

% Column 1: frame number% 
% Column 2 and 3: xloc(nm) and yloc(nm);
% Column 4: Photons
% Column 5 and 6: Localisation error in x (nm) and Localisation error in y (nm)
% Column 7 and 8: dbscan_ID, dbscan_type

% DBSCAN ID corresponds to the cluster ID (i.e.,: 1,2,3,4,...). 
% DBSCAN_Type can be -1,0,1 depending on whether the point corresponds to a
% noise (-1), a border cluster point (0) or a core cluster point (1).

Data_DBSCAN = [frame x y photons sdx sdy dbscan_id dbscan_type]; 
 
% IMPORTANT FILE NAME!
% It is recommended to save the output file with the following format:
% (1) Cell Number or FOV Number
% (2) ROI number
% (3) Channel ID (i.e. pseudo-colour,protein type if having multi-colour data). 
% For example: Cell1_ROI1_Ch1 instead of Cell1_Ch1_ROI1. This is important for NND analysis. 

preDBSCANfilename  = "Cell1_ROI1_Ch1";
clip = "_DBSCAN_E_P"; %please input the Eps (E) and MinPts (P) parameters you used to run DBSCAN
file_extension   = ".txt";
 
dlmwrite(strcat(preDBSCANfilename,clip,file_extension),Data_DBSCAN);

%% Step 4 (Extra, NOT needed for qPAINT analysis)

% You can use the following code to plot DBSCAN output to ease
% visualisation of clustering results.

clear all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Data_DBSCAN = load("Cell1_ROI1_Ch1_DBSCAN_E15_P10.txt"); %Enter the name of your DBSCAN output file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pointsize = 3; 
figure(1)
scatter(Data_DBSCAN(:,2),Data_DBSCAN(:,3), pointsize,Data_DBSCAN(:,7), 'filled')


