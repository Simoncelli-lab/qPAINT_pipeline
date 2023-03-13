clear all     
%% Select ROI's in each cell
%DO STEPS 1-4 ONCE; STEPS 5 ONWARDS PER CHANNEL
%(1) Go to APPS (top) and select RegionFinder
%(2) Change: a. Region Size = 2000 nm
%            b. Colummns: x->1, y->2, ChannelID = 0, FrameID-> 4
%            c. Image Size: 41 um @1nm/pix
%(3) Choose 3xROIs per cell and save snapshot of the selected regions (use
%PrtScn and PAINT).
%(4) Save coord file with the ROIs coordinates of the regions
%(5) Generate a copy of the aligned.txt file with your DATA coordinates and name as ChID.csv (ie. 1.csv; 2.csv...) 
%(6) Remeber to move the ChID.csv file to a new Analysis forder for that
%ChID and also copy there the config.txt file and the coords.txt file.
%(7) Change in the coord.txt file (ROI's coord) the number of the first
%column to the ChID and the numbers of the second column to 1,2,3,4....
%(5) Run the following code and remember to change the number of Ch ID!!!
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% User definable parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
foldern=['ROIs']; %Enter the subfolder where the ROIs map will be stored in subfolders
Size_ROI_size=3000;%Enter the size of your ROI (in nm), e.g Size_ROI_size=2000 for a 2000nm x 2000nm ROI
 
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
Headers=[{'x,y,sd,frame'}];
mkdir(foldern)
F=[foldern,'/'];

%%
i=0;
 
 
for CurrentTableID=StartingTable:TotalNumberofDataTables;
    CurrentRegionFileName=strcat(num2str(CurrentTableID), strcat('.csv'));
    FileID = fopen(CurrentRegionFileName);
    datac=dlmread(CurrentRegionFileName,',', 1, 0); %read data from row 1 col 0
    datac=[datac(:,3) datac(:,4) datac(:,2) datac(:,5)];
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
            foldername=[foldername1,'/1.txt'] % change here number of channel depeding on the Img
            dlmwrite(foldername, Headers, 'delimiter','');   
            dlmwrite(foldername, data2, '-append');
                 
         end
end

