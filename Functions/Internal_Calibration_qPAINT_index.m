function out = Internal_Calibration_qPAINT_index(Directory,integration_time,N_clusterID, Merged_File_Name,Merged_File_Name_2, NChannels)

% Specify some parameters used for data analyse

    filterbymeanFrames      = 1;   % 1 = yes, 0 = no. Number 1 is default
    min_number_of_darktimes = 10;      % minimum number of 'dark times' for a cluster to be analysed?
    min_darktime            = 10;      % minimum dark time duration (in frames)?
    DarkTimes_per_cluster   = 1;       % 1 = analyse based on mean dark time per cluster;  
    max_Cluster_Distance = 150;        % max distance btw two points in a cluster in nanometers
    clip = "_InternalCalibration_150";
    file_extension   = ".csv";

    files = dir(fullfile(Directory,'*'));
    directoryNames = {files([files.isdir]).name}; % this is N
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));% this is N

    
for jj=1:length(directoryNames)
    files_Subdir = dir(fullfile(Directory,directoryNames{jj},'*.txt'));% improve by specifying the file extension.
    directoryNames_Subdir = {files_Subdir(~[files_Subdir.isdir]).name}; % files in subfolder.
    NumberofROIs = numel(directoryNames_Subdir);
    
    % Open data files in Cell 1, Cell 2, Cell 3,... subdirectories

    for j=1:NumberofROIs; 
        
        Data_DBSCAN = load(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}));
        Directory_Data = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j})
    % Sort data by cluster ID, and remove non-cluster points (N_clusterID =
    % -1).
   
    Data_DBSCAN = sortrows(Data_DBSCAN,N_clusterID); 
    nocluster = find(Data_DBSCAN(:,N_clusterID) == -1);
    Data_DBSCAN(nocluster,:) = [];

   % Count number of points per cluster

   x = unique(Data_DBSCAN(:,N_clusterID));  %indexes the unique clusters identified
   N = numel(x);                %calculates total number of clusters
   
    count = zeros(N,1);          %creates an empty matrix
for k = 1:N;                 %from clusterID = 1 to N
      count(k) = sum(Data_DBSCAN(:,N_clusterID)==x(k)); %sum the number of points belonging to a particular cluster ID
end
   
    clustercounts = [x count ];    %displays clusterID to number of points in cluster
    
    % Identify information per cluster ID: number of frames and x,y
    % positions

    frame_start_idx = 1;                %starting index for extracting frames
    clusterid = clustercounts(1,1);     %starting with clusterid 1

    clear clusterx
    clear clustery
    clear clusterframes
    clear total_frames
    clear frame_stop_idx
    
for y = 1:N
    total_frames = clustercounts(y,2);                                 %how many frames per cluster 
    frame_stop_idx = frame_start_idx + total_frames - 1;               %end index for range of frames to extract
    clusterframes{y,1} = sort(Data_DBSCAN(frame_start_idx:frame_stop_idx,1));   %extract frames into cell array
    
    clusterx{y,1} = (Data_DBSCAN(frame_start_idx:frame_stop_idx,2));   %extract x coordinate into cell array
    clustery{y,1} = (Data_DBSCAN(frame_start_idx:frame_stop_idx,3));   %extract x coordinate into cell array
    
    frame_start_idx = frame_start_idx + total_frames; %add new starting index i.e. from the next cluster
    clusterid = clustercounts(y,1); % moving on to next clusterid
end

    clustercounts = num2cell(clustercounts);    %converts cluster counts to cell array
    az_masterdata = [clustercounts,clusterframes,clusterx,clustery]; %consolidates frames to respective clusters in masterdata array

    %Extracting dark times for each cluster

for z = 1:N         %this 'for' loop runs through each clusterID 
    
    az_masterdata{z,6} = diff(az_masterdata{z,3});   %finds difference between frames i.e. the 'on' and 'off' time
    indices = find(abs(az_masterdata{z,6})<min_darktime);       %Identifies error 'off' times i.e 'off' time too short; value of 4 frames is arbitrary
    az_masterdata{z,6}(indices) = [];                %removes all error 'off' times from above condition
    
    off = az_masterdata{z,6};     % off = the 'off' time in frames 
    n_off = length(off);          % l = the number of 'off' times
    off = sort(off);
    
    frame_avg = mean(az_masterdata{z,3});
    darktime_counts = length(az_masterdata{z,6});
    darktime_mean = mean(az_masterdata{z,6});
    
    dark_x_cdf = off;
    dark_y_cdf = ((1:n_off)-0.5)' ./ n_off;
    
    f = @(b,dark_x_cdf) 1-exp(-dark_x_cdf./b(1));                                  % Objective Function
    muMLE = fminsearch(@(b) norm(dark_y_cdf - f(b,dark_x_cdf)), darktime_mean);    % Estimate Parameters
    
    
    darktime_mean = darktime_mean.*integration_time;
    muMLE =muMLE.*integration_time;
    
    qPAINT_index = 100./muMLE;
    clustercounts = az_masterdata{z,2};
    
    clusterx = az_masterdata{z,4};
    clustery = az_masterdata{z,5};    
    numberOfCoords = length(clusterx);
    maxDistance = zeros(1, numberOfCoords);
    indexOfMax = zeros(1, numberOfCoords, 'int32');

for k = 1 : numberOfCoords
    distances = sqrt((clusterx-clusterx(k)).^2 + (clustery-clustery(k)).^2);
    [maxDistance(k), indexOfMax(k)] = max(distances);
end

    cluster_maxDistance = max(maxDistance);
    
if darktime_counts == 0 || n_off < min_number_of_darktimes || cluster_maxDistance > max_Cluster_Distance ;
   z = z + 1
else 
    
    az_masterdata{z,7} = frame_avg;          % collects avg values into datatable
    az_masterdata{z,8} = darktime_counts;
    az_masterdata{z,9} = darktime_mean;
    az_masterdata{z,10} = muMLE;             % collects muMLE values into datatable
    az_masterdata{z,11} = qPAINT_index;
    az_masterdata{z,12} = cluster_maxDistance;
    
end 
end

    darktime_mean = cell2mat(az_masterdata(:,9));
    qPAINT_index = cell2mat(az_masterdata(:,11));

    %Conslidating information into cell array 'masterdata' and adding labels

    masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance'};
    az_masterdata = [masterdata_labels;az_masterdata];

    % Filtering out Dark Times based on User Parameters

    %by mean frames
if filterbymeanFrames == 1
    
    Mean_Frames    = az_masterdata(2:end,7);
    Mean_Frames    = cell2mat(Mean_Frames);

    pd_meanFrame = fitdist(Mean_Frames,'Normal');

    mu_meanFrame  = pd_meanFrame.mu;
    std_meanFrame = pd_meanFrame.sigma;
    
    lowerbound_meanFrame = mu_meanFrame-std_meanFrame;
    upperbound_meanFrame = mu_meanFrame+std_meanFrame;

    Mean_Frames = [mu_meanFrame; Mean_Frames ];       %adding 1 row to match aa_masterdata

    exclude_by_lowbound_meanFrame              = find(Mean_Frames(:,1) < lowerbound_meanFrame); 
    exclude_by_upperbound_meanFrame            = find(Mean_Frames(:,1) > upperbound_meanFrame); 

    masterdata_length = size(az_masterdata);
    masterdata_length = masterdata_length(1,1);
    
    exclude_cluster_by_Frames = [ exclude_by_lowbound_meanFrame;exclude_by_upperbound_meanFrame];
    az_masterdata(exclude_cluster_by_Frames,:) = [];
end

    masterdata_length = size(az_masterdata);
    masterdata_length = masterdata_length(1,1);

    %extracting number of dark times per cluster
if DarkTimes_per_cluster == 1
    
    %extracting dark times per cluster
    darktimes        = az_masterdata(2:end,10);
    darktimes        = cell2mat(darktimes);

    %Removing clusters with [] and NaN dark time values
    darktimes = rmmissing(darktimes);

    % If we choose to analyse total distribution of dark times
elseif DarkTimes_per_cluster == 0
        
    total_DarkTimes = cell2mat(az_masterdata(2,4));

for g = 3:masterdata_length
    
        tot_DT = cell2mat(az_masterdata(g,4));
        total_DarkTimes = [total_DarkTimes;tot_DT];
    
end

darktimes = total_DarkTimes;
end

    qPAINT_index_filtered = 100./darktimes; 
    loc_percluster = az_masterdata(2:end,2);
    loc_percluster = cell2mat(loc_percluster);

    data_to_export = az_masterdata(:,[1,2,7,8,9,10,11,]);
    writecell(data_to_export, strcat(Directory_Data,clip,file_extension));
    
    end 
end


    for kk=1:length(directoryNames)
    files_Subdir_IntCalibration = dir(fullfile(Directory,directoryNames{kk},'*.csv'));
    directoryNames_Subdir_IntCalibration = {files_Subdir_IntCalibration(~[files_Subdir_IntCalibration.isdir]).name}; 

    if NChannels == 1; 
    % Open Internal Calibration files in Cell 1, Cell 2, Cell 3,... subdirectories

    NumberofROIs=NumberofROIs;
    
    for ji=1:NumberofROIs; 
        Data_IntCalibration = importdata(fullfile(Directory,directoryNames{kk},directoryNames_Subdir_IntCalibration{ji}));
        Data_IntCalibration = Data_IntCalibration.data(:,:); 
        Data_IntCalibration = Data_IntCalibration(sum(isnan(Data_IntCalibration),2)==0,:);  % this importdatas the chosen .txt file 
        qPAINT_index = Data_IntCalibration(:,7);  
        %NumLocperCluster = Data_IntCalibration(:,2);
        
        dlmwrite(fullfile(Directory,Merged_File_Name),qPAINT_index, 'delimiter',',','precision',10, '-append');
        %dlmwrite(fullfile(Directory,Merged_File_Name),NumLocperCluster, 'delimiter',',','precision',10, '-append');
    end
    
    elseif NChannels == 2;
      
    for ji=1:2:NumberofROIs; 
        Data_IntCalibration = importdata(fullfile(Directory,directoryNames{kk},directoryNames_Subdir_IntCalibration{ji}));
        Data_IntCalibration = Data_IntCalibration.data(:,:); 
        Data_IntCalibration = Data_IntCalibration(sum(isnan(Data_IntCalibration),2)==0,:);  % this importdatas the chosen .txt file 
        qPAINT_index = Data_IntCalibration(:,7);  
        %NumLocperCluster = Data_IntCalibration(:,2);
        
        Data_IntCalibration_Ch2 = importdata(fullfile(Directory,directoryNames{kk},directoryNames_Subdir_IntCalibration{ji+1}));
        Data_IntCalibration_Ch2 = Data_IntCalibration_Ch2.data(:,:); 
        Data_IntCalibration_Ch2 = Data_IntCalibration_Ch2(sum(isnan(Data_IntCalibration_Ch2),2)==0,:);  % this importdatas the chosen .txt file 
        qPAINT_index_Ch2 = Data_IntCalibration_Ch2(:,7);  
        %NumLocperCluster_Ch2 = Data_IntCalibration_Ch2(:,2);

        dlmwrite(fullfile(Directory,Merged_File_Name),qPAINT_index, 'delimiter',',','precision',10, '-append');
        %dlmwrite(fullfile(Directory,Merged_File_Name),NumLocperCluster, 'delimiter',',','precision',10, '-append');
        
         dlmwrite(fullfile(Directory,Merged_File_Name_2),qPAINT_index_Ch2, 'delimiter',',','precision',10, '-append');
        %dlmwrite(fullfile(Directory,Merged_File_Name_2),NumLocperCluster_Ch2, 'delimiter',',','precision',10, '-append');
        
    end
end
       
    end
end



