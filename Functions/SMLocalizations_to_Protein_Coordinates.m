function out = SMLocalizations_to_Protein_Coordinates(Directory,integration_time,N_clusterID, Frame_number,qPAINT_index_1BS, qPAINT_index_1BS_Ch2,NChannels)

% Specify some parameters used for data analyse

    filterbymeanFrames  = 1;   % 1 = yes, 0 = no. Number 1 is default
    filterbyframesd = 0;  % 1 = yes, 0 = no. Number 0 is default
    filterbydarktimesd = 0;  % 1 = yes, 0 = no. Number 0 is default
    min_number_of_darktimes = 10;      % minimum number of 'dark times' for a cluster to be analysed?
    min_darktime            = 10;      % minimum dark time duration (in frames)?
    DarkTimes_per_cluster   = 1;       % 1 = analyse based on mean dark time per cluster;  
    max_Cluster_Distance = 100;        % max distance btw two points in a cluster in nanometers
    clip = "_Protein_Maps"; 
    clip2 = "_Protein_Info";
    file_extension   = ".csv";

    files = dir(fullfile(Directory,'*'));
    directoryNames = {files([files.isdir]).name}; % this is N
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));% this is N

    
for jj=1:length(directoryNames)
    files_Subdir = dir(fullfile(Directory,directoryNames{jj},'*.txt'));% improve by specifying the file extension.
    directoryNames_Subdir = {files_Subdir(~[files_Subdir.isdir]).name}; % files in subfolder.
    NumberofROIs = numel(directoryNames_Subdir);
    
    % Open data files in Cell 1, Cell 2, Cell 3,... subdirectories
     
    if NChannels == 1; 
    
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
    
opts = statset('Display', 'final')

for z = 1:N         %this 'for' loop runs through each clusterID 
    
    az_masterdata{z,6} = diff(az_masterdata{z,3});   %finds difference between frames i.e. the 'on' and 'off' time
    indices = find(abs(az_masterdata{z,6})<min_darktime);       %Identifies error 'off' times i.e 'off' time too short; value of 4 frames is arbitrary
    az_masterdata{z,6}(indices) = [];                %removes all error 'off' times from above condition
    
    off = az_masterdata{z,6};             % off = the 'off' time in frames 
    n_off = length(az_masterdata{z,6});     % l = the number of 'off' times
    off = sort(off);                   
    
    frame_avg = mean(az_masterdata{z,3});
    frame_std = std(az_masterdata{z,3});
    darktime_counts = length(az_masterdata{z,6});
    darktime_mean = mean(az_masterdata{z,6});
    darktime_std = std(az_masterdata{z,6});
    
    dark_x_cdf = off;
    dark_y_cdf = ((1:n_off)-0.5)' ./ n_off;
    
    f = @(b,dark_x_cdf) 1-exp(-dark_x_cdf./b(1));                                     % Objective Function
    muMLE = fminsearch(@(b) norm(dark_y_cdf - f(b,dark_x_cdf)), darktime_mean);                  % Estimate Parameters
    
    darktime_mean = darktime_mean.*integration_time;
    muMLE =muMLE.*integration_time;
    qPAINT_index = 100./muMLE;
    clustercounts = az_masterdata{z,2};
 
    number_binding_sites = round(qPAINT_index./qPAINT_index_1BS);
    number_binding_sites(number_binding_sites==0) = 1 ;
   
    
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
    
    %if darktime_counts == 0 || number_binding_sites>clustercounts ;
          if darktime_counts == 0 || n_off < min_number_of_darktimes || number_binding_sites>clustercounts ;
        
        z = z + 1
        
          else 
        
    az_masterdata{z,7} = frame_avg;          % collects avg values into datatable
    az_masterdata{z,8} = darktime_counts;
    az_masterdata{z,9} = darktime_mean;
    az_masterdata{z,10} = muMLE; % collects muMLE values into datatable
    az_masterdata{z,11} = qPAINT_index;
    az_masterdata{z,12} = cluster_maxDistance;
    az_masterdata{z,13} = number_binding_sites; 
   
    
    Cluster_xy = [clusterx clustery];
    [idx,ClusterCentroids, sumd, D] = kmeans(Cluster_xy,number_binding_sites,'Distance','cityblock','Replicates',5,'Options',opts);
        D2 = min(D,[],2);
        D3 = [idx, D2];  Distances= sortrows(D3,1);
        
     az_masterdata{z,14} = ClusterCentroids(:,1);   
     az_masterdata{z,15} = ClusterCentroids(:,2);   
     az_masterdata{z,16} = frame_std;
     az_masterdata{z,17} = darktime_std;
      
    end 
end

    darktime_mean = cell2mat(az_masterdata(:,9));
    qPAINT_index = cell2mat(az_masterdata(:,11));

    %Conslidating information into cell array 'masterdata' and adding labels

     masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance', 'number_binding_sites', 'ClusterCentroidsx','ClusterCentroidsy', 'frame std', 'darktime_std'};
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
   

    %filter by frame standard deviation
if filterbyframesd == 1
   std_frame = cell2mat(az_masterdata(2:end,16));
   minframesd = (Frame_number/100)*10; % filter by 10% of the frames 
   exclude_by_frame_sd = find(std_frame(:,1) < minframesd );
   az_masterdata(exclude_by_frame_sd,:) = [];
end 

    % filter by dark time sd
if filterbydarktimesd == 1 
    std_darktime = cell2mat(az_masterdata(2:end,17));
    maxdarktimesd = (max(std_darktime))-((max(std_darktime)/100)*20);
    exclude_by_darktime_sd = find(std_darktime(:,1) > maxdarktimesd );
    az_masterdata(exclude_by_darktime_sd,:) = [];
end

    %removing empty elements in the cell array
    az_masterdata = az_masterdata(~any(cellfun(@isempty, az_masterdata),2), :);
    
    data_to_export = az_masterdata(:,[1,2,7,8,10,11,12,13]);
    writecell(data_to_export, strcat(Directory_Data,clip2,file_extension));
    
    %recovering protein coordinates obtained via k-means
    
    bsx    = az_masterdata(2:end,14); bsx = cell2mat(bsx);
    bsy    = az_masterdata(2:end,15); bsy = cell2mat(bsy);
    data_to_export2 = [bsx bsy];
    writematrix(data_to_export2, strcat(Directory_Data,clip,file_extension));
    
        end
        
    elseif NChannels == 2;
        
    for j=1:2:NumberofROIs; 
        
        Data_DBSCAN = load(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}));
        Directory_Data = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j});
        Data_DBSCAN = sortrows(Data_DBSCAN,N_clusterID); 
        nocluster = find(Data_DBSCAN(:,N_clusterID) == -1);
        Data_DBSCAN(nocluster,:) = [];
   
        Data_DBSCAN_Ch2 = load(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j+1}));
        Directory_Data_Ch2 = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j+1});
        Data_DBSCAN_Ch2 = sortrows(Data_DBSCAN_Ch2,N_clusterID); 
        nocluster_Ch2 = find(Data_DBSCAN_Ch2(:,N_clusterID) == -1);
        Data_DBSCAN_Ch2(nocluster_Ch2,:) = [];

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
    
opts = statset('Display', 'final')

for z = 1:N         %this 'for' loop runs through each clusterID 
    
    az_masterdata{z,6} = diff(az_masterdata{z,3});   %finds difference between frames i.e. the 'on' and 'off' time
    indices = find(abs(az_masterdata{z,6})<min_darktime);       %Identifies error 'off' times i.e 'off' time too short; value of 4 frames is arbitrary
    az_masterdata{z,6}(indices) = [];                %removes all error 'off' times from above condition
    
    off = az_masterdata{z,6};             % off = the 'off' time in frames 
    n_off = length(az_masterdata{z,6});     % l = the number of 'off' times
    off = sort(off);                   
    
    frame_avg = mean(az_masterdata{z,3});
    frame_std = std(az_masterdata{z,3});
    darktime_counts = length(az_masterdata{z,6});
    darktime_mean = mean(az_masterdata{z,6});
    darktime_std = std(az_masterdata{z,6});
    
    dark_x_cdf = off;
    dark_y_cdf = ((1:n_off)-0.5)' ./ n_off;
    
    f = @(b,dark_x_cdf) 1-exp(-dark_x_cdf./b(1));                                     % Objective Function
    muMLE = fminsearch(@(b) norm(dark_y_cdf - f(b,dark_x_cdf)), darktime_mean);                  % Estimate Parameters
    
    darktime_mean = darktime_mean.*integration_time;
    muMLE =muMLE.*integration_time;
    qPAINT_index = 100./muMLE;
    clustercounts = az_masterdata{z,2};
 
    number_binding_sites = round(qPAINT_index./qPAINT_index_1BS);
    number_binding_sites(number_binding_sites==0) = 1 ;
   
    
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
    
    %if darktime_counts == 0 || number_binding_sites>clustercounts ;
          if darktime_counts == 0 || n_off < min_number_of_darktimes || number_binding_sites>clustercounts ;
        
        z = z + 1
        
          else 
        
    az_masterdata{z,7} = frame_avg;          % collects avg values into datatable
    az_masterdata{z,8} = darktime_counts;
    az_masterdata{z,9} = darktime_mean;
    az_masterdata{z,10} = muMLE; % collects muMLE values into datatable
    az_masterdata{z,11} = qPAINT_index;
    az_masterdata{z,12} = cluster_maxDistance;
    az_masterdata{z,13} = number_binding_sites; 
   
    
    Cluster_xy = [clusterx clustery];
    [idx,ClusterCentroids, sumd, D] = kmeans(Cluster_xy,number_binding_sites,'Distance','cityblock','Replicates',5,'Options',opts);
        D2 = min(D,[],2);
        D3 = [idx, D2];  Distances= sortrows(D3,1);
        
     az_masterdata{z,14} = ClusterCentroids(:,1);   
     az_masterdata{z,15} = ClusterCentroids(:,2);   
     az_masterdata{z,16} = frame_std;
     az_masterdata{z,17} = darktime_std;
      
    end 
end

    darktime_mean = cell2mat(az_masterdata(:,9));
    qPAINT_index = cell2mat(az_masterdata(:,11));

    %Conslidating information into cell array 'masterdata' and adding labels

     masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance', 'number_binding_sites', 'ClusterCentroidsx','ClusterCentroidsy', 'frame std', 'darktime_std'};
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
   

    %filter by frame standard deviation
if filterbyframesd == 1
   std_frame = cell2mat(az_masterdata(2:end,16));
   minframesd = (Frame_number/100)*10; % filter by 10% of the frames 
   exclude_by_frame_sd = find(std_frame(:,1) < minframesd );
   az_masterdata(exclude_by_frame_sd,:) = [];
end 

    % filter by dark time sd
if filterbydarktimesd == 1 
    std_darktime = cell2mat(az_masterdata(2:end,17));
    maxdarktimesd = (max(std_darktime))-((max(std_darktime)/100)*20);
    exclude_by_darktime_sd = find(std_darktime(:,1) > maxdarktimesd );
    az_masterdata(exclude_by_darktime_sd,:) = [];
end

    %removing empty elements in the cell array
    az_masterdata = az_masterdata(~any(cellfun(@isempty, az_masterdata),2), :);
    
    data_to_export = az_masterdata(:,[1,2,7,8,10,11,12,13]);
    writecell(data_to_export, strcat(Directory_Data,clip2,file_extension));
    
    %recovering protein coordinates obtained via k-means
    
    bsx    = az_masterdata(2:end,14); bsx = cell2mat(bsx);
    bsy    = az_masterdata(2:end,15); bsy = cell2mat(bsy);
    data_to_export2 = [bsx bsy];
    writematrix(data_to_export2, strcat(Directory_Data,clip,file_extension));
    
       % Count number of points per cluster

   x_Ch2 = unique(Data_DBSCAN_Ch2(:,N_clusterID));  %indexes the unique clusters identified
   N_Ch2 = numel(x_Ch2);                %calculates total number of clusters
   
    count_Ch2 = zeros(N_Ch2,1);          %creates an empty matrix
for k = 1:N_Ch2;                 %from clusterID = 1 to N
      count_Ch2(k) = sum(Data_DBSCAN_Ch2(:,N_clusterID)==x_Ch2(k)); %sum the number of points belonging to a particular cluster ID
end
   
    clustercounts_Ch2 = [x_Ch2 count_Ch2 ];    %displays clusterID to number of points in cluster
    
    
    % Identify information per cluster ID: number of frames and x,y
    % positions

    frame_start_idx_Ch2 = 1;                %starting index for extracting frames
    clusterid_Ch2 = clustercounts_Ch2(1,1);     %starting with clusterid 1

    clear clusterx_Ch2
    clear clustery_Ch2
    clear clusterframes_Ch2
    clear total_frames_Ch2
    clear frame_stop_idx_Ch2
    
for y = 1:N_Ch2
    total_frames_Ch2 = clustercounts_Ch2(y,2);                                 %how many frames per cluster 
    frame_stop_idx_Ch2 = frame_start_idx_Ch2 + total_frames_Ch2 - 1;               %end index for range of frames to extract
    clusterframes_Ch2{y,1} = sort(Data_DBSCAN_Ch2(frame_start_idx_Ch2:frame_stop_idx_Ch2,1));   %extract frames into cell array
    
    clusterx_Ch2{y,1} = (Data_DBSCAN_Ch2(frame_start_idx_Ch2:frame_stop_idx_Ch2,2));   %extract x coordinate into cell array
    clustery_Ch2{y,1} = (Data_DBSCAN_Ch2(frame_start_idx_Ch2:frame_stop_idx_Ch2,3));   %extract x coordinate into cell array
    
    frame_start_idx_Ch2 = frame_start_idx_Ch2 + total_frames_Ch2; %add new starting index i.e. from the next cluster
    clusterid_Ch2 = clustercounts_Ch2(y,1); % moving on to next clusterid
end

    clustercounts_Ch2 = num2cell(clustercounts_Ch2);    %converts cluster counts to cell array
    az_masterdata_Ch2 = [clustercounts_Ch2,clusterframes_Ch2,clusterx_Ch2,clustery_Ch2]; %consolidates frames to respective clusters in masterdata array

    %Extracting dark times for each cluster
    
opts = statset('Display', 'final')

for z = 1:N_Ch2         %this 'for' loop runs through each clusterID 
    
    az_masterdata_Ch2{z,6} = diff(az_masterdata_Ch2{z,3});   %finds difference between frames i.e. the 'on' and 'off' time
    indices_Ch2 = find(abs(az_masterdata_Ch2{z,6})<min_darktime);       %Identifies error 'off' times i.e 'off' time too short; value of 4 frames is arbitrary
    az_masterdata_Ch2{z,6}(indices_Ch2) = [];                %removes all error 'off' times from above condition
    
    off_Ch2 = az_masterdata_Ch2{z,6};             % off = the 'off' time in frames 
    n_off_Ch2 = length(az_masterdata_Ch2{z,6});     % l = the number of 'off' times
    off_Ch2 = sort(off_Ch2);                   
    
     frame_avg_Ch2 = mean(az_masterdata_Ch2{z,3});
     frame_std_Ch2 = std(az_masterdata_Ch2{z,3});
    darktime_counts_Ch2 = length(az_masterdata_Ch2{z,6});
    darktime_mean_Ch2 = mean(az_masterdata_Ch2{z,6});
    darktime_std_Ch2 = std(az_masterdata_Ch2{z,6});
    
   dark_x_cdf_Ch2 = off_Ch2;
   dark_y_cdf_Ch2 = ((1:n_off_Ch2)-0.5)' ./ n_off_Ch2;
    
    f = @(b,dark_x_cdf_Ch2) 1-exp(-dark_x_cdf_Ch2./b(1));                                     % Objective Function
    muMLE_Ch2 = fminsearch(@(b) norm(dark_y_cdf_Ch2 - f(b,dark_x_cdf_Ch2)), darktime_mean_Ch2);                  % Estimate Parameters
    
    darktime_mean_Ch2 = darktime_mean_Ch2.*integration_time;
   muMLE_Ch2 =muMLE_Ch2.*integration_time;
    qPAINT_index_Ch2 = 100./muMLE_Ch2;
    clustercounts_Ch2 = az_masterdata_Ch2{z,2};
 
    number_binding_sites_Ch2 = round(qPAINT_index_Ch2./qPAINT_index_1BS_Ch2);
    number_binding_sites_Ch2(number_binding_sites_Ch2==0) = 1 ;
   
    
    clusterx_Ch2 = az_masterdata_Ch2{z,4};
    clustery_Ch2 = az_masterdata_Ch2{z,5};    
    numberOfCoords_Ch2 = length(clusterx_Ch2);
    maxDistance_Ch2 = zeros(1, numberOfCoords_Ch2);
    indexOfMax_Ch2 = zeros(1, numberOfCoords_Ch2, 'int32');

    for k = 1 : numberOfCoords_Ch2
    distances_Ch2 = sqrt((clusterx_Ch2-clusterx_Ch2(k)).^2 + (clustery_Ch2-clustery_Ch2(k)).^2);
    [maxDistance_Ch2(k), indexOfMax_Ch2(k)] = max(distances_Ch2);
    end
    cluster_maxDistance_Ch2 = max(maxDistance_Ch2);
    
    %if darktime_counts_Ch2 == 0 || number_binding_sites>clustercounts_Ch2 ;
          if darktime_counts_Ch2 == 0 || n_off_Ch2 < min_number_of_darktimes || number_binding_sites_Ch2>clustercounts_Ch2 ;
        
        z = z + 1
        
          else 
        
    az_masterdata_Ch2{z,7} =  frame_avg_Ch2;          % collects avg values into datatable
    az_masterdata_Ch2{z,8} = darktime_counts_Ch2;
    az_masterdata_Ch2{z,9} = darktime_mean_Ch2;
    az_masterdata_Ch2{z,10} =muMLE_Ch2; % collects muMLE values into datatable
    az_masterdata_Ch2{z,11} = qPAINT_index_Ch2;
    az_masterdata_Ch2{z,12} = cluster_maxDistance_Ch2;
    az_masterdata_Ch2{z,13} = number_binding_sites_Ch2; 
   
    
    Cluster_xy_Ch2 = [clusterx_Ch2 clustery_Ch2];
    [idx,ClusterCentroids, sumd, D] = kmeans(Cluster_xy_Ch2,number_binding_sites_Ch2,'Distance','cityblock','Replicates',5,'Options',opts);
        D2 = min(D,[],2);
        D3 = [idx, D2];  Distances= sortrows(D3,1);
        
     az_masterdata_Ch2{z,14} = ClusterCentroids(:,1);   
     az_masterdata_Ch2{z,15} = ClusterCentroids(:,2);   
     az_masterdata_Ch2{z,16} =  frame_std_Ch2;
     az_masterdata_Ch2{z,17} = darktime_std_Ch2;
      
    end 
end

    darktime_mean_Ch2 = cell2mat(az_masterdata_Ch2(:,9));
    qPAINT_index_Ch2 = cell2mat(az_masterdata_Ch2(:,11));

    %Conslidating information into cell array 'masterdata' and adding labels

     masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance', 'number_binding_sites', 'ClusterCentroidsx','ClusterCentroidsy', 'frame std', 'darktime_std'};
     az_masterdata_Ch2 = [masterdata_labels;az_masterdata_Ch2];

    % Filtering out Dark Times based on User Parameters

    %by mean frames
if filterbymeanFrames == 1
    
    Mean_Frames_Ch2   = az_masterdata_Ch2(2:end,7);
    Mean_Frames_Ch2   = cell2mat(Mean_Frames_Ch2);

    pd_meanFrame_Ch2 = fitdist(Mean_Frames_Ch2,'Normal');

    mu_meanFrame_Ch2  = pd_meanFrame_Ch2.mu;
    std_meanFrame_Ch2 = pd_meanFrame_Ch2.sigma;
    
    lowerbound_meanFrame_Ch2 = mu_meanFrame_Ch2-std_meanFrame_Ch2;
    upperbound_meanFrame_Ch2 = mu_meanFrame_Ch2+std_meanFrame_Ch2;

    Mean_Frames_Ch2= [mu_meanFrame_Ch2; Mean_Frames_Ch2];       %adding 1 row to match aa_masterdata

    exclude_by_lowbound_meanFrame_Ch2              = find(Mean_Frames_Ch2(:,1) < lowerbound_meanFrame_Ch2); 
    exclude_by_upperbound_meanFrame_Ch2            = find(Mean_Frames_Ch2(:,1) > upperbound_meanFrame_Ch2); 

    masterdata_length_Ch2 = size(az_masterdata_Ch2);
    masterdata_length_Ch2 = masterdata_length_Ch2(1,1);
    
    exclude_cluster_by_Frames_Ch2 = [exclude_by_lowbound_meanFrame_Ch2;exclude_by_upperbound_meanFrame_Ch2];
    az_masterdata_Ch2(exclude_cluster_by_Frames_Ch2,:) = [];
end

    masterdata_length_Ch2 = size(az_masterdata_Ch2);
    masterdata_length_Ch2 = masterdata_length_Ch2(1,1);

    %extracting number of dark times per cluster
if DarkTimes_per_cluster == 1
    
    %extracting dark times per cluster
    darktimes_Ch2        = az_masterdata_Ch2(2:end,10);
    darktimes_Ch2        = cell2mat(darktimes_Ch2);

    %Removing clusters with [] and NaN dark time values
    darktimes_Ch2 = rmmissing(darktimes_Ch2);

    % If we choose to analyse total distribution of dark times
elseif DarkTimes_per_cluster == 0
        
    total_darktimes_Ch2 = cell2mat(az_masterdata_Ch2(2,4));

for g = 3:masterdata_length
    
        tot_DT_Ch2 = cell2mat(az_masterdata_Ch2(g,4));
        total_darktimes_Ch2 = [total_darktimes_Ch2;tot_DT_Ch2];
    
end

darktimes_Ch2 = total_darktimes_Ch2;
end

    qPAINT_index_filtered = 100./darktimes_Ch2; 
   

    %filter by frame standard deviation
if filterbyframesd == 1
   std_frame_Ch2 = cell2mat(az_masterdata_Ch2(2:end,16));
   minframesd_Ch2 = (Frame_number/100)*10; % filter by 10% of the frames 
   exclude_by_frame_sd_Ch2 = find(std_frame_Ch2(:,1) < minframesd_Ch2);
   az_masterdata_Ch2(exclude_by_frame_sd_Ch2,:) = [];
end 

    % filter by dark time sd
if filterbydarktimesd == 1 
    std_darktime_Ch2 = cell2mat(az_masterdata_Ch2(2:end,17));
    maxdarktimesd_Ch2 = (max(std_darktime_Ch2))-((max(std_darktime_Ch2)/100)*20);
    exclude_by_darktime_sd_Ch2 = find(std_darktime_Ch2(:,1) > maxdarktimesd_Ch2 );
    az_masterdata_Ch2(exclude_by_darktime_sd_Ch2,:) = [];
end

    %removing empty elements in the cell array
    az_masterdata_Ch2 = az_masterdata_Ch2(~any(cellfun(@isempty, az_masterdata_Ch2),2), :);
    
    data_to_export2 = az_masterdata_Ch2(:,[1,2,7,8,10,11,12,13]);
    writecell(data_to_export2, strcat(Directory_Data_Ch2,clip2,file_extension));
    
    %recovering protein coordinates obtained via k-means
    
    bsx_Ch2    = az_masterdata_Ch2(2:end,14); bsx_Ch2 = cell2mat(bsx_Ch2);
    bsy_Ch2    = az_masterdata_Ch2(2:end,15); bsy_Ch2 = cell2mat(bsy_Ch2);
    data_to_export2 = [bsx_Ch2 bsy_Ch2];
    writematrix(data_to_export2, strcat(Directory_Data_Ch2,clip,file_extension));
   
        end 
       
end
end
end


