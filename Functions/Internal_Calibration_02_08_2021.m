%% Converting clustered data into localizations and dark times

clear all
%% User Parameters

% Specify here parameters on how to analyse your data


img_concentration       = 1;        % concentration of your imager buffer solution (in nM)
integration_time        = 0.1;      % integration time used for image acquisition (in seconds)
min_loc                 = 10;       % minimum number of localisations for a cluster to be analysed
min_number_of_darktimes = 10;        % minimum number of 'dark times' for a cluster to be analysed?
min_darktime            = 10;        % minimum dark time duration (in frames)?
filterbymeanFrames      = 0;        % 1 = yes, 0 = no.
FAB_DNA_stoichioimetry  = 1;        % FAB-DNA stoichiometry. If unknown, input 1.
DarkTimes_per_cluster   = 1;        % 1 = analyse based on mean dark time per cluster; 
                                     % 0 = analyse total dark time distribution, regardless of cluster
max_Cluster_Distance = 200;         % max distance btw two points in a cluster in nanometers
                                     
N_clusterID = 7;                     %column number with cluster ID

%% Sorting out  DBSCAN data

% insert .txt filename here
PALMsieverdata =  importdata('Day_3_C1_pMHC+ICAM_pZAP70_FOV3_ROI2_DBSCAN_E20_P10.csv');     % this loads the chosen .txt file 

%In the DBSCAN .txt file, the columns are:
%[frames;xloc;yloc;area;bg;width;dbscan_id;dbscan_type]

PALMsieverdata = sortrows(PALMsieverdata,N_clusterID); %sort data by clusterIDs 
nocluster = find(PALMsieverdata(:,N_clusterID) == -1); %identify non-cluster points (this corresponds to ID = -1 in DBSCAN)
PALMsieverdata(nocluster,:) = []; %remove non-clustered points

%border = find(PALMsieverdata(:,7) == 0);  % use this if you want to remove
%border points of cluster (this corresponds to ID = 0 in DBSCAN)
%PALMsieverdata(border,:) = []; %remove border points

%%  Counting how many points are in each cluster

   x = unique(PALMsieverdata(:,N_clusterID));  %indexes the unique clusters identified
   N = numel(x);                %calculates total number of clusters
   
    count = zeros(N,1);          %creates an empty matrix
    for k = 1:N;                 %from clusterID = 1 to N
      count(k) = sum(PALMsieverdata(:,N_clusterID)==x(k)); %sum the number of points belonging to a particular cluster ID
    end
   
clustercounts = [ x count ];    %displays clusterID to number of points in cluster
%avg_localizations = mean(count);

%% finding the frames and x,y localisations for each clusterID
    

frame_start_idx = 1;                %starting index for extracting frames
clusterid = clustercounts(1,1);     %starting with clusterid 1

for y = 1:N
    
    total_frames = clustercounts(y,2);                                 %how many frames per cluster 
    frame_stop_idx = frame_start_idx + total_frames - 1;               %end index for range of frames to extract
    clusterframes{y,1} = sort(PALMsieverdata(frame_start_idx:frame_stop_idx,1));   %extract frames into cell array
    
    clusterx{y,1} = (PALMsieverdata(frame_start_idx:frame_stop_idx,2));   %extract x coordinate into cell array
    clustery{y,1} = (PALMsieverdata(frame_start_idx:frame_stop_idx,3));   %extract x coordinate into cell array
    
    frame_start_idx = frame_start_idx + total_frames; %add new starting index i.e. from the next cluster
    clusterid = clustercounts(y,1); % moving on to next clusterid

end

clustercounts = num2cell(clustercounts);    %converts cluster counts to cell array
az_masterdata = [clustercounts,clusterframes,clusterx,clustery]; %consolidates frames to respective clusters in masterdata array

%% Extracting dark times for each cluster

%opts = statset('Display', 'final')

for z = 1:N         %this 'for' loop runs through each clusterID 
    
    az_masterdata{z,6} = diff(az_masterdata{z,3});   %finds difference between frames i.e. the 'on' and 'off' time
    indices = find(abs(az_masterdata{z,6})<min_darktime);       %Identifies error 'off' times i.e 'off' time too short; value of 4 frames is arbitrary
    az_masterdata{z,6}(indices) = [];                %removes all error 'off' times from above condition
    
    off = az_masterdata{z,6};             % off = the 'off' time in frames 
    n_off = length(off);     % l = the number of 'off' times
    off = sort(off);
    
    frame_avg = mean(az_masterdata{z,3});
    darktime_counts = length(az_masterdata{z,6});
    darktime_mean = mean(az_masterdata{z,6});
    
    dark_x_cdf = off;
    dark_y_cdf = ((1:n_off)-0.5)' ./ n_off;
    
    f = @(b,dark_x_cdf) 1-exp(-dark_x_cdf./b(1));                                     % Objective Function
    muMLE = fminsearch(@(b) norm(dark_y_cdf - f(b,dark_x_cdf)), darktime_mean);                  % Estimate Parameters
    
    
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
    
    %if darktime_counts == 0 || number_binding_sites>clustercounts ;
          if darktime_counts == 0 || n_off < min_number_of_darktimes || cluster_maxDistance > max_Cluster_Distance ;
        
        z = z + 1
        
          else 
    
    az_masterdata{z,7} = frame_avg;          % collects avg values into datatable
    az_masterdata{z,8} = darktime_counts;
    az_masterdata{z,9} = darktime_mean;
    az_masterdata{z,10} = muMLE; % collects muMLE values into datatable
    az_masterdata{z,11} = qPAINT_index;
    az_masterdata{z,12} = cluster_maxDistance;
    
    end 
end

darktime_mean = cell2mat(az_masterdata(:,9));
qPAINT_index = cell2mat(az_masterdata(:,11));

%conslidating information into cell array 'masterdata' and adding labels

masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance'};
az_masterdata = [masterdata_labels;az_masterdata];

%% Filtering out Dark Times based on User Parameters

%by mean frames
if filterbymeanFrames == 1
    
    Mean_Frames    = az_masterdata(2:end,7);
    Mean_Frames    = cell2mat(Mean_Frames);

    pd_meanFrame = fitdist(Mean_Frames,'Normal');

    mu_meanFrame  = pd_meanFrame.mu
    std_meanFrame = pd_meanFrame.sigma
    
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

loc_percluster        = az_masterdata(2:end,2);
loc_percluster         = cell2mat(loc_percluster);


%% Plotting the histogram of qPAINT_index 
nbis = 50;

rowsToDelete = qPAINT_index_filtered > 10;
qPAINT_index_filtered(rowsToDelete) = [];

rowsToDelete = loc_percluster > 1000 ;
loc_percluster(rowsToDelete) = [];

%%
nbis = 200;
%rowsToDelete = BSinternalcalibrationdarktimes > 5;
%BSinternalcalibrationdarktimes(rowsToDelete) = [];

figure(1), hist(qPAINT_index_filtered, nbis);
%figure(2), hist(darktimes, nbis);
%figure(3), hist(loc_percluster,nbis)

histObject = histogram(qPAINT_index_filtered,nbis);
counts = histObject.Values;
bins = histObject.BinEdges + histObject.BinWidth/2; n=length(bins);bins =bins(1:1:n-1);  

%% Saving data
data_to_export = az_masterdata(:,[1,2,7,8,9,10,11,]);
writecell(data_to_export,'Day_3_C3_aCD3+aCD28_pZAP70_FOV3_ROI1_Internal_calibration_200nm.txt')
%%

hist_bins = [0:2:500]; hist_bins=hist_bins';
h = hist(darktimes, hist_bins);h=h';


