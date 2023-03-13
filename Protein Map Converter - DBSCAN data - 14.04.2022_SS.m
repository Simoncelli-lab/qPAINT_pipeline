%% Converting clustered data into localizations and dark times

clear all
%% User Parameters

% Specify here parameters on how to analyse your data


img_concentration       = 0.1;        % concentration of your imager buffer solution (in nM)
integration_time        = 0.1;      % integration time used for image acquisition (in seconds)
min_loc                 = 10;       % minimum number of localisations for a cluster to be analysed
min_number_of_darktimes = 10;        % minimum number of 'dark times' for a cluster to be analysed?
min_darktime            = 10;        % minimum dark time duration (in frames)?
filterbymeanFrames      = 0;        % 1 = yes, 0 = no.
filterbyframesd         = 0;        %localisations with small standard deviation in frame number are removed 
Frame_number            = 15000;    % number of frames aquired
filterbydarktimesd      = 0;        % clusters with high standard deviation in the dark time will be filtered out 
FAB_DNA_stoichioimetry  = 1;        % FAB-DNA stoichiometry. If unknown, input 1.
DarkTimes_per_cluster   = 0;        % 1 = analyse based on mean dark time per cluster; 
                                     % 0 = analyse total dark time distribution, regardless of cluster
                                     
%darktime_1BS = 58.945;                  % if knonw dark time for 1 DNA docking site in seconds
NumberLocalizations_1BS = 36;        % if knonw number of localisations per cluster for 1 DNA docking site
qPAINT_index_1BS = 0.905;
                                     
N_clusterID = 5;                     %column number with cluster ID


%% Sorting out  DBSCAN data
% insert .txt filename here
PALMsieverdata =  importdata('data_channel2_DBSCAN_E10_P10_correct.txt');     % this loads the chosen .txt file 

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


%% extracting dark times for each cluster

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
 
    %number_binding_sites = round(darktime_1BS./(muMLE));
    number_binding_sites = round(qPAINT_index./qPAINT_index_1BS);
    %number_binding_sites = round(clustercounts./NumberLocalizations_1BS);
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

%conslidating information into cell array 'masterdata' and adding labels

masterdata_labels = {'clusterID','locs per cluster','frames','framesdiff','clusterx', 'clustery', 'mean frames','darktimecounts', 'darktime_mean','muMLE', 'qPAINT_index','maxClusterDistance', 'number_binding_sites', 'ClusterCentroidsx','ClusterCentroidsy', 'frame std', 'darktime_std'};
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

%% Saving data

data_to_export = az_masterdata(:,[1,2,7,8,10,11,12,13]);
writecell(data_to_export,'slide1_2022.06.24_C4_FOV3_ROI3.txt')

%data_to_export2 = dimer_distance;
%writematrix(data_to_export2,'Week3_Chamber3_FOV2_ROI3_DimerDistances_qINDEX_1.23_noframefilter.txt')
    
%% Plotting BS

%PALMsieverdata2 = readtable("P2Y2_DBSCAN_E10_P10.txt");
%PALMsieverdata2 = table2array(PALMsieverdata2);

%locx = PALMsieverdata(:,2);
%locy = PALMsieverdata(:,3);

bsx    = az_masterdata(2:end,14);
bsx     = cell2mat(bsx);

bsy    = az_masterdata(2:end,15);
bsy     = cell2mat(bsy);


sz = 20; tz =2; 
    figure(1)
    
    %scatter(PALMsieverdata2(:,2),PALMsieverdata2(:,3));pbaspect([1 1 1])
    %hold on
    %scatter (locx, locy,tz,'b','filled');pbaspect([1 1 1])
    %hold on
    scatter (bsx, bsy ,sz,'g','filled');pbaspect([1 1 1])
    %hold on
    %scatter(brushedData(:,1),brushedData(:,2),sz,'g','filled');pbaspect([1 1 1])
    
%% Saving Protein Maps x,y protein map data.
l = length(bsx)
sdummy = 20.*ones(l,1);
A = [bsx bsy sdummy];

writematrix(A,'slide1_2022.06.24_C4_FOV3_ROI3_qPAINT_1.036.csv')
