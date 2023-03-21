function out = Protein_Clustering_Quantification(Directory,size_ROI,N1,N2,NChannels);

    
    files = dir(fullfile(Directory,'*'));
    directoryNames = {files([files.isdir]).name}; % this is N
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));% this is N
        
for jj=1:length(directoryNames)
    files_Subdir = dir(fullfile(Directory,directoryNames{jj},'*Protein_Info.csv'));
    directoryNames_Subdir = {files_Subdir(~[files_Subdir.isdir]).name}; % files in subfolder.
    NumberofROIs = numel(directoryNames_Subdir);
    
    if NChannels ==1;
    % Open data files in Cell 1, Cell 2, Cell 3,... subdirectories

    for j=1:NumberofROIs; 
        
        Protein_Info =importdata(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}));
        Protein_Info = Protein_Info.data(:,:); 
        Protein_Info = Protein_Info(sum(isnan(Protein_Info),2)==0,:);  % this importdatas the chosen .txt file 
        ProteinperCluster = Protein_Info(:,8);  
        MaxDistanceCluster = Protein_Info(:,7); 
        
        Directory_Data = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j})
  
   % Calculating Parameters associated to Proteining Clustering
  
    % Total number of Proteins per ROI or per micrometer^2;
    TotalProteins = sum(ProteinperCluster);
    TotalProteinsperum = sum(ProteinperCluster)./(size_ROI*size_ROI);

    % Definning a cluster as more than 3 proteins:
    ProteinperCluster(ProteinperCluster<3)=[];
    
    % Total number of clusters with > 3 molecules per ROI or micrometer^2 (cluster density); 
    Cluster_Number = length(ProteinperCluster);
    Cluster_Density = Cluster_Number./(size_ROI*size_ROI); % in micrometer^2
    
    % Percentange of clusters with > 3 proteins
    ProteinsinClusters = sum(ProteinperCluster);
    PercetangeClusteredProteins = ProteinsinClusters/TotalProteins*100;
           
    % Extra Information:
    % Number of small clusters, i.e., defined as composed btw 3 to N1 molecules
    SmallClusters = sum(ProteinperCluster(:) < N1+1);
   
    % Number of medium clusters, i.e., defined as composed btw N1 to N2 molecules
    MediumClusters = sum(ProteinperCluster(:) < N2+1); MediumClusters = MediumClusters-SmallClusters;
   
    % Number of large clusters, i.e., defined as having more than N2 molecules
    LargeClusters = Cluster_Number - MediumClusters-SmallClusters;
   
    %Grouping together number of small, medium and large clusters.
    Cluster_Size_Info = [SmallClusters MediumClusters LargeClusters];  
    Cluster_Size_Info_perum2 = Cluster_Size_Info./(size_ROI*size_ROI);
    
    dlmwrite(fullfile(Directory,'TotalProteins.txt'),[TotalProteins TotalProteinsperum], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumberofClusters.txt'),[Cluster_Number Cluster_Density], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'PercetangeClusteredProteins.txt'),PercetangeClusteredProteins, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterSize.txt'),MaxDistanceCluster, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumProteinsperCluster.txt'),ProteinperCluster, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterTypes.txt'),[Cluster_Size_Info Cluster_Size_Info_perum2], 'delimiter',',','precision',10, '-append');

    end 

    elseif NChannels ==2;
       
    for j=1:2:NumberofROIs; 
        
        Protein_Info =importdata(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}));
        Protein_Info = Protein_Info.data(:,:); 
        Protein_Info = Protein_Info(sum(isnan(Protein_Info),2)==0,:);  % this importdatas the chosen .txt file 
        ProteinperCluster = Protein_Info(:,8);  
        MaxDistanceCluster = Protein_Info(:,7); 
        
        Directory_Data = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j})
        
        Protein_Info_Ch2 =importdata(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j+1}));
        Protein_Info_Ch2 = Protein_Info_Ch2.data(:,:); 
        Protein_Info_Ch2 = Protein_Info_Ch2(sum(isnan(Protein_Info_Ch2),2)==0,:);  % this importdatas the chosen .txt file 
        ProteinperCluster_Ch2 = Protein_Info_Ch2(:,8);  
        MaxDistanceCluster_Ch2 = Protein_Info_Ch2(:,7); 
        
        Directory_Data__Ch2 = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j+1})
  
   % Calculating Parameters associated to Proteining Clustering
  
    % Total number of Proteins per ROI or per micrometer^2;
    TotalProteins = sum(ProteinperCluster);
    TotalProteinsperum = sum(ProteinperCluster)./(size_ROI*size_ROI);
    TotalProteins_Ch2 = sum(ProteinperCluster_Ch2);
    TotalProteinsperum_Ch2 = sum(ProteinperCluster_Ch2)./(size_ROI*size_ROI);

    % Definning a cluster as more than 3 proteins:
    ProteinperCluster(ProteinperCluster<3)=[];
    ProteinperCluster_Ch2(ProteinperCluster_Ch2<3)=[];

    % Total number of clusters with > 3 molecules per ROI or micrometer^2 (cluster density); 
    Cluster_Number = length(ProteinperCluster);
    Cluster_Density = Cluster_Number./(size_ROI*size_ROI); % in micrometer^2
    Cluster_Number_Ch2 = length(ProteinperCluster_Ch2);
    Cluster_Density_Ch2 = Cluster_Number_Ch2./(size_ROI*size_ROI); % in micrometer^2
    
    % Percentange of clusters with > 3 proteins
    ProteinsinClusters = sum(ProteinperCluster);
    PercetangeClusteredProteins = ProteinsinClusters/TotalProteins*100;
    ProteinsinClusters_Ch2 = sum(ProteinperCluster_Ch2);
    PercetangeClusteredProteins_Ch2 = ProteinsinClusters_Ch2/TotalProteins_Ch2*100;
      
    % Extra Information:
    % Number of small clusters, i.e., defined as composed btw 3 to N1 molecules
    SmallClusters = sum(ProteinperCluster(:) < N1+1);
    SmallClusters_Ch2 = sum(ProteinperCluster_Ch2(:) < N1+1);
   
    % Number of medium clusters, i.e., defined as composed btw N1 to N2 molecules
    MediumClusters = sum(ProteinperCluster(:) < N2+1); MediumClusters = MediumClusters-SmallClusters;
    MediumClusters_Ch2 = sum(ProteinperCluster_Ch2(:) < N2+1); MediumClusters_Ch2 = MediumClusters_Ch2-SmallClusters_Ch2;

    % Number of large clusters, i.e., defined as having more than N2 molecules
    LargeClusters = Cluster_Number - MediumClusters-SmallClusters;
    LargeClusters_Ch2 = Cluster_Number_Ch2 - MediumClusters_Ch2-SmallClusters_Ch2;

    %Grouping together number of small, medium and large clusters.
    Cluster_Size_Info = [SmallClusters MediumClusters LargeClusters];  
    Cluster_Size_Info_perum2 = Cluster_Size_Info./(size_ROI*size_ROI);
    Cluster_Size_Info_Ch2 = [SmallClusters_Ch2 MediumClusters_Ch2 LargeClusters_Ch2];  
    Cluster_Size_Info_perum2_Ch2 = Cluster_Size_Info_Ch2./(size_ROI*size_ROI);

    dlmwrite(fullfile(Directory,'TotalProteins.txt'),[TotalProteins TotalProteinsperum], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumberofClusters.txt'),[Cluster_Number Cluster_Density], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'PercetangeClusteredProteins.txt'),PercetangeClusteredProteins, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterSize.txt'),MaxDistanceCluster, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumProteinsperCluster.txt'),ProteinperCluster, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterTypes.txt'),[Cluster_Size_Info Cluster_Size_Info_perum2], 'delimiter',',','precision',10, '-append');

    dlmwrite(fullfile(Directory,'TotalProteins_Ch2.txt'),[TotalProteins_Ch2 TotalProteinsperum_Ch2], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumberofClusters_Ch2.txt'),[Cluster_Number_Ch2 Cluster_Density_Ch2], 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'PercetangeClusteredProteins_Ch2.txt'),PercetangeClusteredProteins_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterSize_Ch2.txt'),MaxDistanceCluster_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NumProteinsperCluster_Ch2.txt'),ProteinperCluster_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'ClusterTypes_Ch2.txt'),[Cluster_Size_Info_Ch2 Cluster_Size_Info_perum2_Ch2], 'delimiter',',','precision',10, '-append');

 
    end 
end
   end 
end


