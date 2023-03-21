function out = NND(Directory,size_ROI,T1,T2);

    files = dir(fullfile(Directory,'*'));
    directoryNames = {files([files.isdir]).name}; % this is N
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));% this is N

for jj=1:length(directoryNames)
    files_Subdir = dir(fullfile(Directory,directoryNames{jj},'*Protein_Maps.csv'));
    directoryNames_Subdir = {files_Subdir(~[files_Subdir.isdir]).name}; % files in subfolder.
    NumberofROIs = numel(directoryNames_Subdir);
    
    % Open data files in Cell 1, Cell 2, Cell 3,... subdirectories

    for j=1:2:NumberofROIs; 
        
        Protein_Maps_Ch1 =importdata(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}));
        Protein_Maps_Ch2 =importdata(fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j+1}));
        Directory_Data = fullfile(Directory,directoryNames{jj},directoryNames_Subdir{j}); 
    
    % Total number of Proteins per ROI for each Channel 
     Protein_Number_Ch1 = length(Protein_Maps_Ch1);
     Protein_Number_Ch2 = length(Protein_Maps_Ch2);

    % Generating random positions for proteins Ch1 and Ch2. 
     Protein_Maps_Ch1_random = size_ROI*1000.*rand(Protein_Number_Ch1,2);
     Protein_Maps_Ch2_random = size_ROI*1000.*rand(Protein_Number_Ch2,2);
  
     % Calculating NND for random distributions:
     
     NND_Ch1_Ch2_Rnd = pdist2(Protein_Maps_Ch1_random,Protein_Maps_Ch2_random ,'euclidean','Smallest',1); NND_Ch1_Ch2_Rnd=NND_Ch1_Ch2_Rnd';
     NND_Ch2_Ch1_Rnd = pdist2(Protein_Maps_Ch2_random,Protein_Maps_Ch1_random,'euclidean','Smallest',1); NND_Ch2_Ch1_Rnd=NND_Ch2_Ch1_Rnd';
     
     NND_Ch1_Rnd = pdist2(Protein_Maps_Ch1_random,Protein_Maps_Ch1_random,'euclidean','Smallest',2); NND_Ch1_Rnd=NND_Ch1_Rnd';NND_Ch1_Rnd=NND_Ch1_Rnd(:,2);
     NND_Ch2_Rnd = pdist2(Protein_Maps_Ch2_random,Protein_Maps_Ch2_random,'euclidean','Smallest',2); NND_Ch2_Rnd=NND_Ch2_Rnd';NND_Ch2_Rnd=NND_Ch2_Rnd(:,2);
    
     % Calculating NND for experimental data:

     NND_Ch1_Ch2 = pdist2(Protein_Maps_Ch1,Protein_Maps_Ch2,'euclidean','Smallest',1); NND_Ch1_Ch2=NND_Ch1_Ch2';
     NND_Ch2_Ch1 = pdist2(Protein_Maps_Ch2,Protein_Maps_Ch1,'euclidean','Smallest',1); NND_Ch2_Ch1=NND_Ch2_Ch1';
 
     NND_Ch1 = pdist2(Protein_Maps_Ch1,Protein_Maps_Ch1,'euclidean','Smallest',2); NND_Ch1=NND_Ch1';  NND_Ch1=NND_Ch1(:,2);
     NND_Ch2 = pdist2(Protein_Maps_Ch2,Protein_Maps_Ch2,'euclidean','Smallest',2); NND_Ch2=NND_Ch2'; NND_Ch2=NND_Ch2(:,2);
 
 
    % Calculatig mean and mediams for all the distributions
 
        % For Random distributions:
        Avg_NND_Ch1_Rnd = mean(NND_Ch1_Rnd); Median_NND_Ch1_Rnd = median(NND_Ch1_Rnd);
        Avg_NND_Ch2_Rnd = mean(NND_Ch2_Rnd); Median_NND_Ch2_Rnd = median(NND_Ch2_Rnd);

        Avg_NND_Ch1_Ch2_Rnd = mean(NND_Ch1_Ch2_Rnd); Median_NND_Ch1_Ch2_Rnd = median(NND_Ch1_Ch2_Rnd);
        Avg_NND_Ch2_Ch1_Rnd = mean(NND_Ch2_Ch1_Rnd); Median_NND_Ch2_Ch1_Rnd = median(NND_Ch2_Ch1_Rnd);
        
        % For experimental data:
        Avg_NND_Ch1 =mean(NND_Ch1); Median_NND_Ch1 = median(NND_Ch1);
        Avg_NND_Ch2= mean(NND_Ch2); Median_NND_Ch2 = median(NND_Ch2);
  
        Avg_NND_Ch1_Ch2 = mean(NND_Ch1_Ch2); Median_NND_Ch1_Ch2 = median(NND_Ch1_Ch2);
        Avg_NND_Ch2_Ch1 = mean(NND_Ch2_Ch1); Median_NND_Ch2_Ch1 = median(NND_Ch2_Ch1);
    
        % Calculating the % of distances below certain threshold (T1, T2)
 
        % For Random distributions:

        NND_T1_Ch1_Ch2_Rnd = sum(NND_Ch1_Ch2_Rnd<T1)/Protein_Number_Ch2*100;
        NND_T2_Ch1_Ch2_Rnd = sum(NND_Ch1_Ch2_Rnd<T2)/Protein_Number_Ch2*100;
 
        NND_T1_Ch2_Ch1_Rnd = sum(NND_Ch2_Ch1_Rnd<T1)/Protein_Number_Ch1*100;
        NND_T2_Ch2_Ch1_Rnd = sum(NND_Ch2_Ch1_Rnd<T2)/Protein_Number_Ch1*100;
        
        NND_T1_Ch2_Rnd = sum(NND_Ch2_Rnd<T1)/Protein_Number_Ch2 *100;
        NND_T2_Ch2_Rnd = sum(NND_Ch2_Rnd<T2)/Protein_Number_Ch2 *100;
 
        NND_T1_Ch1_Rnd = sum(NND_Ch1_Rnd<T1)/Protein_Number_Ch1*100;
        NND_T2_Ch1_Rnd = sum(NND_Ch1_Rnd<T2)/Protein_Number_Ch1*100;

 
        % For experimental data:

        NND_T1_Ch1_Ch2 = sum(NND_Ch1_Ch2<T1)/Protein_Number_Ch2*100;
        NND_T2_Ch1_Ch2 = sum(NND_Ch1_Ch2<T2)/Protein_Number_Ch2*100;
 
        NND_T1_Ch2_Ch1 = sum(NND_Ch2_Ch1<T1)/Protein_Number_Ch1*100;
        NND_T2_Ch2_Ch1 = sum(NND_Ch2_Ch1<T2)/Protein_Number_Ch1*100;
 
        NND_T1_Ch2 = sum(NND_Ch2<T1)/Protein_Number_Ch2*100;
        NND_T2_Ch2 = sum(NND_Ch2<T2)/Protein_Number_Ch2*100;
 
        NND_T1_Ch1 = sum(NND_Ch1<T1)/Protein_Number_Ch1*100;
        NND_T2_Ch1 = sum(NND_Ch1<T2)/Protein_Number_Ch1*100;
    
 % Mergind results into one matrix 
 
    Results_Ch1 = [Avg_NND_Ch1 Avg_NND_Ch1_Rnd Median_NND_Ch1 Median_NND_Ch1_Rnd NND_T1_Ch1 NND_T1_Ch1_Rnd NND_T2_Ch1 NND_T2_Ch1_Rnd];
    Results_Ch2 = [Avg_NND_Ch2 Avg_NND_Ch2_Rnd Median_NND_Ch2 Median_NND_Ch2_Rnd NND_T1_Ch2 NND_T1_Ch2_Rnd NND_T2_Ch2 NND_T2_Ch2_Rnd];
    Results_Ch1_Ch2 = [Avg_NND_Ch1_Ch2 Avg_NND_Ch1_Ch2_Rnd Median_NND_Ch1_Ch2 Median_NND_Ch1_Ch2_Rnd NND_T1_Ch1_Ch2 NND_T1_Ch1_Ch2_Rnd NND_T2_Ch1_Ch2 NND_T2_Ch1_Ch2_Rnd];
    Results_Ch2_Ch1 = [Avg_NND_Ch2_Ch1 Avg_NND_Ch2_Ch1_Rnd Median_NND_Ch2_Ch1 Median_NND_Ch2_Ch1_Rnd NND_T1_Ch2_Ch1 NND_T1_Ch2_Ch1_Rnd NND_T2_Ch2_Ch1 NND_T2_Ch2_Ch1_Rnd];

% Saving results:

    % NND for all the x,y proteins positions in data and random:
    dlmwrite(fullfile(Directory,'NND_Ch1_Ch2_Rnd.txt'),NND_Ch1_Ch2_Rnd, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch2_Ch1_Rnd.txt'),NND_Ch2_Ch1_Rnd, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch1_Rnd.txt'),NND_Ch1_Rnd,'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch2_Rnd.txt'),NND_Ch2_Rnd,'delimiter',',','precision',10, '-append');
    
    dlmwrite(fullfile(Directory,'NND_Ch1_Ch2.txt'),NND_Ch1_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch2_Ch1.txt'),NND_Ch2_Ch1, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch1.txt'),NND_Ch1,'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'NND_Ch2.txt'),NND_Ch2,'delimiter',',','precision',10, '-append');

    % Summary descriptors of distributions for Data and Random distributions:  
    dlmwrite(fullfile(Directory,'Summary_NND_Ch1.txt'),Results_Ch1, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'Summary_NND_Ch2.txt'),Results_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'Summary_NND_Ch1_Ch2.txt'),Results_Ch1_Ch2, 'delimiter',',','precision',10, '-append');
    dlmwrite(fullfile(Directory,'Summary_NND_Ch2_Ch1.txt'),Results_Ch2_Ch1, 'delimiter',',','precision',10, '-append');
    
    end 
end
 

