A = [frame x y photons sdx sdy dbscan_id dbscan_type]; 
 
preDBSCANfilename  = "data_channel1";
clip = "_DBSCAN_10_10"; %please input the Eps (E) and MinPts (P) parameters you used
file_extension   = ".txt";
 
dlmwrite(strcat(preDBSCANfilename,clip,file_extension),A);
