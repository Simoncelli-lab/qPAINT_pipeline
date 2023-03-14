
A = [a3 a1 a2 a4 a5 a6 dbscan_id dbscan_type]; 

preDBSCANfilename  = "YOUR_FILE_NAME";
clip = "_DBSCAN_EX_PX"; %please input the Eps (E) and MinPts (P) parameters you used
file_extension   = ".txt";

dlmwrite(strcat(preDBSCANfilename,clip,file_extension),A);
