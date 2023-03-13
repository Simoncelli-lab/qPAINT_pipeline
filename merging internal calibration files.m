%% Concatenating files

clear all
%% merging the internal calibration files

ROI1 =  importdata('Day_2_C1_pMHC+ICAM_pZAP70_FOV1_ROI1_Internal_calibration_200nm.txt');  
ROI2 =  importdata('Day_2_C3_aCD3+aCD28_pZAP70_FOV1_ROI1_Internal_calibration_200nm.txt');  
ROI3 =  importdata('Day_2_C3_aCD3+aCD28_pZAP70_FOV2_ROI1_Internal_calibration_200nm.txt');  
ROI4 =  importdata('Day_3_C1_pMHC+ICAM_pZAP70_FOV3_ROI1_Internal_calibration_200nm.txt');  
ROI5 =  importdata('Day_3_C1_pMHC+ICAM_pZAP70_FOV3_ROI2_Internal_calibration_200nm.txt');  
ROI6 =  importdata('Day_3_C1_pMHC+ICAM_pZAP70_FOV4_ROI1_Internal_calibration_200nm.txt');  
ROI7 =  importdata('Day_3_C3_aCD3+aCD28_pZAP70_FOV1_ROI1_Internal_calibration_200nm.txt');  
ROI8 =  importdata('Day_3_C3_aCD3+aCD28_pZAP70_FOV1_ROI2_Internal_calibration_200nm.txt');  
ROI9 =  importdata('Day_3_C3_aCD3+aCD28_pZAP70_FOV3_ROI1_Internal_calibration_200nm.txt');  
%ROI10 =  importdata('chamber_2_CD3CD28_FOV_2_pZAP70_ROI4_Internal_calibration_200nm.txt');  
%ROI11 =  importdata('chamber_2_CD3CD28_FOV_2_pZAP70_ROI5_Internal_calibration_200nm.txt');  
%ROI12 =  importdata('chamber_3_CD3_FOV_1_pZAP70_ROI1_Internal_calibration_200nm.txt');  
%ROI13 =  importdata('chamber_3_CD3_FOV_1_pZAP70_ROI2_Internal_calibration_200nm.txt');  
%ROI14 =  importdata('chamber_3_CD3_FOV_1_pZAP70_ROI3_Internal_calibration_200nm.txt');  
%ROI15 =  importdata('chamber_3_CD3_FOV_1_pZAP70_ROI4_Internal_calibration_200nm.txt');  
%ROI16 =  importdata('chamber_3_CD3_FOV_1_pZAP70_ROI5_Internal_calibration_200nm.txt');  
%ROI17 =  importdata('chamber_4_CD3_FOV_1_pZAP70_ROI1_Internal_calibration_200nm.txt');  
%ROI18 =  importdata('chamber_4_CD3_FOV_1_pZAP70_ROI2_Internal_calibration_200nm.txt');  
%ROI19 =  importdata('chamber_4_CD3_FOV_1_pZAP70_ROI3_Internal_calibration_200nm.txt');  


All_ROIs = [ROI1.data' ROI2.data' ROI3.data' ROI4.data' ROI5.data' ROI6.data' ROI7.data' ROI8.data' ROI9.data'];All_ROIs = All_ROIs';
dlmwrite('All_ROIS_pZAP70_internal_calibration_200nm.txt',All_ROIs);

%% merging the internal calibration files

ROI1 =  importdata('Day_1_C1_pMHC+ICAM_pCD3_FOV1_ROI1_Internal_calibration_200nm.txt');  
ROI2 =  importdata('Day_1_C1_pMHC+ICAM_pCD3_FOV2_ROI1_Internal_calibration_200nm.txt');  
ROI3 =  importdata('Day_1_C1_pMHC+ICAM_pCD3_FOV3_ROI1_Internal_calibration_200nm.txt');  

All_ROIs = [ROI1.data' ROI2.data' ROI3.data' ];All_ROIs = All_ROIs';
%dlmwrite('All_ROIS_pCD3_internal_calibration_200nm.txt',All_ROIs);
