%% Converting hdf5 file to txt and cutting ROIs to specif size (i.e. 2X2, 3x3 or 5x5).
clear all

%this section of code converts hdf5 files to .txt format

filename = 'C1_pMHC+ICAM_pCD3_FOV1_ROI1.hdf5';

info = h5info(filename);
data = h5read(filename,'/locs');
pixelsize = 130;
frame = data.frame;
xloc = data.x; xloc = xloc.*pixelsize;
yloc = data.y; yloc = yloc.*pixelsize;
sx = data.lpx; sx = sx.*pixelsize;
sy = data.lpy; sy = sy.*pixelsize;
sd = (sx + sy)./2; sd=sd.*pixelsize;
photons = data.photons;

% Cutting regions of interest of 3 x 3= micrometers

%xmin = min(xloc); xloc_rel = xloc - xmin;
%ymin = min(yloc); yloc_rel = yloc - ymin;

listLocalizations = [xloc yloc frame photons sx sy];

%rowsToDelete_x = xloc_rel > 3000; % this number relates to the region size (i.e 3 micrometer)
%listLocalizations(rowsToDelete_x,:) = [];

%rowsToDelete_y = listLocalizations(:,2) > 3000; % this number relates to the region size (i.e 3 micrometer)
%listLocalizations(rowsToDelete_y,:) = [];

dlmwrite('C1_pMHC+ICAM_pCD3_FOV1_ROI1.txt',listLocalizations);
