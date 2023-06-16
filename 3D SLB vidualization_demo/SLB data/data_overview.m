

classdef data_overview
    methods(Static)


%% DPPC+cholesterol data
function [locFileName,imageName,angleEstName,center_SLB,radius_SLB,dataIDN,frameN_per_stack,bkgName, offsetName] = DPPC_chol_dataSource(dataN)
switch dataN


    case 4
%---------------#3------------------
fileNameAll = 'F:\OneDrive - Washington University in St. Louis\github\Deep-SMOLM3D\3D SLB vidualization\SLB data\';
locFileName = [fileNameAll '20230208_SLB\data26_31_FoV1\locList_FoV1.mat'];
angleEstName = [fileNameAll '20230208_SLB\data26_31_FoV1\angleList_FoV1.mat'];

imageName = [fileNameAll '20230208_SLB\data26_31_FoV1\data26_centerY_y389_x_163_FoV91_91_1th_FoV.tif'];
bkgName = [fileNameAll '20230208_SLB\data26_31_FoV1\data26_bkg_centerY_y389_x_163_FoV91_91_1th_FoV.mat'];
offsetName = [fileNameAll '20230208_SLB\data26_31_FoV1\data_offset_centerY_y389_x_163_FoV91_91_1th_FoV.mat'];


dataIDN = [26:30];
center_SLB = [0,0,1000];
radius_SLB = 1000;
frameN_per_stack = 2000;


    


end


    end
    end
end