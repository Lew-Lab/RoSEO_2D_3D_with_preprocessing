
%% set the data want to analyse
onWorkstation = 0;

sheetFolder = 'D:\Experimental_data\20230523 Poly-A poly-C SYTOX NR optimize_fromYXQ\processed data RoSEO\RoSEO_dataSheet\';   %E:\Experimental_data
%sheetFolder = 'D:\Box Sync\data\20210504 3D lipid membrane\processed data RoSEO\RoSEO_dataSheet\';
sheetName = {'20230524_data56_57_FoV1.xlsx'};

% sheetFolder = 'E:\Experimental_data\20220924_A1LCD\processed data RoSEO\RoSEO_dataSheet\';
% sheetName = {'20220924_RoIsheet_data2_6_FoV1_MC540.xlsx'
%     '20220924_RoIsheet_data15_19_FoV1_NB.xlsx'
% };

%E:\Experimental_data\20220814 A1-LCD\processed data RoSEO\RoSEO_dataSheet\
%'20220814_RoIsheet_data19_23_FoV1.xlsx'
%     '20220814_RoIsheet_data19_23_FoV2.xlsx',...
%     '20220814_RoIsheet_data19_23_FoV3.xlsx'
% '20220814_RoIsheet_data38_42_FoV1.xlsx'

%'E:\Experimental_data\20220917 A1-LCD\processed data RoSEO\RoSEO_dataSheet\'
%'20220917_RoIsheet_data110_120_FoV1.xlsx'
%'20220917_RoIsheet_data141_147_FoV1.xlsx'
%'20220917_RoIsheet_data141_147_FoV2.xlsx'
%'20220917_RoIsheet_data141_147_FoV3.xlsx'
%'20220917_RoIsheet_data129_133_FoV1.xlsx',...
% '20220917_RoIsheet_data129_133_FoV2.xlsx',...
% '20220917_RoIsheet_data129_133_FoV3.xlsx',...
% '20220917_RoIsheet_data129_133_FoV4.xlsx',...

%E:\Experimental_data\20220924_A1LCD\processed data RoSEO\RoSEO_dataSheet\
%     '20220924_RoIsheet_data29_34_FoV1.xlsx',...
%     '20220924_RoIsheet_data29_34_FoV2.xlsx',...
%     '20220924_RoIsheet_data29_34_FoV3.xlsx',...
%     '20220924_RoIsheet_data29_34_FoV4.xlsx'
% '20220924_RoIsheet_data36_40_FoV1.xlsx',...
%     '20220924_RoIsheet_data36_40_FoV2.xlsx',...
%     '20220924_RoIsheet_data36_40_FoV3.xlsx',...
%     '20220924_RoIsheet_data36_40_FoV1_BF.xlsx',...
%      '20220924_RoIsheet_data36_40_FoV2_BF.xlsx',...
%     '20220924_RoIsheet_data36_40_FoV3_BF.xlsx'
%'20220924_RoIsheet_data2_6_FoV1_MC540.xlsx'


% sheetFolder = 'E:\Experimental_data\20221123 A1-LCD\processed data RoSEO\RoSEO_dataSheet\';
% sheetName = {
%'20221123_RoIsheet_data92_105_FoV1_NB.xlsx',...
% '20221123_RoIsheet_data92_105_FoV1_NR.xlsx'
%'20221123_RoIsheet_data92_105_FoV2_NB.xlsx',...
% '20221123_RoIsheet_data92_105_FoV2_NR.xlsx'
% '20221123_RoIsheet_data92_105_FoV3_NB.xlsx',...
% '20221123_RoIsheet_data92_105_FoV3_NR.xlsx'
%'20221123_RoIsheet_data42_54_FoV1_NB.xlsx',...
% '20221123_RoIsheet_data42_54_FoV1_NR.xlsx',...
% '20221123_RoIsheet_data42_54_FoV2_NB.xlsx',...
% '20221123_RoIsheet_data42_54_FoV2_NR.xlsx',...
% '20221123_RoIsheet_data42_54_FoV3_NB.xlsx',...
% '20221123_RoIsheet_data42_54_FoV3_NR.xlsx'
% '20221123_RoIsheet_data20_40_FoV1_NB.xlsx',...
% '20221123_RoIsheet_data20_40_FoV1_NR.xlsx',...
% '20221123_RoIsheet_data20_40_FoV2_NB.xlsx',...
% '20221123_RoIsheet_data20_40_FoV2_NR.xlsx',...
%'20221123_RoIsheet_data20_40_FoV3_NB.xlsx',...
%'20221123_RoIsheet_data20_40_FoV3_NR.xlsx',...
% };

% 'E:\Experimental_data\20221201 A1-LCD\processed data RoSEO\RoSEO_dataSheet\'
%'20221201_RoIsheet_data13_29_FoV1_NB.xlsx',...
% '20221201_RoIsheet_data13_29_FoV1_NR.xlsx',...
% '20221201_RoIsheet_data13_29_FoV2_NB.xlsx',...
% '20221201_RoIsheet_data13_29_FoV2_NR.xlsx',...
% '20221201_RoIsheet_data13_29_FoV3_NB.xlsx',...
% '20221201_RoIsheet_data13_29_FoV3_NR.xlsx',...
% '20221201_RoIsheet_data13_29_FoV4_NB.xlsx',...
% '20221201_RoIsheet_data13_29_FoV4_NR.xlsx',...
%'20221201_RoIsheet_data30_41_FoV1_NB',...
%'20221201_RoIsheet_data30_41_FoV1_NR',...
%'20221201_RoIsheet_data30_41_FoV2_NB',...
%'20221201_RoIsheet_data30_41_FoV2_NR'


% 'E:\Experimental_data\20221006 A1-LCD\processed data RoSEO\RoSEO_dataSheet\';
%     '20221006_RoIsheet_data2_6_FoV1_MC540.xlsx'
%     '20221006_RoIsheet_data8_12_FoV1_NB.xlsx'

%%
for ii = 1:numel(sheetName)
%% read the sheet info and create necessary folder
opts = detectImportOptions([sheetFolder,sheetName{ii}]);
opts = setvartype(opts,'char');  % or 'string'
sheet_cur= readtable([sheetFolder,sheetName{ii}],opts);
if onWorkstation == 0
sheet_cur = string(table2array(sheet_cur(:,2:3)));
else 
sheet_cur = string(table2array(sheet_cur(:,[2,4])));
end
para_data = f_table2struct(sheet_cur);

%

% check if the save fileExit, create new one
time = datestr(now,'yyyymmdd_HH_MM');
if isfolder(strcat(para_data.saveData_folderName))==0
    mkdir(strcat(para_data.saveData_folderName));
end

if isfolder(strcat(para_data.SMLM_image_save_folder))==0
    mkdir(strcat(para_data.SMLM_image_save_folder));
end



if sheet_cur(sheet_cur(:,1)=='inter_results_save_folder',2)==""
    inter_results_save_folder = strcat(para_data.saveData_folderName,string(time),'results_save');
    mkdir(inter_results_save_folder);
else
   inter_results_save_folder = sheet_cur(sheet_cur(:,1)=='inter_results_save_folder',2);
end
save(strcat(inter_results_save_folder,'/sheet_cur.mat'),'para_data','sheet_cur');
%% create tform using the imaging sample
if sheet_cur(sheet_cur(:,1)=='do_tform',2)=="1"
pixel_sz = 1; % thurderstorm camera pixel size setting
photonThred = 100;
frame_thred = 20000;
ratio_y2x = 1.14;
filter_close_emitter = 0;
xy_ch_list_together = 1;                                                                                                                                                                                                                                                                                                                                                                                                 

f_tform_use_sample(para_data.locSheetFolder_for_tfrom,...
    para_data.frame_width,...
    pixel_sz,...
    photonThred,...
    frame_thred,...
    ratio_y2x,...
    para_data.tform_ROI_center,...
    para_data.tform_ROI_size_each,...
    para_data.tform_ROI_size_whole,...
    para_data.tform_data,...
    para_data.tform_N_Fov,...
    filter_close_emitter,...
    inter_results_save_folder, ...
    xy_ch_list_together,...
    para_data.SMLM_image_save_folder,...
    para_data.locSheetFolder_for_tfrom_beads,...
    para_data.tform_data_beads,...
    para_data.tform_data_pre_name);
end
%% crop image using the above tform
if sheet_cur(sheet_cur(:,1)=='do_crop',2)=="1"
f_crop_image_usig_whole_pixels(...
    para_data.crop_data,...
    para_data.DatafileFolder,...
    inter_results_save_folder,...
    para_data.SMLM_image_save_folder,...
    para_data.crop_ROI_center,...
    para_data.frame_width,...
    para_data.crop_N_frame,...
    para_data.crop_ROI_all,...
    para_data.crop_N_FoV,...
    para_data.crop_ROI_each,...
    para_data.tform_ROI_size_each, ...
    para_data.crop_data_pre_name);
end
%% crop offset image
if sheet_cur(sheet_cur(:,1)=='do_crop',2)=="1"
f_crop_offset_using_whole_pixels(...
    para_data.offSetName,...
    inter_results_save_folder,...
    para_data.SMLM_image_save_folder,...
    para_data.crop_ROI_center,...
    para_data.frame_width,...
    para_data.crop_ROI_all,...
    para_data.crop_N_FoV,...
    para_data.crop_ROI_each,...
    para_data.tform_ROI_size_each);

end
%% estimate background
if sheet_cur(sheet_cur(:,1)=='do_backgroundEst',2)=="1"
PSF_sz_half = 7;
gaussian_fit_size = 5;
pixel_sz = 1;
xy_ch_list_together = 1;
intensity_thred = 0;

f_background_estimation_by_subtrubct_SMs_small_FoV(...
    para_data.locSheetFolder_for_tfrom,...
    xy_ch_list_together,...
    pixel_sz,...
    para_data.crop_data,...
    para_data.SMLM_image_save_folder,...
    inter_results_save_folder,...
    para_data.crop_ROI_center,...
    para_data.frame_width,...
    para_data.crop_N_frame,...
    para_data.crop_ROI_all,...
    para_data.crop_N_FoV,...
    para_data.crop_ROI_each,...
    para_data.tform_ROI_size_each,...
    gaussian_fit_size,...
    PSF_sz_half,...
    intensity_thred, ...
    para_data.crop_data_pre_name)

end

close all;

%% estimate using RoSEO

addpath(genpath('F:\OneDrive - Washington University in St. Louis\github\RoSEO_in_use\RoSE-O-master_download20210115'));
RoSEOPara = f_table2struct_RoSEOPara(sheet_cur);

if sheet_cur(sheet_cur(:,1)=='do_RoSEO',2)=="1"

if isfolder(strcat(RoSEOPara.estResults_save_folder))==0
    mkdir(strcat(RoSEOPara.estResults_save_folder));
end
thred = 30;

f_RoSEO_input(...
    RoSEOPara,...
    para_data.SMLM_image_save_folder,...
    inter_results_save_folder,...
    RoSEOPara.dataID_analysed,...
    RoSEOPara.Nframe_analysed,...
    para_data.crop_ROI_each,...
    RoSEOPara.FoV_analysed,...
    para_data.crop_ROI_center,...
    RoSEOPara.estResults_save_folder,...
    thred)

end
%% reconstruaction for validation


%%

close all;
end

