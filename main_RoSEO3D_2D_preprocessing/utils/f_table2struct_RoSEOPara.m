function para = f_table2struct_RoSEOPara(sheet_cur)

para = struct(...
    'emitter_wavelength',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_emitter_wavelength',2)),...
    'refractiveIndx',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_refractiveIndx',2)),...
    'sampleRefractiveIndx',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_sampleRefractiveIndx',2)),...
    'left_to_right_trans_ratio',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_left_to_right_trans_ratio',2)),...
    'calibrated_mask',(sheet_cur(sheet_cur(:,1)=='ROSEO_calibrated_mask',2)),...
    'maskName',(sheet_cur(sheet_cur(:,1)=='ROSEO_maskName',2)),...
    'toPhoton',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_toPhoton',2)),...
    'reg',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_reg',2)),...
    'pix_sizex',str2num(sheet_cur(sheet_cur(:,1)=='ROSEO_pix_sizex',2)),...
    'FoV_analysed',str2num(sheet_cur(sheet_cur(:,1)=='FoV_analysed',2)),...
    'dataID_analysed',str2num(sheet_cur(sheet_cur(:,1)=='dataID_analysed',2)),...
    'Nframe_analysed',str2num(sheet_cur(sheet_cur(:,1)=='Nframe_analysed',2)),...
    'estResults_save_folder',(sheet_cur(sheet_cur(:,1)=='estResults_save_folder',2)));
    
end