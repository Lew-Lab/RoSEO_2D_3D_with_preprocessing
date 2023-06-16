
function f_crop_offset_using_whole_pixels(offsetFolder,inter_results_save_folder,saveData_folder_Name,ROI_centerY,W,FoV,N_FoV,FoV_each,tform_FoV_each)
W  = W/2;

%%
load(offsetFolder); offset_wholeFoV = offset;
H = size(offset,1);
%% define the cropping info


Nimg = 1;


center_x = FoV(1)/N_FoV(1)/2*[-N_FoV(1)+1:2:N_FoV(1)-1];
center_y = FoV(2)/N_FoV(2)/2*[-N_FoV(2)+1:2:N_FoV(2)-1];
[center_X,center_Y] = meshgrid(center_x,center_y);
center_X = center_X(:);
center_Y = center_Y(:);



for ii = 1:length(center_X)

%count = count+1;
range = round(-(FoV_each-1)/2):1:round((FoV_each-1)/2);

SMLM_save_Nmae = ['data_offset_centerY_y',num2str(ROI_centerY(1)),'_x_',num2str(ROI_centerY(2)),'_','FoV',num2str(FoV(1)),'_',num2str(FoV(2)),'_',num2str(ii),'th_FoV','.mat'];



ROI_centerY_cur = ROI_centerY+[center_x(ii),center_y(ii)];
ROI_centerY_cur_old = ROI_centerY_cur;
ROI_centerY_cur = ROI_centerY_cur+0.5;  % as in thunderstorm [N.5,N.5] is the center of a pixel, so I want to use center of a pixel for registration
ROI_Y_all_cur = [ROI_centerY_cur(1)-range.',ROI_centerY_cur(2)+range.'];
[ROI_Y_all_curX,ROI_Y_all_curY] = meshgrid(ROI_centerY_cur(1)-range.',ROI_centerY_cur(2)+range.');
ROI_Y_all_curX = reshape(ROI_Y_all_curX,[],1);
ROI_Y_all_curY = reshape(ROI_Y_all_curY,[],1);


load(strcat(saveData_folder_Name, '\tformx2y_y_center_', num2str(ROI_centerY_cur_old(1)),'_', num2str(ROI_centerY_cur_old(2)), '_FoV_', num2str(tform_FoV_each),'.mat'));
ROI_centerX_cur = round(transformPointsInverse(tformx2y,[W,0]+[-ROI_centerY_cur(1),ROI_centerY_cur(2)])+[W,0]);
ROI_X_all_cur = (transformPointsInverse(tformx2y,[W,0]+[-ROI_Y_all_curX,ROI_Y_all_curY])+[W,0]);

ROI_centerY_cur = round(ROI_centerY_cur+0.5);  ROI_centerX_cur = round(ROI_centerX_cur+0.5);  % +0.5 is used for compensate the coodinate difference between thunderstorm and matlab
ROI_X_all_cur = ROI_X_all_cur+0.5;



ROI_X_all_curX = ROI_X_all_cur(:,1);
ROI_X_all_curY =ROI_X_all_cur(:,2);
ROI_X_all_cur_intX = round(ROI_X_all_curX);
ROI_X_all_cur_intY = round(ROI_X_all_curY);
linear_idx = sub2ind([H,W*2],ROI_X_all_cur_intY,ROI_X_all_cur_intX);
% 


SMLM_img_ROIy = offset_wholeFoV(ROI_centerY_cur(2)+range,ROI_centerY_cur(1)-range); % flip operator is in -
SMLM_img_ROIx = offset_wholeFoV(linear_idx);
SMLM_img_ROIx = reshape(SMLM_img_ROIx ,FoV_each,FoV_each);
offset = [SMLM_img_ROIx,SMLM_img_ROIy];
    

save(strcat(saveData_folder_Name,SMLM_save_Nmae),'offset');

end
end
