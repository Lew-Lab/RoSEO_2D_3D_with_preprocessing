

fileFolder = ['E:\Experimental_data\20220814 A1-LCD\'];
offsetFolder = [fileFolder,'\processed data RoSEO\offSet.mat'];
load([fileFolder,'\processed data RoSEO\data38-42\tformx2y_y_center_351_162_FoV_170.mat']);
load(offsetFolder); offset_wholeFoV = offset;
tform = tformx2y;
SMLM_save_Nmae_pre = [fileFolder 'processed data RoSEO\data38-42\'];


%% define the cropping info
ROI_centerY = [351,162];
W = 1748/2;
ROI_centerX = transformPointsInverse(tformx2y,[W,0]+[-ROI_centerY(1),ROI_centerY(2)])+[W,0];

Nimg = 1;
FoV = [81,81];  %the size of whole FoV
N_FoV = [1,1]; %seperate the whole FoV into # of sub-FoV
FoV_each = 81;  % the size of each sub-FoV

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
ROI_centerY_cur = ROI_centerY_cur+0.5;  % as in thunderstorm [N.5,N.5] is the center of a pixel, so I want to use center of a pixel for registration
ROI_Y_all_cur = [ROI_centerY_cur(1)-range.',ROI_centerY_cur(2)+range.'];
[ROI_Y_all_curX,ROI_Y_all_curY] = meshgrid(ROI_centerY_cur(1)-range.',ROI_centerY_cur(2)+range.');
ROI_Y_all_curX = reshape(ROI_Y_all_curX,[],1);
ROI_Y_all_curY = reshape(ROI_Y_all_curY,[],1);


% if distance1<distance2
%     tformx2y=tform1;
% else
%     tformx2y=tform2;
% end
%load(strcat('E:\Experimental_data\20220530 amyloid fibril\','\processed data8 data9_16\saved_beads_loc_for_tform_v3\tformx2y_y_center_',num2str(ROI_centerY_cur(1)),'_',num2str(ROI_centerY_cur(2)),'_FoV_150.mat'));
ROI_centerX_cur = round(transformPointsInverse(tformx2y,[W,0]+[-ROI_centerY_cur(1),ROI_centerY_cur(2)])+[W,0]);
ROI_X_all_cur = (transformPointsInverse(tformx2y,[W,0]+[-ROI_Y_all_curX,ROI_Y_all_curY])+[W,0]);

ROI_centerY_cur = round(ROI_centerY_cur+0.5);  ROI_centerX_cur = round(ROI_centerX_cur+0.5);  % +0.5 is used for compensate the coodinate difference between thunderstorm and matlab
ROI_X_all_cur = ROI_X_all_cur+0.5;



% ROI_centerX_cur = round(transformPointsInverse(tformx2y,[W,0]+[-ROI_centerY_cur(1),ROI_centerY_cur(2)])+[W,0]);
% ROI_X_all_cur = (transformPointsInverse(tformx2y,[W,0]+[-ROI_Y_all_cur(:,1),ROI_Y_all_cur(:,2)])+[W,0]);
% [ROI_X_all_curX,ROI_X_all_curY] = meshgrid(ROI_X_all_cur(end:-1:1,1),ROI_X_all_cur(:,2));
%ROI_X_all_curX = reshape(ROI_X_all_cur(:,1),FoV_each,FoV_each);
%ROI_X_all_curY = reshape(ROI_X_all_cur(:,2),FoV_each,FoV_each);
ROI_X_all_curX = ROI_X_all_cur(:,1);
ROI_X_all_curY =ROI_X_all_cur(:,2);
ROI_X_all_cur_intX = round(ROI_X_all_curX);
ROI_X_all_cur_intY = round(ROI_X_all_curY);
linear_idx = sub2ind([400,W*2],ROI_X_all_cur_intY,ROI_X_all_cur_intX);
% 




SMLM_img_ROIy = offset_wholeFoV(ROI_centerY_cur(2)+range,ROI_centerY_cur(1)-range); % flip operator is in -
SMLM_img_ROIx = offset_wholeFoV(linear_idx);
SMLM_img_ROIx = reshape(SMLM_img_ROIx ,FoV_each,FoV_each);
offset = [SMLM_img_ROIx,SMLM_img_ROIy];
    

save([SMLM_save_Nmae_pre,SMLM_save_Nmae],'offset');

end
