
function f_background_estimation_by_subtrubct_SMs_small_FoV(locSheetFolder_for_tfrom,xy_ch_list_together,pixel_sz,dataN,fileFolder,inter_results_save_folder,ROI_centerY,W,Nimg,FoV,N_FoV,FoV_each,tform_FoV_each,gaussian_fit_size,PSF_sz_half,intensity_thred,crop_data_pre_name)
W = W/2;


for kk = 1:length(dataN)
 
% load list data

if xy_ch_list_together==0
fileName1 = strcat(locSheetFolder_for_tfrom,'data',num2str(dataN(kk)), '_xch', char(crop_data_pre_name), '.csv');
data = readtable(fileName1);

%
x_X = data.x_nm_./pixel_sz;
y_X = data.y_nm_./pixel_sz;
signal_X = data.intensity_photon_;
frameN_X = data.frame;
dataX_store = [frameN_X, x_X, y_X, signal_X];
dataX_store(dataX_store(:,4)<=intensity_thred,:) = [];

fileName1 = [[fileFolder,'BKG_list_data1\data'],  num2str(dataN(kk)), '_ych', char(crop_data_pre_name), '.csv'];
data = readtable(fileName1);

%
x_Y = data.x_nm_./pixel_sz;
y_Y = data.y_nm_./pixel_sz;
signal_Y = data.intensity_photon_;
frameN_Y = data.frame;

dataY_store = [frameN_Y, x_Y, y_Y, signal_Y];
dataY_store(dataY_store(:,4)<=intensity_thred,:) = [];

else 

fileName1 = strcat(locSheetFolder_for_tfrom,'data',num2str(dataN(kk)), '_xych', char(crop_data_pre_name), '.csv');
data = readtable(fileName1);

%
x_XY = data.x_nm_./pixel_sz;
y_XY = data.y_nm_./pixel_sz;
signal_XY = data.intensity_photon_;
frameN_XY = data.frame;
dataXY_store = [frameN_XY, x_XY, y_XY, signal_XY];


indX = x_XY>W;
dataX_store = dataXY_store(indX,:);  dataX_store(:,2) = dataX_store(:,2)-W;
dataY_store = dataXY_store(~indX,:);  dataY_store(:,2) = W-dataY_store(:,2);

dataX_store(dataX_store(:,4)<=intensity_thred,:) = [];
dataY_store(dataY_store(:,4)<=intensity_thred,:) = [];
end

%% load the crop data


center_x = FoV(1)/N_FoV(1)/2*[-N_FoV(1)+1:2:N_FoV(1)-1];
center_y = FoV(2)/N_FoV(2)/2*[-N_FoV(2)+1:2:N_FoV(2)-1];
[center_X,center_Y] = meshgrid(center_x,center_y);
center_X = center_X(:);
center_Y = center_Y(:);



for ii = 1:length(center_X)

%caculate the current center

ROI_centerY_cur = ROI_centerY+[center_x(ii),center_y(ii)];
ROI_centerY_cur_old = ROI_centerY_cur;
ROI_centerY_cur = ROI_centerY_cur+0.5; 

tformLoad = load(strcat(fileFolder, '\tformx2y_y_center_', num2str(ROI_centerY_cur_old(1)),'_', num2str(ROI_centerY_cur_old(2)), '_FoV_', num2str(tform_FoV_each),'.mat'));
tformx2y = tformLoad.tformx2y;
ROI_centerX_cur = round(transformPointsInverse(tformx2y,[W,0]+[-ROI_centerY_cur(1),ROI_centerY_cur(2)])+[W,0]);
ROI_centerY_cur = round(ROI_centerY_cur+0.5);  ROI_centerX_cur = round(ROI_centerX_cur+0.5);  % +0.5 is used for compensate the coodinate difference between thunderstorm and matlab

ROI_centerX_cur = ROI_centerX_cur-[W,0]; ROI_centerY_cur(1) = W-ROI_centerY_cur(1);
ROI_centerX_lu_corner = ROI_centerX_cur-FoV_each/2;
ROI_centerY_lu_corner = ROI_centerY_cur-FoV_each/2;

% loc_list in the FoV
dataX = dataX_store; dataY = dataY_store;
dataX(:,2:3) = dataX(:,2:3)-ROI_centerX_lu_corner; dataY(:,2:3) = dataY(:,2:3)-ROI_centerY_lu_corner;
indxX = abs(dataX(:,2)-FoV_each/2)<(FoV_each/2+10) &   abs(dataX(:,3)-FoV_each/2)<(FoV_each/2+10) ;  % +10 is including some edge emitters
indxY = abs(dataY(:,2)-FoV_each/2)<(FoV_each/2+10) &   abs(dataY(:,3)-FoV_each/2)<(FoV_each/2+10) ;  % +10 is including some edge emitters
dataX = dataX(indxX,:); dataY = dataY(indxY,:);

frameN_X = dataX(:,1); x_X = dataX(:,2); y_X=dataX(:,3); signal_X=dataX(:,4);
frameN_Y = dataY(:,1); x_Y = dataY(:,2); y_Y=dataY(:,3); signal_Y=dataY(:,4);
%%

SMLMName = ['data',num2str(dataN(kk)),'_centerY_y',num2str(ROI_centerY(1)),'_x_',num2str(ROI_centerY(2)),'_','FoV',num2str(FoV(1)),'_',num2str(FoV(2)),'_',num2str(ii),'th_FoV','.tif'];
offset_Name= ['data_offset_centerY_y',num2str(ROI_centerY(1)),'_x_',num2str(ROI_centerY(2)),'_','FoV',num2str(FoV(1)),'_',num2str(FoV(2)),'_',num2str(ii),'th_FoV','.mat'];

offsetLoad = load(strcat(fileFolder,offset_Name));


SM_img_offset = offsetLoad.offset;
SMLMR = Tiff(strcat(fileFolder,SMLMName),'r');
for i=1:Nimg
    setDirectory(SMLMR,i);
    SM_img(:,:,i) = double(SMLMR.read);
end
imgSzx = size(SM_img,2)/2;
imgSzy = size(SM_img,1);


SM_img = SM_img-SM_img_offset;

%%
SM_img_subtractX = SM_img(:,(1:imgSzx),:);
SM_img_subtractY = SM_img(:,(1:imgSzx)+imgSzx,:); 

%%
count = 0; firstImage=1;
SM_img_subtractX_old = SM_img_subtractX;
for i = 1:Nimg


    for j = 1:sum(frameN_X==i)
        count = count+1;
        rangex = [max(1,-PSF_sz_half+round(x_X(count))):min(PSF_sz_half+round(x_X(count)),imgSzx)];
        rangey = [max(1,-PSF_sz_half+round(y_X(count))):min(PSF_sz_half+round(y_X(count)),imgSzy)];
        SM_img_subtractX(rangey,rangex,i) = nan; 
    end

    if firstImage==1
        Fig1 = figure('Visible','off','Units','inches','InnerPosition',[1,1,10,8]);
        subplot(1,2,1);
        imagesc(SM_img_subtractX_old(:,:,i)); axis image;  colorbar; hold on;
        scatter(x_X(frameN_X==i),y_X(frameN_X==i),10,'r+');
        subplot(1,2,2);
        imagesc(SM_img_subtractX(:,:,i)); axis image; colorbar; sgtitle('background est exampleFrameX');
        exportgraphics(Fig1,strcat(inter_results_save_folder,'\',num2str(dataN(kk)),'background_est_exampleFrameX.jpg'),'Resolution',600);
        firstImage=0;
    end
end
%%
count = 0; firstImage=1;
SM_img_subtractY_old = SM_img_subtractY;
for i = 1:Nimg
    
    for j = 1:sum(frameN_Y==i)
        count = count+1;
        rangex = [max(1,-PSF_sz_half+round(x_Y(count))):min(PSF_sz_half+round(x_Y(count)),imgSzx)];
        rangey = [max(1,-PSF_sz_half+round(y_Y(count))):min(PSF_sz_half+round(y_Y(count)),imgSzy)];
        SM_img_subtractY(rangey,rangex,i) = nan; 
       
    end
    

        if firstImage==1
        Fig1 = figure('Visible','off','Units','inches','InnerPosition',[1,1,10,8]);
        subplot(1,2,1);
        imagesc(SM_img_subtractY_old(:,:,i)); axis image;  colorbar; hold on;
        scatter(x_Y(frameN_Y==i),y_Y(frameN_Y==i),10,'r+');
        subplot(1,2,2);
        imagesc(SM_img_subtractY(:,:,i)); axis image; colorbar; sgtitle('background est exampleFrameY');
        exportgraphics(Fig1,strcat(inter_results_save_folder,'\',num2str(dataN(kk)),'background_est_exampleFrameY.jpg'),'Resolution',600);
        firstImage=0;
        end

end
%% background in X channel
% save bkg frame every 50 frames

count = 0;
[HImage,WImage,L]=size(SM_img_subtractY);
center_size = 10;
range = round([-center_size/2,center_size/2]);
for i=1:50:Nimg
    count = count+1;
    indx_start = max(1,i-50);
    indx_end = min(i+50,Nimg);
    back_cur = nanmean(SM_img_subtractX(:,:,indx_start:indx_end),3);
    back_cur(back_cur<0)=0;
    %back_cur(isnan(back_cur))=nanmean(back_cur(round(HImage/2)+range,round(WImage/2)+range),'all');
    back_cur(isnan(back_cur))=nanmax(back_cur,[],'all');
    SMLM_bkgX(:,:,count) = imgaussfilt(back_cur(:,:,:),gaussian_fit_size);
    

end
%clear SM_img_subtractX
%% background in Y channel
count = 0;
for i=1:50:Nimg
    count = count+1;
    indx_start = max(1,i-50);
    indx_end = min(i+50,Nimg);
    back_cur = nanmean(SM_img_subtractY(:,:,indx_start:indx_end),3);
    back_cur(back_cur<0)=0;
    %back_cur(isnan(back_cur))=nanmean(back_cur(round(HImage/2)+range,round(WImage/2)+range),'all');
    back_cur(isnan(back_cur))=nanmax(back_cur,[],'all');
    SMLM_bkgY(:,:,count) = imgaussfilt(back_cur(:,:,:),gaussian_fit_size);

end
%clear SM_img_subtractY
SMLM_bkg = [SMLM_bkgX,SMLM_bkgY];

Fig1 = figure('Visible','off','Units','inches','InnerPosition',[1,1,10,8]);
imagesc(SMLM_bkg(:,:,1)); axis image;  colorbar; hold on; title('estimated_backgraound');
exportgraphics(Fig1,strcat(inter_results_save_folder,'\',num2str(dataN(kk)),'estimated_background.jpg'),'Resolution',600);

%% crop FoV; copy from the the crop_save_image.m code

SMLM_save_Nmae = ['data',num2str(dataN(kk)),'_bkg_centerY_y',num2str(ROI_centerY(1)),'_x_',num2str(ROI_centerY(2)),'_','FoV',num2str(FoV(1)),'_',num2str(FoV(2)),'_1th_FoV','.mat'];
save(strcat(fileFolder,SMLM_save_Nmae),'SMLM_bkg')


end
end


end


