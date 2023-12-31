function [SMLM_img,BKG_img] = img_pre_processing(options,indx_stack)

pixSize = options.pixSize;
registershifty = options.registershifty;
registershiftx = options.registershiftx;
imgSize = options.imgSize;
imgCenterY = options.imgCenter;
H = options.imgSizeH;
W = options.imgSizeW;
file_address = options.source_dir;
img_fileName = options.img_stack_name;
backg_fileName = options.backg_stack_name;
tformX2Y = options.tformR2L;
toPhoton = options.toPhoton;

[transformed_centerROI(2), transformed_centerROI(1)] = transformPointsInverse(tformX2Y, imgCenterY(2)*pixSize, imgCenterY(1)*pixSize);
imgCenterX = round(transformed_centerROI/pixSize);
imgCenterX = imgCenterX+[registershifty,registershiftx];


imgR = Tiff([file_address,img_fileName,'.tif'],'r');
count = 1:length(indx_stack);
SMLM_img = nan(H,W,length(count));
indx_start = indx_stack(1)-1;
    
for i=count
    
    setDirectory(imgR,i+indx_start);
    SMLM_img(:,:,i) = double(imgR.read);

end

region = -(imgSize-1)/2:(imgSize-1)/2;
SMLM_imgx = SMLM_img(:,W/2+1:end,:);
SMLM_imgx = SMLM_imgx(imgCenterX(1)+region,imgCenterX(2)+region,:);
SMLM_imgy = flip(SMLM_img(:,1:W/2,:),2);
SMLM_imgy = SMLM_imgy(imgCenterY(1)+region,imgCenterY(2)+region,:);
SMLM_img = [SMLM_imgx,SMLM_imgy];
SMLM_img = SMLM_img*toPhoton;


count = 1:length(indx_stack);
BKG_img = nan(H,W,length(count));
indx_start = indx_stack(1)-1;
    
for i=count
    imgR = Tiff([file_address,backg_fileName,num2str(i+indx_start)],'r');
    setDirectory(imgR,1);
    BKG_img(:,:,i) = double(imgR.read);

end

BKG_imgx = BKG_img(:,W/2+1:end,:);
BKG_imgx = BKG_imgx(imgCenterX(1)+region,imgCenterX(2)+region,:);
BKG_imgy = flip(BKG_img(:,1:W/2,:),2);
BKG_imgy = BKG_imgy(imgCenterY(1)+region,imgCenterY(2)+region,:);
BKG_img = [BKG_imgx,BKG_imgy];
BKG_img = BKG_img*toPhoton;






end