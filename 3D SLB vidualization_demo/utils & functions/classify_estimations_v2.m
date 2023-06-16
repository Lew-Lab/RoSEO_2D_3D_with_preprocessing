%% 2021/05/17
% based on estimation to classify the emitters at similar location and sum
% their PSF images together

%%
shift_off_set = [-138,+17,-173];%x,y,z
Z = 900; phi = -45; R = 1000;
centers = [0,0,Z]; %in unit of nm
%centers = [0,0,Z]; %in unit of nm
box = [5000,5000,50]; %the assemble box centered at the center points, in unit of nm
Name_estimation = '_est_retrieved_1.1_v5.mat';
Name_raw_data = '_beads';
Name_raw_data_bkg = '_3Dlipid_bkg_v2';
cropped_image_size = 41;
raw_image_size = 101;


f_forwardModel = @(x,G_re,b) (G_re*x)+b;
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));
imgPara.img_sizex=101;
imgPara.img_sizey=101;
%% find the points located close to the center
cropped_img = zeros(cropped_image_size,cropped_image_size*2);
N_SMs = 0;
for ii = [2:10]
    indx_data = num2str(ii);
    load([indx_data,Name_estimation]);
    load(['F:\data\20210524 lipid on 3D beads\offset_subtracted\1-10\',indx_data,Name_raw_data_bkg,'.mat'])
    distance_all = abs([SM_est_save_all(:,2)+shift_off_set(1)-centers(1),SM_est_save_all(:,3)+shift_off_set(2)-centers(2),SM_est_save_all(:,4)+shift_off_set(3)-centers(3)]);
    candicate_indx = distance_all(:,1)<box(1)&distance_all(:,2)<box(2)&distance_all(:,3)<box(3)&abs(Angle_save(:,3)-phi)<10;
    N_SMs = N_SMs+sum(candicate_indx);
%crop the image

fileName = ['F:\data\20210524 lipid on 3D beads\offset_subtracted\1-10\',indx_data,Name_raw_data,'.tif'];


imagR = Tiff(fileName,'r');
for kk = 1:size(SM_est_save_all,1)
if candicate_indx(kk)
    
    est_location = SM_est_save_all(kk,2:4);
    center_pixel = round([est_location(1)-shift_off_set(1),est_location(2)-shift_off_set(2)]/58.5+(raw_image_size+1)/2); %in unit of pixels, the z axial position won't be used. 
    radius = (cropped_image_size-1)/2;
    range = [center_pixel(2)+[-radius:radius];center_pixel(1)+[-radius:radius]];

    
    setDirectory(imagR,SM_est_save_all(kk,1));
    SM_img = double(imagR.read);
    backg_l = SMLM_bkg(range(1,:),range(2,:),SM_est_save_all(kk,1));
    backg_r = SMLM_bkg(range(1,:),range(2,:)+raw_image_size,SM_est_save_all(kk,1));
    %backg_l(:,:)=0;
    %backg_r(:,:)=0;
    cropped_img_cur = [(SM_img(range(1,:),range(2,:)))*1.1426,(SM_img(range(1,:),range(2,:)+raw_image_size))]*1/3.4683;
    cropped_img = cropped_img+[(SM_img(range(1,:),range(2,:)))*1.1426,(SM_img(range(1,:),range(2,:)+raw_image_size))]*1/3.4683;
    %figure(); imagesc(cropped_img_cur); axis image

    
    %% reconstruct image
    frame = SM_est_save_all(kk,1);
    SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);
%    while isempty(SM_est_cur)
%        frame = frame+1;
%         SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);
%    end
    if isempty(SM_est_cur)
       I = reshape(SMLM_bkg(:,:,frame)*toPhoton,[],1);
    else
     [gamma,loc] = SM_est2gamma(SM_est_cur,imgPara);
     N = size(loc,2);
     gamma = reshape(gamma,[],1);
     b = reshape(SMLM_bkg(:,:,frame)*toPhoton,[],1);
     [G_re,loc_re_new,~] = update_basisMatrix(N,gamma,loc,imgPara);
     I = f_forwardModel(gamma,G_re,b);
    end
     I_rec=reshape(I,101,202);
     reconstructed_img_cur = [(I_rec(range(1,:),range(2,:))),(I_rec(range(1,:),range(2,:)+raw_image_size))];
    plotI(reconstructed_img_cur,cropped_img_cur,15,15);
    Angle_save(kk,2:end)
end
end
end





N_SMs
figure(); imagesc(cropped_img); axis image


function plotI(I_raw,I_rec,Imax1,Imax2)
temp = hot(300); map1(:,1) = min(temp(:,3),1); 
map1(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1);
map1(:,3) = min(temp(:,1)*1,1);

imgSz = size(I_raw,1);
Fig = figure('Color',[0.1,0.1,0.1]);

set(Fig, 'Units','centimeters','InnerPosition', [5 5 4*2 8*2]);
I = I_raw;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);

Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.5 .4 .4]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax1]);  axis off; colormap(ax1,'hot');  %colorbar('Location','north');
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.5 0.5 .4 .4]); 
imagesc([Iy]); axis image
axis off; axis image;  caxis([0,Imax1]);  colormap(ax2,map1); %colorbar('Location','north');
line([0 0],[0 imgSz],'Color',[196, 192, 192]./255,'LineWidth',2); 
line([imgSz-10.2564*1.667-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2);
%set(gcf, 'InvertHardcopy', 'off');


I = I_rec;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);
%set(Fig, 'Units','centimeters','InnerPosition', [12 12 4 4]);
Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.2 .4 .4]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax2]);  axis off; colormap(ax1,'hot'); % colorbar;
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Iy]); axis image
axis off; axis image;  caxis([0,Imax2]);  colormap(ax2,map1);
line([0 0],[0 imgSz],'Color',[196, 192, 192]./255,'LineWidth',2); 
line([imgSz-10.2564*1.667-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2);
set(gcf, 'InvertHardcopy', 'off');
end
