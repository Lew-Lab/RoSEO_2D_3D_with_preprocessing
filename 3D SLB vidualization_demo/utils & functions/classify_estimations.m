%% 2021/05/17
% based on estimation to classify the emitters at similar location and sum
% their PSF images together

%%
shift_off_set = [-43,113,-43];%x,y,z
Z = 950; phi =-135; R = 1000; r = sqrt(R^2-(R-Z)^2); X = r*cosd(phi); Y = r*sind(phi);
centers = [X,Y,Z]; %in unit of nm
%centers = [0,0,Z]; %in unit of nm
box = [50,50,100]; %the assemble box centered at the center points, in unit of nm
Name_estimation = '_est_retrieved_1.1_v26.mat';
Name_raw_data = '_3Dlipid';
cropped_image_size = 31;
raw_image_size = 101;

shiftxy = [2,0]; %[y,x] in unit of pixel


%% find the points located close to the center
count = 0;
N_SMs = 0;
load('R_shift_each_stack.mat');
cropped_img = zeros(cropped_image_size,cropped_image_size*2);
for ii = [15:24]

    count = count+1;
    shift_off_set = -[R_save(count,1:2),-50];
    center_pixel = round([centers(1)-shift_off_set(1),+centers(2)-shift_off_set(2)]/58.5+(raw_image_size+1)/2); %in unit of pixels, the z axial position won't be used. 
    center_nm = centers-shift_off_set;
    radius = (cropped_image_size-1)/2;
    range = [center_pixel(2)+[-radius:radius];center_pixel(1)+[-radius:radius]];
    


    indx_data = num2str(ii);
    load([indx_data,Name_estimation]);
    load(['raw_data\',indx_data,Name_raw_data,'_bkg_v4.mat'])
    distance_all = abs([SM_est_save_all(:,2)-center_nm(1),SM_est_save_all(:,3)-center_nm(2),SM_est_save_all(:,4)-center_nm(3)]);
    candicate_indx = distance_all(:,1)<box(1)&distance_all(:,2)<box(2)&distance_all(:,3)<box(3);
    N_SMs = N_SMs+sum(candicate_indx);
%crop the image

fileName = ['raw_data\',indx_data,Name_raw_data,'.tif'];


imagR = Tiff(fileName,'r');
for kk = 1:size(SM_est_save_all,1)
if candicate_indx(kk)
    
    setDirectory(imagR,SM_est_save_all(kk,1));
    SM_img = double(imagR.read);
    backg_l = SMLM_bkg(range(1,:),range(2,:),SM_est_save_all(kk,1));
    backg_r = SMLM_bkg(range(1,:)+shiftxy(1),range(2,:)+shiftxy(2)+raw_image_size,SM_est_save_all(kk,1));
    
    cropped_img = cropped_img+[(SM_img(range(1,:),range(2,:))-backg_l)*1.1426,(SM_img(range(1,:)+shiftxy(1),range(2,:)+shiftxy(2)+raw_image_size)-backg_r)]*1/3.4683;
    
    %axisOFF = 'on';
    %plotI([(SM_img(range(1,:),range(2,:))-backg_l)*1.1426,(SM_img(range(1,:)+shiftxy(1),range(2,:)+shiftxy(2)+raw_image_size)-backg_r)]*1/3.4683,axisOFF);
end
end
end

%%
N_SMs
axisOFF = 'off';
isScaleBar = 'off';
%plotI([(SM_img(range(1,:),range(2,:))-backg_l)*1.1426,(SM_img(range(1,:)+shiftxy(1),range(2,:)+shiftxy(2)+raw_image_size)-backg_r)]*1/3.4683,axisOFF);
plotI(cropped_img,axisOFF,isScaleBar);
%%

function plotI(I_raw,axisoff,isScaleBar)
Imax1 = max(I_raw,[],'all');
temp = hot(300); map1(:,1) = min(temp(:,3),1); 
map1(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1);
map1(:,3) = min(temp(:,1)*1,1);

imgSz = size(I_raw,1);
%Fig = figure('Color',[0.1,0.1,0.1]);
Fig = figure();
set(Fig, 'Units','centimeters','InnerPosition', [12 12 4*0.8 8*0.8]);
I = I_raw;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);

Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.5 .40 .20]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax1]);   colormap(ax1,'hot');  %colorbar('Location','north'); 
if strcmp(axisoff, 'off')    
    axis off;
end
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.51 0.5 .40 .20]); 
imagesc([Iy]); axis image
 axis image;  caxis([0,Imax1]);  colormap(ax2,map1); %colorbar('Location','north');  axis off;
 if strcmp(axisoff, 'off')
    axis off;
end
%line([1 1],[1 imgSz],'Color',[196, 192, 192]./255,'LineWidth',2); 
if strcmp(isScaleBar, 'on')
line([imgSz-8.5470-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2);
end
%set(gcf, 'InvertHardcopy', 'off');

end


