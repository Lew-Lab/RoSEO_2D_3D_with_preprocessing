% dataN = [1];
% fig_saveFolder = 'F:\OneDrive - Washington University in St. Louis\github\Deep-SMOLM3D\3D SLB vidualization\SLB data\DPPC_chol_Figs\';
% readFunction = 'data_overview.DPPC_chol_dataSource';
% generate_diffuse_check_frame = 0;
% readFunction = str2func(readFunction);


% dataN = [1];
% fig_saveFolder = 'F:\OneDrive - Washington University in St. Louis\github\Deep-SMOLM3D\3D SLB vidualization\SLB data\STORM_cell_Figs\';
% readFunction = 'data_overview.STROM_cell';
% generate_diffuse_check_frame = 0;
% readFunction = str2func(readFunction);


dataN = [4];
fig_saveFolder = '###\3D SLB vidualization_demo\SLB data\DPPC_chol_Figs\';
readFunction = 'data_overview.DPPC_chol_dataSource';
 readFunction = str2func(readFunction);


if ~exist([fig_saveFolder,'reconstrunction_check\'], 'dir')
   mkdir([fig_saveFolder,'reconstrunction_check\'])
end


%%
addpath(genpath('###\RoSEO3D_in_use\'));
generateMicroscope = 1; %code for generate the forward model related parameter; set to be 0, if you already has this file in the above folder
%

if generateMicroscope==1
RoSEOPara = struct('NFP',-500, ...
                    'imgSize',91,...
                    'emitter_wavelength',610, ...
                    'sampleRefractiveIndx',1.314, ...
                    'refractiveIndx',1.515, ...
                    'left_to_right_trans_ratio',1.1426, ...
                    'maskName','pixOL_v12_com', ...
                    'calibrated_mask','pmask_retrieved_data12_data0602.mat', ...
                    'toPhoton',0.29, ...
                    'pix_sizex',58.5, ...
                    'pix_sizez',50, ...
                    'axial_grid_points',[-200:50:1300], ...
                    'twoD_fitRange', [0,900]);
[imgPara,imgPara_B] = f_generate_microscope(RoSEOPara,RoSEOPara.imgSize);
save([fig_saveFolder,'reconstrunction_check\data',num2str(dataN),'_Microscope_data.mat'],'imgPara','RoSEOPara','imgPara_B');
else
load([fig_saveFolder,'reconstrunction_check\data',num2str(dataN),'_Microscope_data.mat']);
end
%% start the reconstruction

f_forwardModel = @(x,G_re,b) ((G_re*x))+b;
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));

imgSize = RoSEOPara.imgSize;


%% read the estimation/reconstruction related data

[locFileName,imageName,angleEstName,center_SLB,radius_SLB,dataIDN,frameN_per_stack,bkgName, offsetName] = readFunction(dataN);
load(locFileName);
SM_est_save_all = locList;
Angle_est_all = angleList;

% read background and raw image
load(bkgName);
load(offsetName);
for ii = 1:frameN_per_stack
    readT = Tiff(imageName,'r');
    setDirectory(readT,ii);
    SMLM_img(:,:,ii) = read(readT);
end
SMLM_img= double(SMLM_img);
%%
frameN_demo = 500;
frame = 0;

clear temp temp1
temp = []; temp1 = [];

for ii=1:frameN_demo
   
   SM_est_cur = [];
%    while isempty(SM_est_cur)
%        frame = frame+1;
%         SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);
%         
%    end


   while isempty(SM_est_cur)
       frame = frame+1;
        SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);
        angle_cur = Angle_est_all(SM_est_save_all(:,1)==frame,:);
   end

   %if SM_est_cur(:,4)>=700 &  SM_est_cur(:,4)<=800
%temp = [temp;SM_est_cur];
%temp1 = [temp1;angle_cur];

       frame = frame+1;
       SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);

   if isempty(SM_est_cur)
       SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==1,:)*0;
   else

        
   end
   [gamma,loc] = SM_est2gamma(SM_est_cur,imgPara);
   N = size(loc,2);
   gamma = reshape(gamma,[],1);
   b = reshape(SMLM_bkg(:,:,ceil(frame/50)),[],1)*RoSEOPara.toPhoton;
   [G_re,loc_re_new,~] = update_basisMatrix(N,gamma,loc,imgPara,imgPara_B);
   %I = f_forwardModel(gamma,G_re,b*2.2759)/2.2759;
   I = f_forwardModel(gamma,G_re,b);
   I_rec=reshape(I,imgSize,imgSize*2);
   I_raw=(SMLM_img(:,:,frame)-offset)*RoSEOPara.toPhoton;
   
   Fig=plotI_scaled(I_raw,I_rec,30,30,3);

   exportgraphics(Fig,[fig_saveFolder,'reconstrunction_check\data' num2str(dataN),'_frame',num2str(frame),'.jpg'],"Resolution",150,'BackgroundColor','black');
   close(Fig);
   %end

end
% close(v);










%%

function Fig = plotI_scaled(I_raw,I_rec,Imax1,Imax2,scaled_size)
I_raw = imresize(I_raw,scaled_size,'nearest');
I_rec = imresize(I_rec,scaled_size,'nearest');

temp = hot(300); map1(:,1) = min(temp(:,3),1); 
map1(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1);
map1(:,3) = min(temp(:,1)*1,1);

imgSz = size(I_raw,1);
Fig = figure('Color',[0.1,0.1,0.1]);


set(Fig, 'Units','centimeters','InnerPosition', [12 2 4*3 8*3]);
I = I_raw;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);

Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.5 .4 .4]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax1]);  axis off; colormap(ax1,'hot');  
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.5 0.5 .4 .4]); 
imagesc([Iy]); axis image
axis off; axis image;  caxis([0,Imax1]);  colormap(ax2,map1); 
line([0 0],[0 imgSz],'Color',[196, 192, 192]./255,'LineWidth',2); 
%line([imgSz-10.2564*1.667-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2);
%set(gcf, 'InvertHardcopy', 'off');


I = I_rec;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);
%set(Fig, 'Units','centimeters','InnerPosition', [12 12 4 4]);
Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.2 .4 .4]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax2]);  axis off; colormap(ax1,'hot');  colorbar('Position',[0.13,0.25,0.3,0.02],'Location','south','Color','w','FontSize',11);
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Iy]); axis image
axis off; axis image;  caxis([0,Imax2]);  colormap(ax2,map1);
line([0 0],[0 imgSz],'Color',[196, 192, 192]./255,'LineWidth',4); 
line([imgSz-(10.2564*1.667)*3-4*3,imgSz-4*3],[imgSz-5*3 imgSz-5*3],'Color',[254, 254, 254]./255,'LineWidth',2); colorbar('Position',[0.55,0.25,0.3,0.02],'Location','south','Color','w','FontSize',11);
set(gcf, 'InvertHardcopy', 'off');

annotation('textbox', [0.45, 0.88, 0, 0], 'string', 'Raw','Color','w','FontSize',18);
annotation('textbox', [0.35, 0.56, 0, 0], 'string', 'Reconstructed','Color','w','FontSize',18);
annotation('textbox', [0.25, 0.84, 0, 0], 'string', 'x-pol','Color','w','FontSize',15);
annotation('textbox', [0.62, 0.84, 0, 0], 'string', 'y-pol','Color','w','FontSize',15);

end



