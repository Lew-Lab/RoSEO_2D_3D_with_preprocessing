%%
addpath(genpath('C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\RoSEO3D_v2_for_OSF\'));
f_forwardModel = @(x,G_re,b) ((G_re*x))+b;
f_loss = @(Iobs,Iest) sum(Iest-Iobs.*log(Iest+10^-16));

imgSize = 91;
imgPara.img_sizex=imgSize;
imgPara.img_sizey=imgSize;


shift = [0,0];
range = [1:91];


%frame = 1792;
frame = 0;
frame2 = 0;
for ii=1:200
   
   SM_est_cur = [];
   while isempty(SM_est_cur)
       frame2 = frame2+1;
       frame = frame+1;
        SM_est_cur = SM_est_save_all(SM_est_save_all(:,1)==frame,:);
        
   end
   [gamma,loc] = SM_est2gamma(SM_est_cur,imgPara);
   N = size(loc,2);
   gamma = reshape(gamma,[],1);
   b = reshape(SM_bkg(:,:,frame2),[],1)*toPhoton;
    [G_re,loc_re_new,~] = update_basisMatrix(N,gamma,loc,imgPara,imgPara_B);
    I = f_forwardModel(gamma,G_re,b);
    I_rec=reshape(I,imgSize,imgSize*2);
    I_rec = [I_rec(range+shift(1),range+shift(2)),I_rec(range+shift(1),range+shift(2)+imgSize)];
    %plotI(I_rec,18);
    %saveas(gca,['C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210414\reconst_image\' num2str(ii),',jpg'],'jpg');
    I_raw=reshape(SM_img(:,:,frame2)*toPhoton,imgSize,imgSize*2);
    I_raw = [I_raw(range,range),I_raw(range,range+imgSize)];
    Fig=plotI(I_raw,I_rec,30,30);
    %plotI(zeros(101,202),zeros(101,202),20,20);
    saveas(Fig,['C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210504\15-24\reconstructed_v2\' num2str(ii),'.png'],'png');
    close all;
    

end


function Fig = plotI(I_raw,I_rec,Imax1,Imax2)
temp = hot(300); map1(:,1) = min(temp(:,3),1); 
map1(:,2) = min(temp(:,1)*0.50+temp(:,2)*0.50,1);
map1(:,3) = min(temp(:,1)*1,1);

imgSz = size(I_raw,1);
Fig = figure('Color',[0.1,0.1,0.1]);


set(Fig, 'Units','centimeters','InnerPosition', [12 12 4 8]);
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
line([imgSz-10.2564*1.667-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2);
%set(gcf, 'InvertHardcopy', 'off');


I = I_rec;
%set(Fig, 'Units','centimeters','InnerPosition', [8 8 4 4]);
%set(Fig, 'Units','centimeters','InnerPosition', [12 12 4 4]);
Ix = I(:,1:imgSz); Iy = I(:,imgSz+1:imgSz*2);
%Ix= imgaussfilt(Ix,1); Iy = imgaussfilt(Iy,1);
ax1 = axes('Position',[0.1 0.2 .4 .4]);
%ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Ix]);  axis image; caxis([0,Imax2]);  axis off; colormap(ax1,'hot');  colorbar('Position',[0.13,0.25,0.3,0.02],'Location','south','Color','w');
%ax1 = axes('Position',[0.1 0.2 .4 .4]);
ax2 = axes('Position',[0.5 0.2 .4 .4]); 
imagesc([Iy]); axis image
axis off; axis image;  caxis([0,Imax2]);  colormap(ax2,map1);
line([0 0],[0 imgSz],'Color',[196, 192, 192]./255,'LineWidth',2); 
line([imgSz-10.2564*1.667-4,imgSz-4],[imgSz-5 imgSz-5],'Color',[254, 254, 254]./255,'LineWidth',2); colorbar('Position',[0.55,0.25,0.3,0.02],'Location','south','Color','w');
set(gcf, 'InvertHardcopy', 'off');

annotation('textbox', [0.3, 0.91, 0, 0], 'string', 'Experimental','Color','w','FontSize',8);
annotation('textbox', [0.3, 0.59, 0, 0], 'string', 'Reconstructed','Color','w','FontSize',8);
annotation('textbox', [0.17, 0.87, 0, 0], 'string', 'x-pol','Color','w','FontSize',8);
annotation('textbox', [0.56, 0.87, 0, 0], 'string', 'y-pol','Color','w','FontSize',8);

end
