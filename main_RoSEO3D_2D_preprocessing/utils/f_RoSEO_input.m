function f_RoSEO_input(RoSEOPara,dataFileFolder,inter_results_save_folder,ID,Nimg,imgSize,FoV_analysed,crop_ROI_center,resultsSaveFolder,thred)
%% readme:
% add the data folder, and the algorithm folder into the path
% change the fileName in Line 39,40

%% parameters of the microscope

%----------------------------------------------
% create nanoscope object analysis
emitter_wavelength = RoSEOPara.emitter_wavelength; %nm
refractiveIndx=RoSEOPara.refractiveIndx;% objective refractive index
sampleRefractiveIndx=RoSEOPara.sampleRefractiveIndx;% sample refractive index

left_to_right_trans_ratio = RoSEOPara.left_to_right_trans_ratio; %1.1426 for 488/561+523/610
zerothOrder = [0,0];%zerothOrder_RL; 
phasemaskpara.zeroorder = zerothOrder; %usually used for PSFs with large footprint, and some light leak into the center of PSF
maskName = fullfile('phasemask', 'standard');
phasemaskpara.maskname = char(RoSEOPara.maskName);
if RoSEOPara.calibrated_mask~=""
phasemaskpara.calibrated_mask = char(RoSEOPara.calibrated_mask);   % if use calibrated pmask, give the retreved mask filed name
pmask = phasemaskpara.calibrated_mask;
else 
    pmask = [];
end
toPhoton = RoSEOPara.toPhoton;  %gain of the camera electron/ADU
imgPara.pix_sizex = RoSEOPara.pix_sizex; %the griding unit for RoSEO, in unit of nm
reg = RoSEOPara.reg;

%%
% build microscope object used in the RoSEO

n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

[FPSFx, FPSFy] = n1.createPSFstruct(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio);  
%
basis_matrix = [reshape([n1.XXxBasis,n1.XXyBasis],[],1),...
               reshape([n1.YYxBasis,n1.YYyBasis],[],1),...
               reshape([n1.ZZxBasis,n1.ZZyBasis],[],1),...
               reshape([n1.XYxBasis,n1.XYyBasis],[],1),...
               reshape([n1.XZxBasis,n1.XZyBasis],[],1),...
               reshape([n1.YZxBasis,n1.YZyBasis],[],1)];

h = visualizeBases(n1);
exportgraphics(h,strcat(inter_results_save_folder,'\','BasisImage_analyzed_by_RoSEO.jpg'),'Resolution',600);
close(h);
%% parameter of the data analyse process

for mm = 1:length(ID)

    for hh = 1:length(FoV_analysed)
       
%% simulate images
FoV_cur = FoV_analysed(hh);
note = '';
suffix_name = ['_centerY_y', num2str(crop_ROI_center(1)), '_x_', num2str(crop_ROI_center(2)), '_FoV', num2str(imgSize),'_', num2str(imgSize),'_', num2str(FoV_cur),'th_FoV'];

offSetname = [char(dataFileFolder),'data_offset',suffix_name,'.mat']; load(offSetname);
SMLMName = [char(dataFileFolder),'data',num2str(ID(mm)),suffix_name,'.tif'];
bkg = load([char(dataFileFolder),'data',num2str(ID(mm)),'_bkg',suffix_name,'.mat']);

Nimg_perStack = 32*8;   %# of images been analysed in each parallel loop
stacks = ceil(Nimg/Nimg_perStack);

SM_est_save = {};
NLL_save = {};
SMLMR = Tiff(SMLMName,'r');
SMLM_bkg_temp = bkg.SMLM_bkg;  
for ii = 1:stacks
    count = 0;
    for jj = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg)
        count = count+1;
        SMLMR = Tiff(SMLMName,'r');
        setDirectory(SMLMR,jj);
        SM_img(:,:,count) = double(SMLMR.read)-offset;
        SM_bkg(:,:,count) = SMLM_bkg_temp(:,:,ceil(jj/50))+1;
    end
 
parfor kk = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg)
   
    count = rem(kk,Nimg_perStack);
    if count ==0
        count = Nimg_perStack;
    end
   
SMLM_img = ([SM_img(:,:,count)])*toPhoton;
SMLM_bkg1= [SM_bkg(:,:,count)]*toPhoton;

%SMLM_bkg1 = ones(size(SMLM_bkg1))*9;
%SMLM_bkg1 = [ones(imgSize,imgSize)*mean(SMLM_bkg1(:,ceil(1:imgSize/4)),'all'),...
%             ones(imgSize,imgSize)*mean(SMLM_bkg1(:,ceil(1:imgSize/4)+imgSize),'all')];
SMLM_img(SMLM_img<0)=0;

[~,~,est] = RoSEO(n1, SMLM_img, SMLM_bkg1, FPSFx, FPSFy,thred,'regval', reg);

SM_est_save{kk} = est;


end
save([char(resultsSaveFolder), 'est_temp.mat'],'SM_est_save');


frame_oder = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg);
for framekk = [1:5]
Fig1 = figure('Visible','off','Units','inches','InnerPosition',[1,1,10,8]);
imagesc([SM_img(:,:,framekk)]*toPhoton); axis image; hold on;  colorbar; caxis([0,100]); axis off;
est = SM_est_save{frame_oder(framekk)};
if length(est)>0
scatter(est(:,2)/58.5+imgSize/2-1,est(:,3)/58.5+imgSize/2-1,'r+','LineWidth',1.5);scatter(est(:,2)/58.5+imgSize/2-1+imgSize,est(:,3)/58.5+imgSize/2-1,'r+','LineWidth',1.5);
exportgraphics(Fig1,strcat(inter_results_save_folder,'\','data',num2str(ID(mm)),'FoV',num2str(FoV_cur),'frame',num2str(frame_oder(framekk)),'RoSEO_est_results_example.jpg'),'Resolution',400);
end
end
close all

end





%% project the 6 moments estimation to 3 angulular space

mean_bkg = 9;%mean(mean(mean(mean(SM_bkg))))*toPhoton;

Angle_save_cell={};
SM_est_save_all_cell = {};
Angle_save = [];
SM_est_save_all = [];

parfor kk = 1:Nimg
    SM_est = SM_est_save{kk};

    Angle_frame = [];
    SM_est_frame = [];

if isempty(SM_est)==0
    SM_est(SM_est(:,4)<0,:) = [];
    
    %
    saveAngle = [];
    saveOren = [];

    Angle_frame = [];
    SM_est_frame = [];
    for ll = 1:size(SM_est,1)
        [mux,muy,muz,rotMobil] = secondM2SymmCone_for_experiment(double(SM_est(ll,5:10)),double(SM_est(ll,4)),mean_bkg,basis_matrix);
        if muz<=0
            mux = -mux;
            muy = -muy;
            muz = -muz;
        end
        saveOren(ll,:) = [mux,muy,muz,rotMobil];
        [thetaD, phiD, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
        saveAngle = [thetaD,phiD,alphaD,rotMobil,3*pi-sqrt(rotMobil*8*pi^2+pi^2)];

        Angle_frame(ll,:) = [Angle_save;kk,saveAngle];
        SM_est_frame(ll,:) = [SM_est_save_all;double(SM_est(ll,:))];
    end
end

Angle_save_cell{kk} = Angle_frame;
SM_est_frame_save_cell{kk} = SM_est_frame;

end

%%

for kk = 1:length(Angle_save_cell)
    Angle_save = [Angle_save;Angle_save_cell{kk}];
    SM_est_save_all = [SM_est_save_all;SM_est_frame_save_cell{kk}];

end

Fig1 = figure('Visible','off','Units','inches','InnerPosition',[1,1,10,8]);
scatter(SM_est_save_all(:,2),SM_est_save_all(:,3),4,'filled'); axis image;
exportgraphics(Fig1,strcat(inter_results_save_folder,'\','data',num2str(ID(mm)),'FoV',num2str(FoV_cur),'RoSEO_est_results_all.jpg'),'Resolution',400);

close all;

save([char(resultsSaveFolder), num2str(ID(mm)) '_est_FoV',num2str(FoV_cur),'.mat'],'SM_est_save','Angle_save','SM_est_save_all','left_to_right_trans_ratio','pmask','reg','note');
end
%%

end
end

