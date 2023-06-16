%% readme:
% add the data folder, and the algorithm folder into the path
% change the fileName in Line 39,40

%% parameters of the microscope

%----------------------------------------------
% create nanoscope object analysis
emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1.1426;
zerothOrder = [0,0];%zerothOrder_RL; 
phasemaskpara.zeroorder = zerothOrder; %usually used for PSFs with large footprint, and some light leak into the center of PSF
maskName = fullfile('phasemask', 'pixOL_v12_com');
phasemaskpara.maskname = maskName;
phasemaskpara.calibrated_mask = 'pixOL_com_retrieved_data12_data0602.mat';   % if use calibrated pmask, give the retreved mask filed name
pmask = phasemaskpara.calibrated_mask;
NFP=-500;  % normal focal plane of the microscope when imaging
toPhoton = 0.66; %gain of the camera electron/ADU



%% parameter of the data analyse process

imgPara.pix_sizex = 58.5; %the griding unit for RoSEO3D, in unit of nm
imgPara.pix_sizez = 50; % griding size along z direction
imgPara.axial_grid_points = [-6:1:28]*imgPara.pix_sizez;  % build the z slices 

scaling_factor = [100,100,100];   % used to balence the weights on orientation estimation and location estimation, experimental trier number
imgPara.scaling_factor = scaling_factor;

imgSize = 91;  % image size = max(imagesizeX, imagesizeY)
Nimg = 2000; % total number of frames in your stack

%--------------------------------------------------
% data folder name information
fileName = [pwd '\example_data\20210504 3D DPPC_with_chol_withOffset\'];
ID = [15]; %prefix of the data
suffix_name = '_beads_w_offset.tif'; % suffix of name of the SM image 

reg = 0.22; % regulization for step 1 in RoSEO
note = ''; % for adding additional notes

offset = imread(['offSet.tif']); % offset name 
%% build microscope object used in the RoSEO3D step2 and step3

n1 = Nanoscope('imageSize', imgSize*2+41,... 
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

%
[PSFx, PSFy] = n1.createPSFstruct3D(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',NFP,...
    'axial_grid_points',imgPara.axial_grid_points,...
    'scaling_factor',scaling_factor); %distance is is nm

% define image parameters

Bx = cat(4,PSFx.XXx,PSFx.YYx,PSFx.ZZx,PSFx.XYx,PSFx.XZx,PSFx.YZx,...
     PSFx.XXxdx,PSFx.YYxdx,PSFx.ZZxdx,...
     PSFx.XXxdy,PSFx.YYxdy,PSFx.ZZxdy,...
     PSFx.XXxdz,PSFx.YYxdz,PSFx.ZZxdz);
By = cat(4,PSFy.XXy,PSFy.YYy,PSFy.ZZy,PSFy.XYy,PSFy.XZy,PSFy.YZy,...
     PSFy.XXydx,PSFy.YYydx,PSFy.ZZydx,...
     PSFy.XXydy,PSFy.YYydy,PSFy.ZZydy,...
     PSFy.XXydz,PSFy.YYydz,PSFy.ZZydz);

imgPara.img_size = imgSize;
imgPara.PSF_size_opt = 21; % the best size for croping the PSF (pixels)
imgPara.PSF_size = imgPara.PSF_size_opt;
imgPara_B.Bx = Bx;
imgPara_B.By = By;
imgPara.PSF_size_opt_min = 11; % the minimum size for croping the PSF ,use 21 for my experimental pixOL data
combine_thred = 21;
%% build microscope object used in the RoSEO3D step2 and step3

n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

[FPSFx, FPSFy] = n1.createPSFstruct(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',NFP,...
    'molecule_plane',[0,900]);   % [z min, z max]; the values are the expected z range for SMs, no need to be accurate

%% initialize a parallel pool

p = gcp('nocreate');
if isempty(p)
    parpool('local', 32)
end


%% simulate images

for mm = 1:length(ID)

SMLMName = [fileName,num2str(ID(mm)),suffix_name];

Nimg_perStack = 32*8;   %# of images been analysed in each parallel loop
stacks = ceil(Nimg/Nimg_perStack);

SM_est_save = {};
NLL_save = {};
SMLMR = Tiff(SMLMName,'r');
for ii = 1:stacks
    count = 0;
    for jj = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg)
        count = count+1;
        SMLMR = Tiff(SMLMName,'r');
        setDirectory(SMLMR,jj);
        SM_img(:,:,count) = double(SMLMR.read)-offset;
        
        SM_bkg(:,:,count) = imread([fileName,'bkg/data',num2str(ID(mm)),'frame',num2str(jj),'.tif']); %loading background
    end
    
for kk = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg)
    count = rem(kk,Nimg_perStack);
    if count ==0
        count = Nimg_perStack;
    end
   
SMLM_img = ([SM_img(:,:,count)])*toPhoton;
SMLM_bkg1= [SM_bkg(:,:,count)]*toPhoton;
SMLM_img(SMLM_img<0)=0;
[SM_est,NLL_cur_out,~] = RoSEO3D(SMLM_img, SMLM_bkg1,n1,FPSFx, FPSFy,imgPara,imgPara_B, 'regval',reg);

SM_est_save{kk} = SM_est;
NLL_save{kk}=NLL_cur_out;
end
end

mean_bkg = mean(mean(mean(mean(SM_bkg))))*toPhoton;

%% project the 6 moments estimation to 3 angulular space
% ***
Angle_save=[];
SM_est_save_all = [];
NLL_save_all = [];
for kk = 2%1:Nimg
    SM_est = SM_est_save{kk};
if isempty(SM_est)==0
    NNL_cur = NLL_save{kk}.';
    NLL_save_all = [NLL_save_all;kk*ones(size(NNL_cur)),NNL_cur];
    SM_est(SM_est(:,4)<0,:) = [];
    
    %
    saveAngle = [];
    saveOren = [];
    for ll = 1:size(SM_est,1)
        [mux,muy,muz,rotMobil] = secondM2SymmCone_RoSEO3D(double(SM_est(ll,:)),mean_bkg,imgPara,imgPara_B);
        if muz<=0
            mux = -mux;
            muy = -muy;
            muz = -muz;
        end
        saveOren(ll,:) = [mux,muy,muz,rotMobil];
        [thetaD, phiD, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
        saveAngle = [thetaD,phiD,alphaD,rotMobil,3*pi-sqrt(rotMobil*8*pi^2+pi^2)];

        Angle_save = [Angle_save;kk,saveAngle];
        SM_est_save_all = [SM_est_save_all;kk,double(SM_est(ll,:))];
    end
end


end

%%
save([num2str(ID(mm)) '_est.mat'],'SM_est_save','NLL_save','Angle_save','SM_est_save_all','NLL_save_all','left_to_right_trans_ratio','NFP','pmask','reg','note','combine_thred');
end
%%
