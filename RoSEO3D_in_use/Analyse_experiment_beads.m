% change fileName in line 27
%% parameters of the microscope

%----------------------------------------------
% create nanoscope object analysis

emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1.1426;
zerothOrder = [0,0];%zerothOrder_RL; 
maskName = fullfile('phasemask', 'pixOL_v12');
phasemaskpara.zeroorder = zerothOrder; %usually used for PSFs with large footprint, and some light leak into the center of PSF
phasemaskpara.maskname = maskName;
phasemaskpara.calibrated_mask = 'pmask_retrieved_data12_data0602.mat'; 
pmask = phasemaskpara.calibrated_mask;
NFP=50;  % normal focal plane of the microscope when imaging
toPhoton = 0.28;%0.28; %gain of the camera electron/ADU



%% parameter of the data analyse process

imgPara.pix_sizex = 58.5; %the griding unit for RoSEO3D, in unit of nm
imgPara.pix_sizez = 50; % griding size along z direction
imgPara.axial_grid_points = [-16:1:16]*imgPara.pix_sizez;  % build the z slices 

scaling_factor = [100,100,100];   % used to balence the weights on orientation estimation and location estimation, experimental trier number
imgPara.scaling_factor = scaling_factor;

imgSize = 71;  % image size = max(imagesizeX, imagesizeY)
Nimg = 319; % total number of frames in your stack

reg = 0.22; % regulization for step 1 in RoSEO
note = ''; % for adding additional notes
fileName = [pwd '\example_data\beads_data\'];
%% build microscope object used in the RoSEO3D step2 and step3

n1 = Nanoscope_beads('imageSize', imgSize*2+21,... 
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
imgPara.PSF_size_opt_min = 11; % the minimum size for croping the PSF
%% build microscope object used in the RoSEO3D step2 and step3

n1 = Nanoscope_beads('imageSize', imgSize,...
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

[FPSFx, FPSFy] = n1.createPSFstruct(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',NFP,...
    'molecule_plane',[-500,500]);   % [z min, z max]; the values are the expected z range for SMs, no need to be accurate

%% initialize a parallel pool

p = gcp('nocreate');
if isempty(p)
    parpool('local', 32)
end


%% simulate images
ID = [6,12];
for mm = 1:length(ID)

load(['data',num2str(ID(mm)),'_beads.mat']);
SM_img = pixOL_beads;
imgSz = size(SM_img,1);
Nimg = 319;

parfor kk = 1:Nimg
    %kk
    %SMLM_img = poissrnd(I);
    %SMLM_bkg = bkg;
SMLM_img = SM_img(:,:,kk)*toPhoton;
pixel_numb = imgSz^2-(imgSz-4)^2;
SMLM_bkg = [ones(imgSz,imgSz)*(sum(SMLM_img(1:imgSz,1:imgSz),'all')-sum(SMLM_img(3:imgSz-2,3:imgSz-2),'all'))/pixel_numb,ones(imgSz,imgSz)*(sum(SMLM_img(1:imgSz,(1:imgSz)+imgSz),'all')-sum(SMLM_img(3:imgSz-2,(3:imgSz-2)+imgSz),'all'))/pixel_numb];

%range = [1:67]+2;
%shift = [0,0];
%SMLM_img = [SM_img(range,range,kk),SM_img(range+shift(1),range+71+shift(2),kk)]*toPhoton;
%SMLM_bkg = [ones(67,67)*(sum(SM_img(1:71,1:71,kk),'all')*toPhoton-sum(SMLM_img(1:67,1:67),'all'))/552,ones(67,67)*(sum(SM_img(1:71,(1:71)+71,kk),'all')*toPhoton-sum(SMLM_img(1:67,(1:67)+67),'all'))/552];

[SM_est,NLL_cur_out] = RoSEO3D_beads(SMLM_img, SMLM_bkg,n1,FPSFx, FPSFy,imgPara,imgPara_B,'regval',reg);

SM_est_save{kk} = SM_est;
NLL_save{kk}=NLL_cur_out;
end
SMLM_img = SM_img(:,:,1)*toPhoton;
pixel_numb = imgSz^2-(imgSz-4)^2;
SMLM_bkg = [ones(imgSz,imgSz)*(sum(SMLM_img(1:imgSz,1:imgSz),'all')-sum(SMLM_img(3:imgSz-2,3:imgSz-2),'all'))/pixel_numb,ones(imgSz,imgSz)*(sum(SMLM_img(1:imgSz,(1:imgSz)+imgSz),'all')-sum(SMLM_img(3:imgSz-2,(3:imgSz-2)+imgSz),'all'))/pixel_numb];

mean_bkg = mean(mean(mean(SMLM_bkg)));

%% project the 6 moments estimation to 3 angulular space
% ***
Angle_save=[];
SM_est_save_all = [];
NLL_save_all = [];
for kk = 1:Nimg
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
save([num2str(ID(mm)) '_est_v4_210602.mat'],'SM_est_save','NLL_save','Angle_save','SM_est_save_all','NLL_save_all','left_to_right_trans_ratio','NFP','pmask','reg','note');
end
%%
