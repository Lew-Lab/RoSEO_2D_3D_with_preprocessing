%% readme:
% add the data folder, and the algorithm folder into the path
% change the fileName in Line 39,40

%% parameters of the microscope

%----------------------------------------------
% create nanoscope object analysis
emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1.1426; %1.1426 for 488/561+523/610
zerothOrder = [0,0];%zerothOrder_RL; 
phasemaskpara.zeroorder = zerothOrder; %usually used for PSFs with large footprint, and some light leak into the center of PSF
maskName = fullfile('phasemask', 'standard');
phasemaskpara.maskname = maskName;
%phasemaskpara.calibrated_mask = 'pixOL_com_retrieved_data12_data0602.mat';   % if use calibrated pmask, give the retreved mask filed name
%pmask = phasemaskpara.calibrated_mask;
pmask = [];
toPhoton = 0.29; %gain of the camera electron/ADU



%% parameter of the data analyse process

imgPara.pix_sizex = 58.5; %the griding unit for RoSEO, in unit of nm


imgSize = 91;  % image size = max(imagesizeX, imagesizeY)
Nimg = 10; % total number of frames in your stack

%--------------------------------------------------
% data folder name information
fileName = ['E:\Experimental_data\20220814 A1-LCD\processed data RoSEO\data43-47\'];
ID = [43:47]; %prefix of the data
suffix_name = '_centerY_y318_x_177_FoV91_91_1th_FoV'; % suffix of name of the SM image 

reg = 0.5; % regulization for step 1 in RoSEO
note = ''; % for adding additional notes

load([fileName 'offset_centerY_y318_x_177_FoV91_91_1th_FoV.mat']); % offset name 


%% build microscope object used in the RoSEO

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

visualizeBases(n1);
%% initialize a parallel pool

% p = gcp('nocreate');
% if isempty(p)
%     parpool('local', 32)
% end


%% simulate images

for mm = 1:length(ID)

SMLMName = [fileName,'data',num2str(ID(mm)),suffix_name,'.tif'];
bkg = load([fileName,'data',num2str(ID(mm)),'_bkg',suffix_name,'.mat']);
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
        SM_bkg(:,:,count) = SMLM_bkg_temp(:,:,ceil(jj/50));
    end
 
parfor kk = Nimg_perStack*(ii-1)+1:min(Nimg_perStack*ii,Nimg)
    count = rem(kk,Nimg_perStack);
    if count ==0
        count = Nimg_perStack;
    end
   
SMLM_img = ([SM_img(:,:,count)])*toPhoton;
SMLM_bkg1= [SM_bkg(:,:,count)]*toPhoton;
SMLM_bkg1 = [ones(imgSize,imgSize)*mean(SMLM_bkg1(:,ceil(1:imgSize/4)),'all'),...
             ones(imgSize,imgSize)*mean(SMLM_bkg1(:,ceil(1:imgSize/4)+imgSize),'all')];
SMLM_img(SMLM_img<0)=0;

[~,~,est] = RoSEO(n1, SMLM_img, SMLM_bkg1, FPSFx, FPSFy,'regval', reg);

SM_est_save{kk} = est;
end
save('est_temp.mat','SM_est_save');
end

mean_bkg = mean(mean(mean(mean(SM_bkg))))*toPhoton;

%% project the 6 moments estimation to 3 angulular space
% ***
Angle_save_cell={};
SM_est_save_all_cell = {};
Angle_save = [];
SM_est_save_all = [];

parfor kk = 1:Nimg
    SM_est = SM_est_save{kk};
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


save([num2str(ID(mm)) '_est.mat'],'SM_est_save','Angle_save','SM_est_save_all','left_to_right_trans_ratio','pmask','reg','note');
end
%%
