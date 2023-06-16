%% Monte-Carlo simulation conditions
x1 = linspace(-1,1,30);
x2 = linspace(-1,1,30);
[X,Y] = meshgrid(x1,x2);
indx  = (X.^2+Y.^2)>=1;
X(indx)=nan;
Y(indx)=nan;

mux_sim = 2*X.*sqrt(1-X.^2-Y.^2);
muy_sim = 2*Y.*sqrt(1-X.^2-Y.^2);
muz_sim = 1-2*(X.^2+Y.^2);
indx = muz_sim<0;
mux_sim(indx) = nan;
muy_sim(indx) = nan;
muz_sim(indx) = nan;


phi_sim = linspace(0,180,36); % in unit of degree
theta_sim = linspace(0,90,10); % in unit of degree
omega_sim = [0,1,2];
z_sim = 700; % in unit of nm
NFP = -580;  % in unit of nm
Nimg_perState = 200;
imgSize = 71;
signal = 2500;
background = 3;
%% parameters of the microscope

%----------------------------------------------
% create nanoscope object analysis

emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1;
zerothOrder = [0,0];%zerothOrder_RL; 
phasemaskpara.zeroorder = zerothOrder; %usually used for PSFs with large footprint, and some light leak into the center of PSF
maskName = fullfile('phasemask', 'pixOL_v12');
phasemaskpara.maskname = maskName;
toPhoton = 1; %gain of the camera electron/ADU



%% parameter of the data analyse process

imgPara.pix_sizex = 58.5; %the griding unit for RoSEO3D, in unit of nm
imgPara.pix_sizez = 50; % griding size along z direction
imgPara.axial_grid_points = [-6:1:28]*imgPara.pix_sizez;  % build the z slices 

scaling_factor = [100,100,100];   % used to balence the weights on orientation estimation and location estimation, experimental trier number
imgPara.scaling_factor = scaling_factor;

%--------------------------------------------------
reg = 0.22; % regulization for step 1 in RoSEO
note = ''; % for adding additional notes

%% build microscope object used in the RoSEO3D step2 and step3

n1 = Nanoscope('imageSize', imgSize*2+41,... 
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);


%% build microscope object used in the RoSEO3D step2 and step3

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


%% build the micrscope for the simulated z plane. for generate the SM images

[PSFx, PSFy] = n1.createPSFstruct3D(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',NFP,...
    'axial_grid_points',[z_sim,z_sim+10],...
    'scaling_factor',scaling_factor); %distance is is nm

Basis_pixOL_X = cat(4,PSFx.XXx,PSFx.YYx,PSFx.ZZx,PSFx.XYx,PSFx.XZx,PSFx.YZx);
Basis_pixOL_Y = cat(4,PSFy.XXy,PSFy.YYy,PSFy.ZZy,PSFy.XYy,PSFy.XZy,PSFy.YZy);
Basis_pixOL_X = squeeze(Basis_pixOL_X(:,:,1,:));
Basis_pixOL_Y = squeeze(Basis_pixOL_Y(:,:,1,:));
Basis_pixOL = reshape(cat(2, Basis_pixOL_X,Basis_pixOL_Y),[],6);
%% initialize a parallel pool

p = gcp('nocreate');
if isempty(p)
    parpool('local', 32)
end


%% simulate images

for ii = 1:size(mux_sim,1)%length(theta_sim)
    for jj = 1:size(mux_sim,2)
        
      for mm = 1:length(omega_sim)
        if isnan(mux_sim(ii,jj))==0
            mux_cur = mux_sim(ii,jj);
            muy_cur = muy_sim(ii,jj);
            muz_cur = muz_sim(ii,jj);
            theta_cur = acos(muz_cur)/pi*180;
            phi_cur = atan2(muy_cur,mux_cur)/pi*180;
            %phi_cur = phi_sim(jj);
            
            omega_cur = omega_sim(mm);
            gamma = 1-3*omega_cur/4/pi+omega_cur^2/8/pi^2;
            % generate image

            [muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_angleD_gamma_to_M(theta_cur,phi_cur,gamma);
            SM_img = reshape(signal*Basis_pixOL*[muxx,muyy,muzz,muxy,muxz,muyz].'+background,imgSize,imgSize*2);
            %SM_img = double(SM_img);
            SM_bkg = ones(imgSize,imgSize*2)*background;
            
            
            
        parfor kk = 1:Nimg_perState

        SMLM_img = poissrnd(SM_img);
        SMLM_bkg1= SM_bkg;

        [SM_est,NLL_cur_out,~,NLL_GT] = RoSEO3D(SM_img,SMLM_img, SMLM_bkg1,n1,FPSFx, FPSFy,imgPara,imgPara_B, 'regval',reg);
        
         if abs(SM_est(1))>100  
             NLL_cur_out_save(kk)=NLL_cur_out;
             NLL_GT_save(kk)=NLL_GT;            
         end
        SM_est_save{kk} = SM_est;
        NLL_save{kk}=NLL_cur_out;
        end


        mean_bkg = mean(mean(mean(mean(SM_bkg))))*toPhoton;

%% project the 6 moments estimation to 3 angulular space
% ***
Angle_save=[];
SM_est_save_all = [];
NLL_save_all = [];
for kk = 1:Nimg_perState
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
        else
            SM_est_save = {nan};
            NLL_save = {nan};
            Angle_save = {nan};
            SM_est_save_all = {nan};
            NLL_save_all = {nan};
        end

SM_est_save_map{ii,jj,mm}=SM_est_save;
NLL_save_map{ii,jj,mm}=NLL_save;
Angle_save_map{ii,jj,mm}=Angle_save;
SM_est_save_all_map{ii,jj,mm}=SM_est_save_all;
NLL_save_all_map{ii,jj,mm}=NLL_save_all;
       
%
        end
    end
end
%

save(['MC_at_z_' num2str(z_sim) '_SNR2500_3_uniformSample_PSFsizeThred1.mat'],'SM_est_save_map','NLL_save_map','Angle_save_map','SM_est_save_all_map','NLL_save_all_map');
