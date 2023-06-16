function [imgPara,imgPara_B] = f_generate_microscope(RoSEOPara,imgSize)
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
imgPara.pix_sizez = RoSEOPara.pix_sizez; %the griding unit for RoSEO, in unit of nm
imgPara.axial_grid_points = RoSEOPara.axial_grid_points;  % build the z slices  
imgPara.twoD_fitRange = RoSEOPara.twoD_fitRange; %
%reg = RoSEOPara.reg;
NFP = RoSEOPara.NFP;
imgPara.scaling_factor = [1,1,1];

%%
% build microscope object used in the RoSEO2D for finding the 
n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', toPhoton,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);

[FPSFx, FPSFy] = n1.createPSFstruct(n1,...
    'ytoxchanneltransratio', left_to_right_trans_ratio,...
    'normal_focal_plane',NFP,...
    'molecule_plane',imgPara.twoD_fitRange);   % [z min, z max]; the values are the expected z range for SMs, no need to be accurate

basis_matrix = [reshape([n1.XXxBasis,n1.XXyBasis],[],1),...
               reshape([n1.YYxBasis,n1.YYyBasis],[],1),...
               reshape([n1.ZZxBasis,n1.ZZyBasis],[],1),...
               reshape([n1.XYxBasis,n1.XYyBasis],[],1),...
               reshape([n1.XZxBasis,n1.XZyBasis],[],1),...
               reshape([n1.YZxBasis,n1.YZyBasis],[],1)];

h = visualizeBases(basis_matrix);


%% build basis image for RoSEO3D

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
    'scaling_factor',imgPara.scaling_factor); %distance is is nm

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
imgPara.img_sizex = imgSize;
imgPara.img_sizey = imgSize;

end


