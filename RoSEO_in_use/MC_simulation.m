%% Necessary information

% place following files in a directory called source_dir
% 1- transmissionRatio_name: a mat file that contains transmit_ratio_L2R as
% a variable
% transmit_ratio_L2R is a scalar indicating the light transmission ratio
% between left and right channels
%
% 2- zerothOrder_name: a mat file containing zerothOrder as a variable
% zerothOrder_RL is a 1*2 array specifying the zeroth order factor
% in the pupil plane according to [x_channel y_channel] or [right channel left
% channel] formt
%
% 3- image stack and background stack in binary format
% You may want to write your image and background stack into binary format
% with names img_stack_name.bin and backg_stack_name.bin accordingly
% use writeSMLMbackg2bin function


%% format data and computational resources

%----------------------------------------------
frames_per_state = 1; % total number of frames in your stack

imgSize = 41;
background = 2;
signal = 1000;


v=linspace(0,0.5,50);
u=linspace(0.50,0.995,30);
phiD_simulate=2*pi*v/pi*180;
thetaD_simulate=acos(2*u-1)/pi*180;
omega_simulate = 0;
gamma_simulate = 1;

%%
% initialize a parallel pool
% p = gcp('nocreate');
% if isempty(p)
%     parpool('local', number_of_workers)
% end
% num_full_stack = floor(num_frames/num_frames_p);


%% analyze via RoSEO

%----------------------------------------------
% create nanoscope object analysis

maskName = fullfile('phasemask', 'pixOL_v12');
emitter_wavelength = 610; %nm
refractiveIndx=1.515;% objective refractive index
sampleRefractiveIndx=1.314;% sample refractive index

left_to_right_trans_ratio = 1;%transmit_ratio_L2R; for Di03-R488/561, 593/46; y channel.x channel
zerothOrder = [0,0];%zerothOrder_RL;
%construct phasemaskpara
phasemaskpara.zeroorder = zerothOrder;
phasemaskpara.maskname = maskName;


n1 = Nanoscope('imageSize', imgSize,...
    'ADcount', 1,...
    'emissWavelength', emitter_wavelength, ...
    'refractiveIndx',refractiveIndx,...
    'sampleRefractiveIndx',sampleRefractiveIndx,...
    'phasemaskpara', phasemaskpara);
% create PSF matrix accounting for channel transmission ratio
[FPSFx, FPSFy] = n1.createPSFstruct(n1, 'ytoxchanneltransratio', left_to_right_trans_ratio);


%% check basis
B = n1.computeBasis_static(n1);
B = B./(sum(B(1:3,:),'all')/3);

%%
RoSEO_h = @(img)RoSEO(n1, img, backg, FPSFx, FPSFy, 'regval', .22);

loc_data_tot = [];
est_all = [];
frame_count = 0;
for ii = 1:length(phiD_simulate)
for jj = 1:length(thetaD_simulate)

    phiD_cur = phiD_simulate(ii);
    thetaD_cur = thetaD_simulate(jj);
    gamma_cur = gamma_simulate; 
    [muxx,muyy,muzz,muxy,muxz,muyz] = Quickly_rotating_matrix_angle_gamma_to_M(thetaD_cur,phiD_cur,gamma_cur);
    M = [muxx,muyy,muzz,muxy,muxz,muyz];
    I_cur = signal*M*B+background;

    for kk = 1:frames_per_state
        I_possion = poissrnd(I_cur);
        I_possion = reshape(I_possion,imgSize,imgSize*2);
        [~,~,est] = RoSEO(n1, I_possion, I_possion*0+background, FPSFx, FPSFy);
        loc_data_tot{kk}=est;
    end


%% save results


for kk = 1:frames_per_state

    frame_count = frame_count+1;

    if ~isempty(loc_data_tot{kk})
        loc_data = loc_data_tot{ii,jj,kk};
        NSMs = size(loc_data,1);
        for ll=1:NSMs
            
            secM_cur(1:6)=loc_data(ll,5:10);
            signal = loc_data(ll,4);
            [mux,muy,muz,rotMobil] = secondM2SymmCone_for_experiment(secM_cur,loc_data(ll,4),background,B.');
            if muz<=0
                mux = -mux;
                muy = -muy;
                muz = -muz;
            end
            [thetaD, phiD, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
            omega = 2*pi*(3/2-sqrt(2*rotMobil+1/4));
            est_cur = [frame_count,signal,loc_data(ll,2:3),thetaD,phiD,rotMobil,loc_data(ll,5:10)];
            est_all = cat(1,est_all,est_cur);
        end
    end
end
end
end


save('RoSEO_MC_simulation_pixOL_SNR1000_2.mat','est_all');