function [gammanew_out,loc_re_new_out,NLL_output] = FIST_optimize_v2(loc_re_init,SMLM_img_re,b,imgPara,imgParaB)
N = length(loc_re_init(:))/3;

%% step 1: psudo + small interation of FISTA to decide good initial points

%---------------------------FISTA optimize parameters---------------------------
MaxIt = 40;
Lmax = 0.0001;
rx = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;
rz = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;


gammaold_save = zeros(15,size(imgParaB.Bx,3));
NLL_cur_save = zeros(1,size(imgParaB.Bx,3));
loc_re_old_save = zeros(3,size(imgParaB.Bx,3));

for  ii= 1:length(imgPara.axial_grid_points)

loc_re_old = [0,0,imgPara.axial_grid_points(ii)].';
[G_re,~,~] = update_basisMatrix(N,zeros(15,1),loc_re_old,imgPara,imgParaB);
gammaold = (pinv(G_re)*(SMLM_img_re-b));
[G_re,loc_re_old,gammaold] = update_basisMatrix(N,gammaold,loc_re_old,imgPara,imgParaB);
gammaold = f_projection(gammaold,rx,rz,imgPara.scaling_factor);
%gammaold = gamma15Togamma9(gammaold,imgPara.scaling_factor);


[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,rx,rz,G_re,MaxIt,Lmax,imgPara,imgParaB); 
gammaold_save(:,ii) = gammanew_out;
loc_re_old_save(:,ii) = loc_re_new_out;
NLL_cur_save(ii) = NLL_cur_out;
end
%---------------------------choose the best initial points from candidates---------------------------

[~,indx] = min(NLL_cur_save);
gammaold = gammaold_save(:,indx);
loc_re_old = loc_re_old_save(:,indx);
%% upsampling
% img_size = sqrt(length(b)/2);
% b = reshape(repelem(reshape(b,img_size,img_size*2),2,2),[],1);
% SMLM_img_re = reshape(repelem(reshape(SMLM_img_re,img_size,img_size*2),2,2),[],1);
% imgPara.pix_sizex = imgPara.pix_sizex/2;
% 
% imgPara.Bx = imgPara.Bx1;
% imgPara.By = imgPara.By1;
% img_size =  img_size*2;
% 
% imgPara.Bx1 = [];
% imgPara.By1 = [];
% 
% imgPara.img_sizex = img_size;
% imgPara.img_sizey = img_size;

%% step 2: FISTA for ROI 

[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_old,imgPara,imgParaB);

%---------------------------FISTA optimize parameters---------------------------

MaxIt = 40;
Lmax = 0.0001;
rx = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;
rz = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;

[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,rx,rz,G_re,MaxIt,Lmax,imgPara,imgParaB);      

%---------------------------output---------------------------
loc_re_new_out = loc_re_new_out+loc_re_init;
NLL_output = NLL_cur_out;
end








