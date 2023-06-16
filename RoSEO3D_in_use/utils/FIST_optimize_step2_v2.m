function [gammanew_out,loc_re_new_out,NLL_cur_out] = FIST_optimize_step2_v2(gamma_init,loc_re_init,SMLM_img_re,b,imgPara,imgParaB)


gamma_init = reshape(gamma_init,[],1);
N = length(gamma_init)/15;

%---------------------------FISTA optimize parameters---------------------------
MaxIt = 100;
gammaold = gamma_init;
Lmax = 0.0001;
loc_re_old = loc_re_init;
rx = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;
rz = imgPara.pix_sizex/imgPara.scaling_factor(1)/2*2;
[G_re,~,~] = update_basisMatrix(N,gammaold,loc_re_old,imgPara,imgParaB);

%--------------------------
[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA(SMLM_img_re,b,gammaold,loc_re_old,rx,rz,G_re,round(MaxIt/2),Lmax,imgPara,imgParaB);  

[G_re,~,~] = update_basisMatrix(N,gammanew_out,loc_re_new_out,imgPara,imgParaB);  
[gammanew_out,loc_re_new_out,NLL_cur_out] = FISTA_nonupdates(SMLM_img_re,b,gammanew_out,loc_re_new_out,rx*2,rz*2,G_re,round(MaxIt/2),Lmax,imgPara);  
 
end




