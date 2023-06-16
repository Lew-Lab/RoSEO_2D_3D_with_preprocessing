
load('Micrscopy_inform.mat')
load('est_retrieval_1.1_v22.mat')
%%
clear sXX sYY sZZ sXY sXZ sYZ XX_dx YY_dx ZZ_dx XX_dy YY_dy ZZ_dy XX_dz YY_dz ZZ_dz SM_thetaD SM_phiD
x_cord = x; y_cord = y;z_cord = z;
%x_cord = x; y_cord = y;z_cord = z-57+14;  %-90

phiD_cord = rem(360-phiD+180,360); phiD_cord(phiD_cord>180)=phiD_cord(phiD_cord>180)-360;
r=100;
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<1100 & z_cord>000 &SM_est_final(:,5)>400;

ii=0;
for jj=1:size(SM_est_final,1)
  
    if indx(jj)==1
        ii = ii+1;
SM_est_cur = SM_est_final(ii,:);

[gamma,loc] = SM_est2gamma(SM_est_cur,imgPara);
signal(ii) = sum(gamma(1:6));
sXX(ii) = gamma(1);
sYY(ii) = gamma(2);
sZZ(ii) = gamma(3);
sXY(ii) = gamma(4);
sXZ(ii) = gamma(5);
sYZ(ii) = gamma(6);

XX_dx(ii) = gamma(7)/sXX(ii)*100;
YY_dx(ii) = gamma(8)/sYY(ii)*100;
ZZ_dx(ii) = gamma(9)/sZZ(ii)*100;

XX_dy(ii) = gamma(10)/sXX(ii)*100;
YY_dy(ii) = gamma(11)/sYY(ii)*100;
ZZ_dy(ii) = gamma(12)/sZZ(ii)*100;

XX_dz(ii) = gamma(13)/sXX(ii)*100;
YY_dz(ii) = gamma(14)/sYY(ii)*100;
ZZ_dz(ii) = gamma(15)/sZZ(ii)*100;

SM_phiD(ii) = phiD(ii);
SM_thetaD(ii) = thetaD(ii);

    end
end

%%
bins = 20;
figure();
subplot(3,3,1);
histogram(XX_dx,bins ); title('XX dx');

subplot(3,3,2);
histogram(YY_dx,bins); title('YY dx');

subplot(3,3,3);
histogram(ZZ_dx,bins); title('ZZ dx');

subplot(3,3,4);
histogram(XX_dy,bins); title('XX dy');

subplot(3,3,5);
histogram(YY_dy,bins); title('YY dy');

subplot(3,3,6);
histogram(ZZ_dy,bins); title('ZZ dy');

subplot(3,3,7);
histogram(XX_dz,bins); title('XX dz');

subplot(3,3,8);
histogram(YY_dz,bins); title('YY dz');

subplot(3,3,9);
histogram(ZZ_dz,bins); title('ZZ dz');

whitebg('w');
set(gcf,'Color','w');

%%
dens = 0.2;
figure();
subplot(3,3,1);
scatter(XX_dx,sXX,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dx'); ylabel('XX signal'); ylim([0,2000]); 

subplot(3,3,2);
scatter(YY_dx,sYY,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dx');  ylabel('YY signal'); ylim([0,2000]);

subplot(3,3,3);
scatter(ZZ_dx,sZZ,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dx'); ylabel('ZZ signal'); ylim([0,2000]);

subplot(3,3,4);
scatter(XX_dy,sXX,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dy'); ylabel('XX signal'); ylim([0,2000]);

subplot(3,3,5);
scatter(YY_dy,sYY,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dy'); ylabel('YY signal'); ylim([0,2000]);

subplot(3,3,6);
scatter(ZZ_dy,sZZ,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dy'); ylabel('ZZ signal'); ylim([0,2000]);

subplot(3,3,7);
scatter(XX_dz,sXX,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dz'); ylabel('XX signal'); ylim([0,2000]);

subplot(3,3,8);
scatter(YY_dz,sYY,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dz'); ylabel('YY signal'); ylim([0,2000]);

subplot(3,3,9);
scatter(ZZ_dz,sZZ,2,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dz'); ylabel('ZZ signal'); ylim([0,2000]);

whitebg('w');
set(gcf,'Color','w');

%%
y = SM_thetaD; ylim_choice = [0,90]; sizeS = 2; ylabelName = 'thetaD';
dens = 0.2;
figure();
subplot(3,3,1);
scatter(XX_dx,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dx'); ylabel(ylabelName); ylim(ylim_choice); 

subplot(3,3,2);
scatter(YY_dx,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dx');  ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,3);
scatter(ZZ_dx,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dx'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,4);
scatter(XX_dy,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dy'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,5);
scatter(YY_dy,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dy'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,6);
scatter(ZZ_dy,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dy'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,7);
scatter(XX_dz,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('XX dz'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,8);
scatter(YY_dz,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('YY dz'); ylabel(ylabelName); ylim(ylim_choice);

subplot(3,3,9);
scatter(ZZ_dz,y, sizeS,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); xlabel('ZZ dz'); ylabel(ylabelName); ylim(ylim_choice);

whitebg('w');
set(gcf,'Color','w');


%%
figure();
subplot(2,3,1);
histogram(SM_est_final(:,end-5)); title('XX');

subplot(2,3,2);
histogram(SM_est_final(:,end-4)); title('YY');

subplot(2,3,3);
histogram(SM_est_final(:,end-3)); title('ZZ');

subplot(2,3,4);
histogram(SM_est_final(:,end-2)); title('XY');

subplot(2,3,5);
histogram(SM_est_final(:,end-1)); title('XZ');

subplot(2,3,6);
histogram(SM_est_final(:,end)); title('YZ');
