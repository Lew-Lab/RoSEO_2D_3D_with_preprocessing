
%%
Angle_save=[];
SM_est_save_all = [];
for kk = 1:2000 
    SM_est = SM_est_save{kk};
if isempty(SM_est)==0
    SM_est(SM_est(:,4)<500,:) = [];
    
    %
    saveAngle = [];
    saveOren = [];
    for ll = 1:size(SM_est,1)
        [mux,muy,muz,rotMobil] = secondM2SymmCone_RoSEO3D(double(SM_est(ll,:)),mean(mean(backg)),imgPara);
        if muz<=0
            mux = -mux;
            muy = -muy;
            muz = -muz;
        end
        saveOren(ll,:) = [mux,muy,muz,rotMobil];
        [thetaD, phiD_cord, alphaD] = symmCone2angle(mux,muy,muz,rotMobil);
        saveAngle = [thetaD,phiD_cord,alphaD,rotMobil,3*pi-sqrt(rotMobil*8*pi^2+pi^2)];

        Angle_save = [Angle_save;kk,saveAngle];
        SM_est_save_all = [SM_est_save_all;kk,double(SM_est(ll,:))];
    end
end


end

%%

indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x1 = SM_est_save_all(indx,2);
y1 = SM_est_save_all(indx,3);
z1 = SM_est_save_all(indx,4);
[R,~]=sphereFit([x1,y1,z1])
R=[0,0,0];
%
indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x = SM_est_save_all(:,2)-R(1);
y = SM_est_save_all(:,3)-R(2);
z = SM_est_save_all(:,4);
thetaD = Angle_save(:,2);
phiD = Angle_save(:,3);
omega = Angle_save(:,6);
Angle_save_final = Angle_save;
SM_est_final = SM_est_save_all;
gamma_est_final = gamma_est_all;
%




%save('est_retrieval_1.2_v2.mat','x','y','z','thetaD','phiD','omega','Angle_save_final','SM_est_final');
%%



indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x1 = SM_est_save_all(indx,2);
y1 = SM_est_save_all(indx,3);
z1 = SM_est_save_all(indx,4);
[R,~]=sphereFit([x1,y1,z1])
R = [0,0,0];


r=30;
indx = abs(SM_est_save_all(:,2))<1300  & abs(SM_est_save_all(:,3))<1300 & SM_est_save_all(:,4)>100 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
indx(:)=1;
x = [x;SM_est_save_all(:,2)-R(1)];
y = [y;SM_est_save_all(:,3)-R(2)];
z = [z;SM_est_save_all(:,4)];
thetaD = [thetaD;Angle_save(:,2)];
phiD = [phiD;Angle_save(:,3)];
omega = [omega;Angle_save(:,6)];
Angle_save_final = [Angle_save_final;Angle_save];
SM_est_final = [SM_est_final;SM_est_save_all];
gamma_est_final = [gamma_est_final;gamma_est_all];
%%
%save('est_v1.mat','x','y','z','thetaD','phiD','omega','Angle_save_final','SM_est_final');
%%
load('colorSpace.mat');

%%
figure();
histogram(SM_est_final(:,5));
xlabel('photons','FontSize',10); ylabel('counts','FontSize',10);
title({['mean:' num2str(mean(SM_est_final(:,5)))],['median:' num2str(median(SM_est_final(:,5)))]});
ax = gca;
ax.FontSize = 10; 

%%
theta_view = linspace(0,20,100);
phi_view = linspace(0,360,100);
%theta_view = 90;
%phi_view = 0;
%x_cord = x; y_cord = y;z_cord = z; 
%phiD_cord = phiD;
%x_cord = x-58; y_cord = y-30;z_cord = z+1000; 
%x_cord = x-175-6+19; y_cord = y+80-46+6;z_cord = z-55; 
%x_cord = x-43; y_cord = y+113;z_cord = z-57+14;  %-90
%x_cord = x-43-36; y_cord = y+113+36;z_cord = z-22;
%x_cord = x; y_cord = y;z_cord = z-57+14;  %-90
x_cord = x; y_cord = y;z_cord = z;
phiD_cord = phiD;
%phiD_cord = rem(360-phiD+180,360); phiD_cord(phiD_cord>180)=phiD_cord(phiD_cord>180)-360;
r=100;
indx = abs(x_cord)<1200 & abs(y_cord)<1200 & z_cord<1200 & z_cord>-100 &SM_est_final(:,5)>400;

[R,C]=sphereFit([x_cord(indx), y_cord(indx),z_cord(indx)])

space = round(linspace(1,6380,4500));
%indx(space)=0;
%indx(1000:end)=0;
%indx = abs(x)<7000000;
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<1000 & z>-200;
%figure('visible','off');
figure('Color',[0.1,0.1,0.1]);
hold on
scatter3(x_cord(indx),y_cord(indx),z_cord(indx),20,thetaD(indx),'filled','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7); axis image; 
%colormap(squeeze(colorSpace));
colormap(turbo);
color = squeeze(colorSpace);
thetaD_dic = linspace(0,90,256);
phiD_dic = linspace(-180,180,403);
%r = 40;
%caxis([-180,180]);
% hold on;
%colorbar;
% for ii=1:length(indx)
%     if indx(ii)==1
%         [~,color_indx] = min(abs(phiD_cord(ii)-phiD_dic));
%         plot3([x_cord(ii)-r*sind(thetaD(ii)).*cosd(phiD_cord(ii)),x_cord(ii),x_cord(ii)+r*sind(thetaD(ii)).*cosd(phiD_cord(ii))].',...
%              [y_cord(ii)-r*sind(thetaD(ii)).*sind(phiD_cord(ii)),y_cord(ii),y_cord(ii)+r*sind(thetaD(ii)).*sind(phiD_cord(ii))].',...
%              [z_cord(ii)+r*cosd(thetaD(ii)),z_cord(ii),z_cord(ii)-r*cosd(thetaD(ii))].',...
%              'Color',color(color_indx,:),'LineWidth',1);
%          %thetaD1 = real(acos((1000-(z_cord))./1000)/pi*180); phiD1 = atan2(y_cord,x_cord)/pi*180; 
%     end
% end
% axis image
%phiD1(thetaD1>90)=phiD1(thetaD1>90)+180; thetaD1(thetaD1>90)=180-thetaD1(thetaD1>90);
%phiD1 = rem(phiD1,360); phiD1(phiD1>180) = -360+phiD1(phiD1>180);
%plot3([x(indx)-r*sind(thetaD1(indx)).*cosd(phiD1(indx)),x(indx),x(indx)+r*sind(thetaD1(indx)).*cosd(phiD1(indx))].',...
%     [y(indx)-r*sind(thetaD1(indx)).*sind(phiD1(indx)),y(indx),y(indx)+r*sind(thetaD1(indx)).*sind(phiD1(indx))].',...
%     [z(indx)+r*cosd(thetaD1(indx)),z(indx),z(indx)-r*cosd(thetaD1(indx))].','r');
 
 %
R=350;
[x_sphere,y_sphere,z_sphere] = sphere;
%surf(x_sphere*1000-120,y_sphere*1000-180,z_sphere*1000+800);
%viscircles([-120,-180],580)
%xlim([100,300]);
xlabel('x(nm)'); ylabel('y(nm)'); zlabel('z(nm)');

%hcb=colorbar; hcb.Title.String = "\phi(\circ)";
ax = gca;
ax.FontSize = 10; 
whitebg('black');
set(gcf, 'InvertHardcopy', 'off');
% for ii=1:length(theta_view)
% view(phi_view(ii),theta_view(ii));
% axis vis3d;
% saveas(gcf,['C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210414\spikes_ROIxyz_bkg\',num2str(ii),'.png'])
% end
%
[thetaD1,phiD1]=sphereLine(x_cord,y_cord,z_cord,1000);
%thetaD1 = real(acos((1000-(z_cord))./1000)/pi*180); phiD1 = atan2(y_cord,x_cord)/pi*180; 
%phiD1(thetaD1>90)=phiD1(thetaD1>90)+180; thetaD1(thetaD1>90)=180-thetaD1(thetaD1>90);
%phiD1 = rem(phiD1,360); phiD1(phiD1>180) = -360+phiD1(phiD1>180);

%% z slice plot
zSlice = linspace(000,1200,13);
%zSlice = linspace(-100,600,8);
%zSlice = linspace(600,1400,9);
%figure('Color',[0.1,0.1,0.1]);

Fig = figure('Color',[0.1,0.1,0.1]);
set(Fig, 'Units','centimeters','InnerPosition', [10 10 35 12]*1.0);
x_offset = [0,sqrt(1000^2-(1000-(zSlice(2:end))))*1.2];


load('pixOL_loc_precision_NFP-700.mat');
sigma_r = sqrt(sigma_x_opt.^2+sigma_y_opt.^2);
for ii=1:2:length(zSlice)-1
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) &SM_est_final(:,5)>400;
space = round(linspace(1,6380,4500));
%indx(space)=0;
%indx(1000:end)=0;
%indx = abs(x)<7000000;
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<1000 & z>-200;
%figure('visible','off');
%ax1 = axes('Position',[0.1/6+0.5/6*(2*ii-2) 0.6 0.9/6 .4]);
hold on
%scatter(x_cord(indx),y_cord(indx),30,phiD_cord(indx),'filled','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7); axis image; 
%colormap(squeeze(colorSpace));


color = squeeze(colorSpace);
%color = squeeze(turbo(403));
thetaD_dic = linspace(0,90,256);
phiD_dic = linspace(-180,180,403);
r = 40;
%caxis([-180,180]);
hold on;
%colorbar;
for jj=1:length(indx)
    if indx(jj)==1
        [~,color_indx] = min(abs(phiD_cord(jj)-phiD_dic));
        plot(sum(x_offset(1:ii+1))+[x_cord(jj)-r*sind(thetaD(jj)).*cosd(phiD_cord(jj)),x_cord(jj),x_cord(jj)+r*sind(thetaD(jj)).*cosd(phiD_cord(jj))].',...
             [y_cord(jj)-r*sind(thetaD(jj)).*sind(phiD_cord(jj)),y_cord(jj),y_cord(jj)+r*sind(thetaD(jj)).*sind(phiD_cord(jj))].',...
             'Color',color(color_indx,:),'LineWidth',1);
         %thetaD1 = real(acos((1000-(z_cord))./1000)/pi*180); phiD1 = atan2(y_cord,x_cord)/pi*180; 
    end
end
axis image
%xlim([-1300,1300]);ylim([-1300,1300]);
%---- plot reference sphere
R_cur = sqrt(1000^2-(1000- zSlice(ii)-50).^2);
R_top = sqrt(1000^2-(1000- zSlice(ii)-100).^2);
R_bottom = sqrt(1000^2-(1000- zSlice(ii)).^2);
if zSlice(ii)+50<=0
    R_cur=0;
    R_top=0;
    R_bottom=0;
end

[~,indx_top]=min(abs(zSlice(ii+1)-z_range*10^9));
sigma_top = sigma_r(indx_top)*10^9;

[~,indx_bottom]=min(abs(zSlice(ii)-z_range*10^9));
sigma_bottom = sigma_r(indx_bottom)*10^9;

R_sigma_top = R_top+sigma_top;
R_sigma_bottom = R_bottom-sigma_bottom; 
if R_sigma_bottom<=0
    R_sigma_bottom=0;
end
%viscircles([sum(x_offset(1:ii+1)),0]+[0,0],R_cur,'LineStyle','--','LineWidth',0.1);
%viscircles([sum(x_offset(1:ii+1)),0]+[0,0],R_sigma_top,'LineStyle','-','LineWidth',0.5);
%viscircles([sum(x_offset(1:ii+1)),0]+[0,0],R_sigma_bottom,'LineStyle','-','LineWidth',0.5);
R=350;
[x_sphere,y_sphere,z_sphere] = sphere;
%---
whitebg('k');
%colorbar('Color','w');
%xlabel('x','Color','w'); ylabel('y','Color','w'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
ax = gca;
ax.FontSize = 10; 
axis off;
whitebg('black');


set(gcf,'Color',[0 0 0]);
text(-700+sum(x_offset(1:ii+1)),1300,[num2str(zSlice(ii)),'nm<z <',num2str(zSlice(ii+1)),'nm'],'Color','w','FontSize',8);
set(gcf, 'InvertHardcopy', 'off');


%view(phi_view(ii),theta_view(ii));
% for kk=1:10
% saveas(gcf,['C:\Users\wu.t\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210509\3-17\soikes_calibrated_v3\',num2str((ii)),'.png'])
% end

end


%% z slice plot
zSlice = linspace(0,1200,13);
%zSlice = linspace(-100,600,8);
%zSlice = linspace(600,1000,5);
%figure('Color',[0.1,0.1,0.1]);


figure();
load('pixOL_loc_precision_NFP-700.mat');
sigma_r = sqrt(sigma_x_opt.^2+sigma_y_opt.^2);
for ii=1:2:length(zSlice)-1
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) &SM_est_final(:,5)>1000;
space = round(linspace(1,6380,4500));

%xlim([-1300,1300]);ylim([-1300,1300]);
%---- plot reference sphere
R_cur = sqrt(1000^2-(1000- zSlice(ii)-50).^2);
R_top = sqrt(1000^2-(1000- zSlice(ii)-100).^2);
R_bottom = sqrt(1000^2-(1000- zSlice(ii)).^2);
if zSlice(ii)+50<=0
    R_cur=0;
    R_top=0;
    R_bottom=0;
end
[~,indx_top]=min(abs(zSlice(ii+1)-z_range*10^9));
sigma_top = sigma_r(indx_top)*10^9;

[~,indx_bottom]=min(abs(zSlice(ii)-z_range*10^9));
sigma_bottom = sigma_r(indx_bottom)*10^9;

subplot(1,(length(zSlice)-1)/2,(ii+1)/2); hold on;

dD = sqrt((x_cord(indx).^2+y_cord(indx).^2))-R_cur;
% pd = fitdist(dD,'Normal');
% x=-200:1:200;
% p = normpdf(x,pd.mu,pd.sigma)*sum(indx)*400;
% plot(x,p); hold on;
edges=-400:30:400;
h = histogram(dD,edges);
h.FaceColor = [128, 255, 0]/255;
h.EdgeColor = 'none';
% h = histfit(dD);
% h(1).FaceColor = 'none';
% h(1).EdgeColor = 'none';

x_ax=-400:1:400;
gauss_pdf = normpdf(x_ax,0,(sigma_top+sigma_bottom)/2);
rec_pdf = zeros(size(x_ax));
rec_pdf(x_ax>(R_bottom-R_cur)&x_ax<(R_top-R_cur)) = 1; 

if sum(rec_pdf)==0
   pdf = gauss_pdf; pdf = pdf/sum(pdf);
else
rec_pdf = rec_pdf/sum(rec_pdf);
pdf = conv(rec_pdf,gauss_pdf,'same'); pdf = pdf/sum(pdf);
end




plot(x_ax,pdf*sum(indx)*30,'Color',[255 255 153]/255, 'LineWidth',1.5);
xlabel('r (nm)'); ylabel('count');
ax = gca;
ax.FontSize = 8; 

end


%%
figure('Color','k');
R = 0.015;
N = 16;
% because there are 31 points through the whole circle, 16 points at one direction, 16 points at opposite direction but one point is common (the middle point)
[rho,th] = meshgrid(linspace(0,R,N),linspace(0,2*pi,403));
xp = rho.*cos(th);
yp = rho.*sin(th);
Tp = squeeze(colorSpace);
pcolor(xp,yp,th/pi*180-180); hold on;
scatter(0,0,100,[0,0,0],'filled');
shading('interp');
%colorbar
colormap(Tp)
axis off; axis image
set(gcf, 'InvertHardcopy', 'off');
%% theta
figure('Color','w')
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),z_cord(indx),4,thetaD(indx),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; 
viscircles([0,1000],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('x(nm)'); ylabel('z(nm)'); 
ylim([-100,1000])
hcb=colorbar; hcb.Title.String = "\theta(\circ)";
colormap(turbo);
%axis off;
whitebg('k');
colorbar('Color','k');
xlabel('x','Color','k'); ylabel('z','Color','k'); 
%axis off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');

figure('Color','w')
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(y_cord(indx),z_cord(indx),4,thetaD(indx),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; 
viscircles([0,1000],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('y(nm)'); ylabel('z(nm)'); 
ylim([-100,1000])
hcb=colorbar; hcb.Title.String = "\theta(\circ)";
colormap(turbo);
%axis off;
whitebg('k');
colorbar('Color','k');
xlabel('y','Color','k'); ylabel('z','Color','k'); 
%axis off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');

figure('Color',[0,0,0])
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),y_cord(indx),4,thetaD(indx),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image;
viscircles([0,0],1000)
ax = gca;
ax.FontSize = 10; 

hcb=colorbar; hcb.Title.String = "\theta(\circ)";
colormap(turbo);
whitebg('k');
colorbar('Color','k');
xlabel('x','Color','k'); ylabel('y','Color','k'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');
%axis off;
% [X,Y,Z] = sphere;
% X(Z>0)=nan;Y(Z>0)=nan;Z(Z>0)=nan;
% X = X*R+20;Y = Y*R+30;Z=Z*R+R;
% surf(X,Y,Z,'FaceAlpha',0.5); %colormap('hot')

%% phi

figure('Color',[0,0,0])
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),z_cord(indx),8,phiD_cord(indx),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; colorbar;
viscircles([0,1000],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('x(nm)'); ylabel('z(nm)'); 
ylim([-100,1000])
hcb=colorbar; hcb.Title.String = "\phi(\circ)";
colorbar('Color','k');
xlabel('x','Color','k'); ylabel('z','Color','k'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');

colormap(squeeze(colorSpace));
%colormap('turbo');
figure('Color',[0,0,0])
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),y_cord(indx),8,phiD_cord(indx),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; colorbar;
viscircles([0,0],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('x(nm)'); ylabel('y(nm)'); 
hcb=colorbar; hcb.Title.String = "\phi(\circ)";
colormap(squeeze(colorSpace));
whitebg('k');
colorbar('Color','k');
xlabel('x','Color','k'); ylabel('y','Color','k'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');
%colormap('turbo');

%% omega

figure();
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),z_cord(indx),8,omega(indx)/pi,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; colorbar;
viscircles([0,1000],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('x(nm)'); ylabel('z(nm)'); 
ylim([-100,1000])
hcb=colorbar; hcb.Title.String = "\Omega(\pi)";
whitebg('k');
colorbar('Color','k');
xlabel('x','Color','k'); ylabel('z','Color','k'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');

figure();
%indx = x>-1400 & x<1140 & y>-1700 & y<950 & z<700 & z>-200;
scatter(x_cord(indx),y_cord(indx),8,omega(indx)/pi,'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); axis image; colorbar;
viscircles([0,0],1000)
ax = gca;
ax.FontSize = 10; 
xlabel('x(nm)'); ylabel('y(nm)'); 
hcb=colorbar; hcb.Title.String = "\Omega(\pi)";
whitebg('k');
colorbar('Color','w');
xlabel('x','Color','w'); ylabel('z','Color','w'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gcf,'Color','none');
%% normal theta
figure();
vector_est = [sind(thetaD(indx)).*cosd(phiD_cord(indx)),sind(thetaD(indx)).*sind(phiD_cord(indx)),cosd(thetaD(indx))];
vector_per = real([sind(thetaD1(indx)).*cosd(phiD1(indx)),sind(thetaD1(indx)).*sind(phiD1(indx)),cosd(thetaD1(indx))]);
angleD = real(acos(sum((vector_est.*vector_per),2))./pi*180);
angleD(angleD>90)=180-angleD(angleD>90);
histogram(angleD,50); %title('normal \theta');
ax = gca;
ax.FontSize = 10; 
whitebg('w');
xlabel('normal \theta(\circ)'); 
ylabel('count');
set(gcf,'Color','w');

%% normal theta
figure();
vector_est = [-sind(thetaD(indx)).*cosd(phiD_cord(indx)),sind(thetaD(indx)).*sind(phiD_cord(indx)),cosd(thetaD(indx))];
vector_perz = real([-sind(thetaD1(indx)).*cosd(phiD1(indx)),sind(thetaD1(indx)).*sind(phiD1(indx)),cosd(thetaD1(indx))]);
vector_perx = real([cos(thetaD1(indx)).*cosd(phiD1(indx)),-cosd(thetaD1(indx)).*sind(phiD1(indx)),sind(thetaD1(indx))]);
vector_pery = real([sind(phiD1(indx)), cosd(phiD1(indx)),phiD1(indx)*0]);
x_per = sum((vector_est.*vector_perx),2);
y_per = sum((vector_est.*vector_pery),2);
z_per = sum((vector_est.*vector_perz),2);
x_per(z_per<0)=-x_per(z_per<0);
y_per(z_per<0)=-y_per(z_per<0);
z_per(z_per<0)=-z_per(z_per<0);

normal_theta = real(acos(z_per))/pi*180;
normal_phi = real(atan2(y_per,x_per))/pi*180;
figure();
histogram(normal_theta,50); %title('normal \theta');
ax = gca;
ax.FontSize = 10; 
whitebg('w');
%xlabel('normal \theta(\circ)'); 
ylabel('count');
%set(gcf,'Color','none');
set(gcf,'Color','w');

figure();
histogram(normal_phi,50); %title('normal \theta');
ax = gca;
ax.FontSize = 10; 
whitebg('w');
xlabel('normal \phi(\circ)'); 
ylabel('count');
set(gcf,'Color','w');
%%
figure();
histogram(omega(indx)/pi,30); title('\Omega (\pi)');
xlim([0,2]);
ax = gca;
ax.FontSize = 10; 
%xlabel('\Omega (\pi)'); 
ylabel('count');
set(gcf,'Color','w');
%%
dens = 0.015;
figure();
scatter(angleD,omega(indx)/pi,8,'filled','MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);

hold on;
scatter(median(angleD),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0, 0, 1],'LineWidth',1);

ax = gca;
ax.FontSize = 10; 
box on;
whitebg('w');
set(gcf,'Color','w');
ylabel('\Omega(\pi)'); xlabel(' \theta perpendicular (\circ)');
%%
dens = 0.015;
hold on;
scatter(angleD,omega(indx)/pi,8,'filled','MarkerFaceColor',	[0.9290, 0.6940, 0.1250],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);
ylabel('\Omega(\pi)'); xlabel(' \theta perpendicular (\circ)');
hold on;
scatter(median(angleD),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'LineWidth',1);
ax = gca;
ax.FontSize = 10; 
box on;
whitebg('w');
set(gcf,'Color','w');

%% bind the data
theta_bin_edges = linspace(0,90,40);
omega_bin_edges = linspace(0,2,40)*pi;
hotN = bin2Ddata(omega(indx),angleD,omega_bin_edges,theta_bin_edges);
figure();
imagesc(theta_bin_edges,omega_bin_edges/pi,hotN);   colormap('hot');
ax = gca;
ax.YDir = 'normal'; %caxis([0,70]);

%% bind the data
phi_bin_edges = linspace(-180,180,40);
omega_bin_edges = linspace(0,2,40)*pi;
hotN = bin2Ddata(omega(indx),phiD_cord(indx),omega_bin_edges,phi_bin_edges);
figure();
imagesc(phi_bin_edges,omega_bin_edges/pi,hotN);  colormap('hot');  
ax = gca;
ax.YDir = 'normal';
%%
dens = 0.1;
figure();
subplot(1,3,1); 
scatter(thetaD(indx),z_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('z(nm)'); xlabel('\theta(\circ)');
%hold on; scatter(thetaD1(indx),z(indx),2,'filled','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); ylabel('z(nm)'); xlabel('theta(\circ)');
ax = gca;
ax.FontSize = 10;  box on;
subplot(1,3,2);
scatter(phiD_cord(indx),x_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('x(nm)'); xlabel('\phi(\circ)');
%hold on; scatter(phiD1(indx),x(indx),2,'filled','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
ax = gca;
ax.FontSize = 10;  box on;
subplot(1,3,3); 
scatter(phiD_cord(indx),y_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('y(nm)'); xlabel('\phi(\circ)');
%hold on; scatter(phiD1(indx),y(indx),2,'filled','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
ax = gca;
ax.FontSize = 10; box on;
whitebg('w');
set(gcf,'Color','w');

figure();
subplot(1,3,1); 
scatter(thetaD(indx),z_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('z(nm)'); xlabel('\theta(\circ)');
hold on; scatter(thetaD1(indx),z_cord(indx),2,'filled','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); ylabel('z(nm)'); xlabel('\theta(\circ)');
ax = gca;
ax.FontSize = 10;  box on;
subplot(1,3,2); 
scatter(phiD_cord(indx),x_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('x(nm)'); xlabel('\phi(\circ)');
hold on; scatter(phiD1(indx),x_cord(indx),2,'filled','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
ax = gca;
ax.FontSize = 10; box on;
subplot(1,3,3); 
scatter(phiD_cord(indx),y_cord(indx),8,omega(indx),'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('y(nm)'); xlabel('\phi(\circ)');
hold on; scatter(phiD1(indx),y_cord(indx),2,'filled','r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
ax = gca; 
ax.FontSize = 10; 
box on;
set(gcf,'Color','w');

%%
figure();
histogram(phiD_cord(indx),150); xlabel('estimated \phi(\circ)'); ylabel('count');
whitebg('w');
set(gcf,'Color','w');

%%
figure();
histogram(y_cord(indx),400); xlabel('estimated \phi(\circ)'); ylabel('count');
%set(gcf,'Color','w');
%%
figure();
scatter(phiD_cord(indx),omega(indx)/pi,8,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('omega'); xlabel('\phi(\circ)');

xlim([-180,180]);
%%

figure(); histogram(SM_est_final(indx,6),40); xlabel('XX');
figure(); histogram(SM_est_final(indx,7),40); xlabel('YY');
figure(); histogram(SM_est_final(indx,8),40); xlabel('ZZ');
figure(); histogram(SM_est_final(indx,9),40); xlabel('XY');
figure(); histogram(SM_est_final(indx,10),40); xlabel('XZ');
figure(); histogram(SM_est_final(indx,11),40); xlabel('YZ');
%%
dens = 0.2;
figure();
scatter(z_cord(indx),omega(indx)/pi,8,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('omega'); xlabel('z');
xlim([0,800]);

%%
figure();
scatter(z_cord(indx),phiD_cord(indx),8,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('phiD'); xlabel('z');
xlim([0,800]);

%%

figure();
scatter(thetaD(indx),z_cord(indx),8,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens); ylabel('z(nm)'); xlabel('\theta(\circ)');

xlim([0,90]);

%%

indx = omega<0.5;

figure(); histogram(z_cord(indx));

figure(); histogram(GT_final(indx,3));

figure(); histogram(z_cord(indx)-GT_final(indx,3));
%%
indx = omega>0.4 & z_cord>300 & z_cord<600;

figure(); 
titleString={'XX','YY','ZZ','XY','XZ','YZ'};

for ii = 1:6
    
subplot(2,3,ii);
histogram(SM_est_final(indx,5+ii)-GT_final(indx,5+ii),30);
xlabel(titleString{ii});
set(gcf,'Color','w');
whitebg('w');
end

%%
gamma_new = gamma2M_dx(gamma_est_final(:,2:16));
loc_int = round(GT_final(:,1:3)./[58.6,59.6,50]).*[58.6,59.6,50];
XXxyz = gamma_new(:,[8,11,14])+loc_int;
YYxyz = gamma_new(:,[8,11,14]+1)+loc_int;
ZZxyz = gamma_new(:,[8,11,14]+2)+loc_int;

loc_GT = GT_final(:,1:3);


XXdxyz = XXxyz-loc_GT;
YYdxyz = YYxyz-loc_GT;
ZZdxyz = ZZxyz-loc_GT;

%%
indx = omega>0.4 & z_cord>300 & z_cord<600;
figure(); 
titleString={'XXdx','XXdy','XXdz';'YYdx','YYdy','YYdz';'ZZdx','ZZdy','ZZdz';};
data_for_plot = cat(3,XXdxyz,YYdxyz,ZZdxyz);

count = 0;
for ii = 1:3
    for jj = 1:3
        count = count+1;
        subplot(3,3,count);
        histogram(data_for_plot(indx,jj,ii),30);
        xlabel(titleString{ii,jj});
        set(gcf,'Color','w');
        whitebg('w');

    end
end


%%
indx = omega>0.4 & z_cord>300 & z_cord<600;
figure(); 
titleString={'dx','dy','dz';};
dx = SM_est_final(:,2)-GT_final(:,1);
dy = SM_est_final(:,3)-GT_final(:,2);
dz = SM_est_final(:,4)-GT_final(:,3);

data_for_plot = [dx,dy,dz];
count = 0;
for ii = 1:3
 
        subplot(1,3,ii);
        histogram(data_for_plot(indx,ii),30);
        xlabel(titleString{ii});
        set(gcf,'Color','w');
        whitebg('w');

end

%%
indx = omega>0.4 & z_cord>300 & z_cord>600;
figure(); 
titleString={'XXdx','XXdy','XXdz';'YYdx','YYdy','YYdz';'ZZdx','ZZdy','ZZdz';};
data_for_plot = cat(3,XXdxyz-[dx,dy,dz],YYdxyz-[dx,dy,dz],ZZdxyz-[dx,dy,dz]);

count = 0;
for ii = 1:3
    for jj = 1:3
        count = count+1;
        subplot(3,3,count);
        histogram(data_for_plot(indx,jj,ii),30);
        xlabel(titleString{ii,jj});
        set(gcf,'Color','w');
        whitebg('w');

    end
end

%%
dens = 0.03;
indx = omega>0.4 & z_cord>300 & z_cord<600;
figure();

titleString={'XXdx','XXdy','XXdz';'YYdx','YYdy','YYdz';'ZZdx','ZZdy','ZZdz';};
data_for_plot = cat(3,XXdxyz,YYdxyz,ZZdxyz);

count = 0;
for ii = 1:3
    for jj = 1:3
        count = count+1;
        subplot(3,3,count);
        scatter(data_for_plot(indx,jj,ii),SM_est_final(indx,5+5)-GT_final(indx,5+5),10,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);
        xlabel(titleString{ii,jj});
        set(gcf,'Color','w');
        whitebg('w');
        ylabel('XZ');

    end
end

%%
dens = 0.03;
indx = omega>0.4 & z_cord>300 & z_cord<600;
figure();

titleString={'dx','dy','dz';};
data_for_plot = [dx,dy,dz];

count = 0;
for ii = 1:3
    
        count = count+1;
        subplot(1,3,ii);
        scatter(data_for_plot(indx,ii),SM_est_final(indx,5+5)-GT_final(indx,5+5),10,'filled','MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);
        xlabel(titleString{ii});
        set(gcf,'Color','w');
        whitebg('w');
        ylabel('XZ');

  
end