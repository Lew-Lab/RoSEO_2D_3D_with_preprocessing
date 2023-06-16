
load('20210504_DPPC_40chol_est_retrieval_1.1_xy_centered_v27.mat');
%
load('colorSpace.mat');
%load('est_xy_centered_v13.mat');


x_cord = x; y_cord = y;z_cord = z;  %-90

mux = sind(phiD); muy = cosd(phiD); 
phiD_cord = atan2(mux,-muy)/pi*180;

indx = abs(x_cord)<1200 & abs(y_cord)<1200 & z_cord<1200 & z_cord>000 &SM_est_final(:,5)>1000;

%plot 3D spikes
figure('Color',[0.1,0.1,0.1]);

hold on

load('turbo.mat');
%colormap(turbo);
%color = squeeze(colorSpace);
color = turbo;
thetaD_dic = linspace(0,90,256);
phiD_dic = linspace(-180,180,403);

count = 0;
r=50;
 for ii=1:length(indx)
    if indx(ii)==1
        count = count+1;
        if rem(count,4)==0
        [~,color_indx] = min(abs(thetaD(ii)-thetaD_dic));
        plot3([x_cord(ii)-r*sind(thetaD(ii)).*cosd(phiD_cord(ii)),x_cord(ii),x_cord(ii)+r*sind(thetaD(ii)).*cosd(phiD_cord(ii))].',...
             [y_cord(ii)-r*sind(thetaD(ii)).*sind(phiD_cord(ii)),y_cord(ii),y_cord(ii)+r*sind(thetaD(ii)).*sind(phiD_cord(ii))].',...
             [z_cord(ii)+r*cosd(thetaD(ii)),z_cord(ii),z_cord(ii)-r*cosd(thetaD(ii))].',...
             'Color',color(color_indx,:),'LineWidth',1);
        end
    end
end
axis image


xlabel('x(nm)'); ylabel('y(nm)'); zlabel('z(nm)');

ax = gca;
ax.FontSize = 10; 
whitebg('black');
set(gcf, 'InvertHardcopy', 'off');
phi_view = linspace(0,180,18);
theta_view = linspace(0,30,5);
%----------save image at different viewing angles---------------------
 for jj=1:length(theta_view)
     xPlane = [-1 1 1 -1]*1200;      % X coordinates of plane corners, ordered around the plane
    yPlane1 = [0,0,0,0];      % Corresponding y coordinates for plane 1
    zPlane = [0 0 1600 1600];  % Z coordinates of plane corners
    hold on;                   % Add to existing plot
    h = patch(xPlane, yPlane1, zPlane, 'k', 'FaceAlpha', 0.8);  % Plot plane 1
    rotate(h,[0,0,1],phi_view(jj));
    rotate(h,[1,0,0],360-theta_view(jj));

     view(phi_view(jj),theta_view(jj));
     set(gcf, 'Color', 'None')
     axis vis3d;
     saveas(gcf,['C:\Users\wutt0\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210504\15-24\spikes\',num2str(jj),'.png'])
     delete(h);
 end

%% z slice scattering plot
zSlice = linspace(-100,1200,14);

for ii=1:length(zSlice)-1
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) &SM_est_final(:,5)>400;
space = round(linspace(1,6380,4500));

figure('Color',[0.1,0.1,0.1]);
hold on

color = squeeze(colorSpace);
thetaD_dic = linspace(0,90,256);
phiD_dic = linspace(-180,180,403);
r = 40;
%caxis([-180,180]);
hold on;
%colorbar;
for jj=1:length(indx)
    if indx(jj)==1
        [~,color_indx] = min(abs(phiD_cord(jj)-phiD_dic));
        plot([x_cord(jj)-r*sind(thetaD(jj)).*cosd(phiD_cord(jj)),x_cord(jj),x_cord(jj)+r*sind(thetaD(jj)).*cosd(phiD_cord(jj))].',...
             [y_cord(jj)-r*sind(thetaD(jj)).*sind(phiD_cord(jj)),y_cord(jj),y_cord(jj)+r*sind(thetaD(jj)).*sind(phiD_cord(jj))].',...
             'Color',color(color_indx,:),'LineWidth',1);
    end
end
axis image
xlim([-1300,1300]);ylim([-1300,1300]);

%---
whitebg('k');
%colorbar('Color','w');
%xlabel('x','Color','w'); ylabel('y','Color','w'); 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
ax = gca;
ax.FontSize = 10; 
axis off;



set(gcf,'Color',[0 0 0]);
text(500,1100,[num2str(zSlice(ii)),'nm<z <',num2str(zSlice(ii+1)),'nm'],'Color','w','FontSize',10);
set(gcf, 'InvertHardcopy', 'off');
%view(phi_view(ii),theta_view(ii));
for kk=1:10
saveas(gcf,['C:\Users\wutt0\OneDrive - Washington University in St. Louis\github\data of PSF-optimization\3Dbeads_experiment\20210504\15-24\spikes3\',num2str(((ii-1)*10+kk)),'.png'])
end

end
