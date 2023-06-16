

save('simulated_pixOL_fixed_0624data.mat','x_RMSE','y_RMSE','z_RMSE','zSlice');
%% scatter plot of the 3D estimation

load('colorSpace.mat');
%load('est_xy_centered_v13.mat');

z_offset =0;
x_cord = x-56.5; y_cord = y-56.5;z_cord = z;
phiD_cord = phiD;



indx = abs(x_cord)<1200 & abs(y_cord)<1200 & z_cord<1200 & z_cord>000 &SM_est_final(:,5)>400;

[R,C]=sphereFit([x_cord(indx), y_cord(indx),z_cord(indx)]);


figure('Color',[0.1,0.1,0.1]);
hold on
scatter3(x_cord(indx),y_cord(indx),z_cord(indx),20,phiD_cord(indx),'filled','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',.7); axis image; 
xlabel('x(nm)'); ylabel('y(nm)'); zlabel('z(nm)');
%colormap(turbo);
color = squeeze(colorSpace);
colormap(color);


ax = gca;
ax.FontSize = 10; 
whitebg('black');
set(gcf, 'InvertHardcopy', 'off');


%% orientation in each z slice

Fig = figure('Color',[0.1,0.1,0.1]); hold on
set(Fig, 'Units','centimeters','InnerPosition', [10 10 35 12]*1.0);

%-----------------------------plotting parameters-----------------------------
r = 70; %line length
zSlice = linspace(000,800,9); d_slice = zSlice(2)-zSlice(1); % plotted slice 
x_offset = [0,sqrt(1000^2-(1000-(zSlice(2:end))))*2.2]; %the distance between each circle plot
load('colorSpace.mat'); color = squeeze(colorSpace); % colormap
phiD_dic = linspace(-180,180,403); %partition the colormap with phi
%-----------------------------

for ii=1:1:length(zSlice)-1
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) &SM_est_final(:,5)>000;

for jj=1:length(indx)
    if indx(jj)==1
        [~,color_indx] = min(abs(phiD_cord(jj)-phiD_dic));
        plot(sum(x_offset(1:ii+1))+[x_cord(jj)-r*sind(thetaD(jj)).*cosd(phiD_cord(jj)),x_cord(jj),x_cord(jj)+r*sind(thetaD(jj)).*cosd(phiD_cord(jj))].',...
             [y_cord(jj)-r*sind(thetaD(jj)).*sind(phiD_cord(jj)),y_cord(jj),y_cord(jj)+r*sind(thetaD(jj)).*sind(phiD_cord(jj))].',...
             'Color',color(color_indx,:),'LineWidth',0.5);
    end
end
axis image
%xlim([-1300,1300]);ylim([-1300,1300]);

%-------- plot reference sphere-----------------------------
% R_cur = sqrt(1000^2-(1000- zSlice(ii)-d_slice/2).^2);
% R_top = sqrt(1000^2-(1000- zSlice(ii)-d_slice).^2);
% R_bottom = sqrt(1000^2-(1000- zSlice(ii)).^2);
% if zSlice(ii)+d_slice<=0
%     if zSlice(ii)+d_slice/2<=0
%         R_cur=0;
%     end
%     R_bottom=0;
% end

%-------- figure setting-----------------------------
set(gca,'xtick',[]);
set(gca,'ytick',[]);
ax = gca;
ax.FontSize = 10; 
axis off;
set(gcf, 'InvertHardcopy', 'off');
whitebg('k');
set(gcf,'Color',[0 0 0]);

text(-700+sum(x_offset(1:ii+1)),1300,[num2str(zSlice(ii)),'nm<h <',num2str(zSlice(ii+1)),'nm'],'Color','w','FontSize',8);


if ii==12
    plot([sum(x_offset(1:ii+1))+600,sum(x_offset(1:ii+1))+600+400],[-1100,-1100],'Color','w','LineWidth',2);
end
end

%% lateral precision and axial precision by calculating the RMSE

load('pixOL_loc_precision_NFP-700.mat');
sigma_r = sqrt(sigma_x_opt.^2+sigma_y_opt.^2);

zSlice = linspace(000,900,10); d_slice = zSlice(2)-zSlice(1); % plotted slice 
x_GT = GT_final(:,1); y_GT = GT_final(:,2); z_GT = GT_final(:,3);
x_est = SM_est_final(:,2); y_est = SM_est_final(:,3); z_est = SM_est_final(:,4);
for ii = 1:1:length(zSlice)-1
    indx = z_GT<zSlice(ii+1) & z_GT>zSlice(ii);
    
    x_RMSE(ii) = sqrt(mean((x_est(indx)-56.5-x_GT(indx)).^2));
    y_RMSE(ii) = sqrt(mean((y_est(indx)-56.5-y_GT(indx)).^2));
    z_RMSE(ii) = sqrt(mean((z_est(indx)-z_GT(indx)).^2));
    
    r_RMSE(ii) = sqrt(mean((y_est(indx)-y_GT(indx)).^2+(x_est(indx)-x_GT(indx)).^2));
end



figure(); hold on;
plot(zSlice(1:end-1),x_RMSE);
plot(zSlice(1:end-1),y_RMSE);
plot(zSlice(1:end-1),z_RMSE);
%plot(zSlice(1:end-1),r_RMSE);
legend('x','y','z');
whitebg('w');
%% lateral location precision quantification
% reference: perfect pixOL precision
load('pixOL_loc_precision_NFP-700.mat');
sigma_r = sqrt(sigma_x_opt.^2+sigma_y_opt.^2);

zSlice = linspace(000,900,10); d_slice = zSlice(2)-zSlice(1); % plotted slice 

figure();
count = 0;
for ii=1:1:length(zSlice)-1
    count = count+1;
    
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) &SM_est_final(:,5)>400;
r_data = sqrt((x_cord(indx).^2+y_cord(indx).^2));

%------------  reference sphere  ----------------
R_cur = sqrt(1000^2-(1000- zSlice(ii)-d_slice/2).^2);
R_top = sqrt(1000^2-(1000- zSlice(ii)-d_slice).^2);
R_bottom = sqrt(1000^2-(1000- zSlice(ii)).^2);
if zSlice(ii)+d_slice<=0
    if zSlice(ii)+d_slice/2<=0
        R_cur=0;
    end
    R_bottom=0;
end

if zSlice(ii)>=1000
   temp1 = max([R_top,R_bottom,R_cur]);
   temp2 = min([R_top,R_bottom,R_cur]);
   R_top = temp1;
   R_bottom = temp2;
end
%------------------------------
[~,indx_top]=min(abs(zSlice(ii+1)-z_range*10^9));
sigma_top = sigma_r(indx_top)*10^9;

[~,indx_bottom]=min(abs(zSlice(ii)-z_range*10^9));
sigma_bottom = sigma_r(indx_bottom)*10^9;

%subplot(1,6,count); hold on;

N_bins = length(r_data)/40;
edges=linspace((prctile(r_data,2)-50),(prctile(r_data,98)+50),N_bins);
edges_distance = edges(2)-edges(1);
h = histogram(r_data,edges);
h.FaceColor = [128, 255, 0]/255;
%h.FaceColor = [0.4660, 0.6740, 0.1880];
h.EdgeColor = 'none';

%-------------calculate the FWHM -----------------

pdf_data = h.Values; center_data = h.BinEdges;
[max_Data,indx_max] = max(pdf_data);
half_max= max_Data/2;
[~,indx_left] = min(abs(pdf_data(1:indx_max)-half_max));
[~,indx_right] = min(abs(pdf_data(indx_max+1:end)-half_max));
data_FWHM(ii) = center_data(indx_right+indx_max)-center_data(indx_left);


x_ax = R_bottom-180:R_top+180;
if ii==12
   x_ax =  R_bottom-270:R_top+270;
end
gauss_pdf = normpdf((1:length(x_ax))-length(x_ax)/2,0,(sigma_top+sigma_bottom)/2);
circ_pdf = zeros(size(x_ax));
circ_pdf = (x_ax+1)./cos(asin((x_ax+1)/1000))-(x_ax)./cos(asin((x_ax)/1000));
circ_pdf(x_ax<(R_bottom)|x_ax>(R_top)) = 0; 
circ_pdf = real(circ_pdf);
circ_pdf(circ_pdf<0) = 0;
if sum(circ_pdf)==0
   pdf = gauss_pdf; pdf = pdf/sum(pdf);
else
circ_pdf = circ_pdf/sum(circ_pdf);
pdf = conv(circ_pdf,gauss_pdf,'same');
pdf = pdf/sum(pdf);
end
pdf = pdf*max(pdf_data)/max(pdf);
plot(x_ax,pdf,'Color',[255 255 153]/255, 'LineWidth',1.5);
%plot(x_ax,pdf*sum(indx)*edges_distance,'Color',[0.9290, 0.6940, 0.1250], 'LineWidth',1.5);
xlabel('r (nm)'); ylabel('count');
ax = gca;
ax.FontSize = 8;
xlim([min(x_ax(1),edges(1)),max(x_ax(end),edges(end))]);

%-------------calculate the FWHM -----------------

%
%loc_data_distr = deconv(pdf_data/sum(pdf_data),gauss_pdf);
pdf_temp = circ_pdf*sum(indx)*edges_distance;
[max_Data,indx_max] = max(pdf_temp);
half_max= max_Data/2;
[~,indx_left] = min(abs(pdf_temp(1:indx_max)-half_max));
[~,indx_right] = min(abs(pdf_temp(indx_max+1:end)-half_max));
pfd_FWHM_sphere(ii) =indx_right+indx_max-indx_left;

%
pdf_temp = pdf;
[max_Data,indx_max] = max(pdf_temp);
half_max= max_Data/2;
[~,indx_left] = min(abs(pdf_temp(1:indx_max)-half_max));
[~,indx_right] = min(abs(pdf_temp(indx_max+1:end)-half_max));
pfd_FWHM_whole(ii) = indx_right+indx_max-indx_left;



set(gcf,'Color','k');
whitebg('k');

end

%% plot the FWHM
pfd_FWHM_sphere(end)=0;
count = 0;
for ii = 2:1:length(zSlice)-1
    count = count+1;
Y{count} = [num2str(zSlice(ii)),'<z <',num2str(zSlice(ii+1))];
end
X = categorical(Y);
X = reordercats(X,Y);

sigma_DPPC = (data_FWHM(2:1:end));
sigma_DPPC(sigma_DPPC<0)=0;
sigma_perfect = (pfd_FWHM_whole(2:1:end));
figure();
plot(X,sigma_DPPC); hold on;
plot(X,sigma_perfect);
legend('DPPC+40%Chol data','perfect pixOL PSF','EdgeColor','none'); ylabel('FWHM(nm)');
set(gcf,'Color','w');
whitebg('w');

%% theta
isColorbar=1; isScaleBar=1;
indx = abs(x_cord)<1200 & abs(y_cord)<1200 & z_cord<1000 & z_cord>00 &SM_est_final(:,5)>1000;
plot_scatter(x_cord(indx),y_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);
plot_scatter(x_cord(indx),z_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);
plot_scatter(y_cord(indx),z_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);


%% phi
isColorbar=0;
plot_scatter(x_cord(indx),y_cord(indx),phiD_cord(indx), squeeze(colorSpace),isColorbar,isScaleBar);
plot_scatter(x_cord(indx),z_cord(indx),phiD_cord(indx), squeeze(colorSpace),isColorbar,isScaleBar);
plot_scatter(y_cord(indx),z_cord(indx),phiD_cord(indx), squeeze(colorSpace),isColorbar,isScaleBar);


%% omega
isColorbar=0;
plot_scatter(x_cord(indx),y_cord(indx),omega(indx)/pi,parula(256),isColorbar,isScaleBar);
plot_scatter(x_cord(indx),z_cord(indx),omega(indx)/pi,parula(256),isColorbar,isScaleBar);
plot_scatter(y_cord(indx),z_cord(indx),omega(indx)/pi,parula(256),isColorbar,isScaleBar);

%% theta perpendicular histogram
[X_ref,Y_ref,Z_ref,thetaD_ref,phiD_ref]=sphereLine(x_cord,y_cord,z_cord,1000);

indx(:)=1;
%
figure();
vector_est = [sind(thetaD(indx)).*cosd(phiD_cord(indx)),sind(thetaD(indx)).*sind(phiD_cord(indx)),cosd(thetaD(indx))];
vector_per = real([sind(thetaD_ref(indx)).*cosd(phiD_ref(indx)),sind(thetaD_ref(indx)).*sind(phiD_ref(indx)),cosd(thetaD_ref(indx))]);
theta_perpendicular = real(acos(sum((vector_est.*vector_per),2))./pi*180);
theta_perpendicular(theta_perpendicular>90)=180-theta_perpendicular(theta_perpendicular>90);
histogram(theta_perpendicular,50); %title('normal \theta');
xlabel('normal \theta(\circ)','Color','k'); 
ax = gca;
ax.FontSize = 10; 
ylabel('count');
set(gcf,'Color','w');
whitebg('w');
%save('data_for_DPPC_v13.mat','indx','x_cord','y_cord','z_cord','thetaD','phiD_cord','omega','sigma_DPPC','sigma_perfect','colorSpace','thetaD_ref','phiD_ref','theta_perpendicular','z_offset');
%% Omega histogram
figure();
histogram(omega(indx)/pi,30); 
xlim([0,2]);
ax = gca;
ax.FontSize = 10; 
xlabel('\Omega (\pi)'); 
%title('\Omega (\pi)');
ylabel('count');
set(gcf,'Color','w');
whitebg('w');
%% theta perpendicular - omega scatter plot 1
dens = 0.05;
figure();
scatter(theta_perpendicular(indx),omega(indx)/pi,8,'filled','MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);

hold on;
scatter(median(theta_perpendicular(indx)),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0, 128, 255]./255,'MarkerEdgeColor',[0, 128, 255]./255,'LineWidth',1);

ax = gca;
ax.FontSize = 10; 
box on;
whitebg('w');
set(gcf,'Color','w');
ylabel('\Omega(\pi)'); xlabel(' \theta perpendicular (\circ)');
%% theta perpendicular - omega scatter plot 2
dens = 0.03;
hold on;
scatter(theta_perpendicular(indx),omega(indx)/pi,8,'filled','MarkerFaceColor',	[0.9290, 0.6940, 0.1250],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);
ylabel('\Omega(\pi)'); xlabel(' \theta perpendicular (\circ)');
hold on;
scatter(median(theta_perpendicular(indx)),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'LineWidth',1);
ax = gca;
ax.FontSize = 10; 
box on;
whitebg('w');
set(gcf,'Color','w');

%% bind the omega theta perpendicular data
theta_bin_edges = linspace(0,90,40);
omega_bin_edges = linspace(0,2,40)*pi;
hotN = bin2Ddata(omega(indx),theta_perpendicular,omega_bin_edges,theta_bin_edges);
figure();
imagesc(theta_bin_edges,omega_bin_edges/pi,hotN);   colormap('hot');
ax = gca;
ax.YDir = 'normal'; %caxis([0,70]);

%%  bind the phi - theta perpendicular data
phi_bin_edges = linspace(-180,180,40);
omega_bin_edges = linspace(0,2,40)*pi;
hotN = bin2Ddata(omega(indx),phiD_cord(indx),omega_bin_edges,phi_bin_edges);
figure();
imagesc(phi_bin_edges,omega_bin_edges/pi,hotN);  colormap('hot');  
ax = gca;
ax.YDir = 'normal';

%% photon distributopn
figure();
histogram(SM_est_final(:,5));
xlabel('photons','FontSize',10); ylabel('counts','FontSize',10);
title({['mean:' num2str(mean(SM_est_final(:,5)))],['median:' num2str(median(SM_est_final(:,5)))]});
ax = gca;
ax.FontSize = 10; 




%% ---------- function --------------
function plot_scatter(x_cord,y_cord,colorPara, colormapC,isColorbar,isScaleBar)
figure()
scatter(x_cord,y_cord,3,colorPara,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); axis image; 
colormap(colormapC);
whitebg('k'); %#ok<WHITEBG>
ax = gca;
ax.XColor = 'k'; ax.YColor = 'k';
ax.FontSize = 10; 
set(gca,'xtick',[]);
set(gca,'ytick',[]);
if isColorbar==1
    hcb = colorbar('Color','k');
    %hcb.Title.String = "\theta(\circ)";
end
if isScaleBar == 1
   hold on;
   plot([max(x_cord)-500,max(x_cord)-500+400],[min(y_cord)+150,min(y_cord)+150],'Color','w','LineWidth',2);
    
end
%viscircles([0,0],1000)
%xlabel('x','Color','k'); ylabel('y','Color','k'); 

end