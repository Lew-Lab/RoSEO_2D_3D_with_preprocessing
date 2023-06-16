%%
addpath(genpath('###\3D SLB vidualization_demo\'))

load('colorSpace.mat');
colorinUse = squeeze(colorSpace);
%%
% 
% dataN = [1];
% fig_saveFolder = 'F:\OneDrive - Washington University in St. Louis\github\Deep-SMOLM3D\3D SLB vidualization\SLB data\DPPC_Figs\';
% readFunction = 'data_overview.DPPC_only_dataSource';


dataN = [4];
fig_saveFolder = '###\3D SLB vidualization_demo\SLB data\DPPC_chol_Figs\';
readFunction = 'data_overview.DPPC_chol_dataSource';

% dataN = [1];
% fig_saveFolder = 'F:\OneDrive - Washington University in St. Louis\github\Deep-SMOLM3D\3D SLB vidualization\SLB data\DPPC_chol_Figs\';
% readFunction = 'data_overview_network.DPPC_chol_dataSource';





generate_diffuse_check_frame = 0;
readFunction = str2func(readFunction);

for mm = dataN

[locFileName,imageName,angleEstName,center_SLB,radius_SLB,dataIDN,frameN_per_stack,bkgName, offsetName] = readFunction(mm);

cdir = dir(imageName);
fileFolder = cdir.folder;
locFolder = '\est_results\';


%% read loc
locList = [];
angleList = [];
R_save = [];
stackCount = 0;

for ii = dataIDN

    stackCount = stackCount+1;
    load([fileFolder,locFolder,num2str(ii), '_est_FoV1_v2.mat']);

    %***------------- add one index column***-------------*-------------
    temp = SM_est_save_all;
    temp(:,1) = Angle_save(:,1);
    temp(:,2:11) = SM_est_save_all;
    SM_est_save_all = temp;
    %*-------------*-------------*-------------*-------------


    %indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
    indx = abs(SM_est_save_all(:,2))<1700  & abs(SM_est_save_all(:,3))<1700 & SM_est_save_all(:,4)>400 & SM_est_save_all(:,4)<1100  & SM_est_save_all(:,5)>2000;
    [R,~]=sphereFit(SM_est_save_all(indx,2:4));
    R_save = [R_save;R]; 
    SM_est_save_all(:,2:3) = SM_est_save_all(:,2:3)-[R(1),R(2)];   %compensatet the lateral drift between different stacks

    locList_cur =  SM_est_save_all;
    angleList_cur = Angle_save;
    locList_cur(:,1) = locList_cur(:,1)+(stackCount-1)*frameN_per_stack;
    angleList_cur(:,1) = angleList_cur(:,1)+(stackCount-1)*frameN_per_stack;
    %****filping operator************ remove this part if no flipping is
    %need
    phiD = angleList_cur(:,3);
    mux = sind(phiD); muy = cosd(phiD); 
    phiD_cord = atan2(mux,-muy)/pi*180;
    angleList_cur(:,3) = phiD_cord;
    %***************
    locList = [locList;locList_cur];
    angleList = [angleList; angleList_cur];

    Fig1 = z_slice_scatterPlot(locList_cur, angleList_cur, colorinUse);
    exportgraphics(Fig1, [fig_saveFolder, 'data', num2str(mm), '_stack_' , num2str(ii),'_z_slices_image.jpg'],'Resolution',500,'BackgroundColor','k');
close all
end
%%
    Fig1 = z_slice_scatterPlot(locList, angleList, colorinUse);
    exportgraphics(Fig1, [fig_saveFolder, 'data', num2str(mm), '_stack_all_z_slices_image.jpg'],'Resolution',500,'BackgroundColor','k');
    %%
    Fig2 = theta_scatter(locList, angleList);
    exportgraphics(Fig2, [fig_saveFolder, 'data', num2str(mm), '_stack_all_theta_image.jpg'],'Resolution',500,'BackgroundColor','k');

    Fig3 = threeD_scatter(locList, angleList,colorinUse);
    exportgraphics(Fig3, [fig_saveFolder, 'data', num2str(mm), '_stack_all_3D_image.jpg'],'Resolution',500,'BackgroundColor','k');
%%
save(locFileName,'locList', 'angleList','center_SLB','radius_SLB','R_save');

end


%%

function Fig2 = theta_scatter(loc,ang)
Fig2 = figure('Color',[0.1,0.1,0.1]); hold on
set(Fig2, 'Units','centimeters','InnerPosition', [10 10 16 12]); sgtitle('\theta estimations');

x_cord = loc(:,2);
y_cord = loc(:,3);
z_cord = loc(:,4);
signal = loc(:,5);
thetaD = ang(:,2);
phiD = ang(:,3);
indx = abs(x_cord)<2000 & abs(y_cord)<2000 & z_cord<1000 & z_cord>00 & signal>1000;

subplot(1,2,1);
isColorbar=1; isScaleBar=1;
plot_scatter(x_cord(indx),z_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);

subplot(1,2,2);
plot_scatter(x_cord(indx),y_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);

end

function plot_scatter(x_cord,y_cord,colorPara, colormapC,isColorbar,isScaleBar)
scatter(x_cord,y_cord,3,colorPara,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); axis image; 
colormap(colormapC);
whitebg('k'); %#ok<WHITEBG>
ax = gca;
ax.XColor = 'k'; ax.YColor = 'k';
set(gca,'Color','k')
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

function Fig = z_slice_scatterPlot(loc,ang, colorUse)
x_cord = loc(:,2);
y_cord = loc(:,3);
z_cord = loc(:,4);
signal = loc(:,5);
thetaD = ang(:,2);
phiD = ang(:,3);
% orientation in each z slice

Fig = figure('Color',[0.1,0.1,0.1]); hold on
set(Fig, 'Units','centimeters','InnerPosition', [10 10 35 12]*1.0);
set(gcf, 'InvertHardcopy', 'off');
%-----------------------------plotting parameters-----------------------------
r = 70; %line length
zSlice = linspace(000,1300,13); d_slice = zSlice(2)-zSlice(1); % plotted slice 
x_offset = [0,sqrt(1000^2-(1000-(zSlice(2:end))))*1.2]; %the distance between each circle plot
phiD_dic = linspace(-180,180,403); %partition the colormap with phi
%-----------------------------

for ii=1:2:length(zSlice)-1
indx = abs(x_cord)<1300 & abs(y_cord)<1300 & z_cord<zSlice(ii+1) & z_cord>zSlice(ii) & signal>1000;

for jj=1:length(indx)
    if indx(jj)==1
        [~,color_indx] = min(abs(phiD(jj)-phiD_dic));
        plot(sum(x_offset(1:ii+1))+[x_cord(jj)-r*sind(thetaD(jj)).*cosd(phiD(jj)),x_cord(jj),x_cord(jj)+r*sind(thetaD(jj)).*cosd(phiD(jj))].',...
             [y_cord(jj)-r*sind(thetaD(jj)).*sind(phiD(jj)),y_cord(jj),y_cord(jj)+r*sind(thetaD(jj)).*sind(phiD(jj))].',...
             'Color',colorUse(color_indx,:),'LineWidth',0.5);
    end
end
axis image
%xlim([-1300,1300]);ylim([-1300,1300]);

%-------- plot reference sphere-----------------------------
R_cur = sqrt(1000^2-(1000- zSlice(ii)-d_slice/2).^2);
R_top = sqrt(1000^2-(1000- zSlice(ii)-d_slice).^2);
R_bottom = sqrt(1000^2-(1000- zSlice(ii)).^2);
if zSlice(ii)+d_slice<=0
    if zSlice(ii)+d_slice/2<=0
        R_cur=0;
    end
    R_bottom=0;
end

%-------- figure setting-----------------------------
set(gca,'xtick',[]);
set(gca,'ytick',[]);
ax = gca;
ax.FontSize = 10; 
axis off;
set(gcf, 'InvertHardcopy', 'off');
whitebg('k');
set(gcf,'Color',[0 0 0]);

%text(-700+sum(x_offset(1:ii+1)),1300,[num2str(round(zSlice(ii))),'nm<h <',num2str(round(zSlice(ii+1))),'nm'],'Color','w','FontSize',8);

text(-700+sum(x_offset(1:ii+1)),1300,[num2str(round(zSlice(ii))),'nm<h <',num2str(round(zSlice(ii+1))),'nm'],'Color','w','FontSize',8);

if ii==12
    plot([sum(x_offset(1:ii+1))+600,sum(x_offset(1:ii+1))+600+400],[-1100,-1100],'Color','w','LineWidth',2);
end
end

end

function Fig2 = threeD_scatter(loc,ang,colorinUse)
Fig2 = figure('Color',[0.1,0.1,0.1]); hold on
set(Fig2, 'Units','centimeters','InnerPosition', [10 10 16 12]); sgtitle('\theta estimations');

x_cord = loc(:,2);
y_cord = loc(:,3);
z_cord = loc(:,4);
signal = loc(:,5);
thetaD = ang(:,2);
phiD = ang(:,3);
gamma = ang(:,5);
%indx = abs(x_cord)<1200 & abs(y_cord)<1200 & z_cord<1000 & z_cord>00 & signal>1000;
indx = signal>500 & z_cord<1300 & z_cord>00;

subplot(1,2,1);
isColorbar=0; isScaleBar=0;
%plot_scatter(x_cord(indx),z_cord(indx),thetaD(indx), turbo(256),isColorbar,isScaleBar);



scatter3(x_cord(indx),y_cord(indx),z_cord(indx),3,gamma(indx) ,'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5); axis image; 
%colormap(colorinUse);
whitebg('k'); %#ok<WHITEBG>
ax = gca;
ax.XColor = 'k'; ax.YColor = 'k';
set(gca,'Color','k')
ax.FontSize = 10; 
set(gca,'xtick',[]);
set(gca,'ytick',[]);

end