zSlice = linspace(000,1200,13);
count = 0;
for ii = 2:1:length(zSlice)-1
    count = count+1;
%Y{count} = [num2str(zSlice(ii)),'nm<z<',num2str(zSlice(ii+1)),'nm'];
Y{count} = [num2str(zSlice(ii)+50)];
end
X = categorical(Y);
X = reordercats(X,Y);
X = [150:100:1150];
figure();
plot(X,sigma_DPPC,'LineWidth',1.5,'Color',[0, 0.4470, 0.7410]); hold on;

%%

plot(X,sigma_DPPC,'LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250]); hold on;
plot(X,sigma_perfect,'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980]);
legend('DPPC+40%chol','DPPC','theory','EdgeColor','none','Color','none','FontSize',8); ylabel('FWHM(nm)');
set(gcf,'Color','w');
whitebg('w');
xlim([150,1150]);
xlabel('h(nm)');
%% -----------omega-theta_perpendicular scatter-----------

dens = 0.1;
%figure();
scatter(theta_perpendicular,omega(indx)/pi,2,'filled','MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);

hold on;
scatter(median(theta_perpendicular),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0, 102, 204]./255,'MarkerEdgeColor',[0, 102, 204]./255,'LineWidth',2);

ax = gca;
ax.FontSize =9; 
box on;
whitebg('w');
set(gcf,'Color','w');
ylabel('\Omega(\pi)'); xlabel('\theta\perp (\circ)');
legend('DPPC+40%chol','','DPPC','','FontSize',8,'EdgeColor','none','Color','none');
%% theta perpendicular - omega scatter plot 2
dens = 0.1;
hold on;
scatter(theta_perpendicular,omega(indx)/pi,2,'filled','MarkerFaceColor',	[0.9290, 0.6940, 0.1250],'MarkerFaceAlpha',dens,'MarkerEdgeAlpha',dens);
ylabel('\Omega(\pi)'); xlabel('\theta\perp  (\circ)');
hold on;
scatter(median(theta_perpendicular),median(omega(indx))/pi,100,'+','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'LineWidth',2);
ax = gca;
ax.FontSize = 9; 
box on;
whitebg('w');
set(gcf,'Color','w');
xticks([0,90])
xlim([0,90]);

