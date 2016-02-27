figure;
phi = [-0.25:0.01:1.4];
epsilon = 0.02;
V = -epsilon*(phi.^3)/2.0+ (phi.^2).*((phi-1).^2)/4.0;
plot(phi,V,'LineWidth',4);

% SEdensity = (dphi.^2)/2.0-epsilon*(phi.^3)/2.0+ (phi.^2).*((phi-1).^2)/4.0;
% plot(r,SEdensity,'LineWidth',4);

set(gcf,'Position',[0,0,500,500],'Color','w');
set(gca,'LineWidth',2,'FontSize',16);
ax = xlabel('$\phi$');
yx = ylabel('$V$');
set(ax,'Interpreter','latex','FontWeight','Bold','FontSize',24);
set(yx,'Interpreter','latex','FontWeight','Bold','FontSize',24);
hold off;