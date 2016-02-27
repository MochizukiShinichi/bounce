%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test region

    % CE = CE./(Phic.^4);
    SE = (exp(-SE).*(Radius.^4));
    CE =CE/norm(CE);
    SE = SE/norm(SE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% action region
    % CE = (CE-mean(CE))/std(CE);
    % SE = (SE-mean(SE))/std(SE);
    % CE =CE/norm(CE);
    % SE = SE/norm(SE);
    figure;
    plot(Epsilon,CE,'LineWidth',4);
    hold on;
    plot(Epsilon,SE,'LineWidth',4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot format
    set(gcf,'Position',[0,0,500,500],'Color','w');
    set(gca,'LineWidth',2,'FontSize',16);
    ax = xlabel('$\epsilon$');
    set(ax,'Interpreter','latex','FontWeight','Bold','FontSize',24);
%     set(yx,'Interpreter','latex','FontWeight','Bold','FontSize',24);
    lh = legend('CE','Tunnelling Amplitude' );
    set(lh,'Interpreter','latex','FontSize',16);
    hold off;