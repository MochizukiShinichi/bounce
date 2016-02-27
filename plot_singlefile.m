V = -epsilon*(phi.^3)/2.0+ (phi.^2).*((phi-1).^2)/4.0;
KE = (dphi.^2)/2.0;
SEdensity = (dphi.^2)/2.0-epsilon*(phi.^3)/2.0+ (phi.^2).*((phi-1).^2)/4.0;
figure;
plot(r,V,'LineWidth',4);
hold on;
plot(r,KE,'LineWidth',4);
plot(r,SEdensity,'LineWidth',4);
plot(r,KE-V,'LineWidth',4);
legend('V','KE','SE','\rho');