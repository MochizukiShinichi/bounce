%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import files
dos('ls epsilon*_phic*.mat > epiphic.csv');
fn=importdata('epiphic.csv'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize storage variables
Lf = length(fn);
CE = zeros(Lf,1);
SE = zeros(Lf,1);
Epsilon = zeros(Lf,1);
Phic = zeros(Lf,1);
Radius = zeros(Lf,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Lf
    load(char(fn(i)));
    Epsilon(i) = epsilon;
    Phic(i) = phi(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Action
    rho = (dphi.^2)/2.0+epsilon*(phi.^3)/2.0- (phi.^2).*((phi-1).^2)/4.0;
    SEdensity = (dphi.^2)/2.0-epsilon*(phi.^3)/2.0+ (phi.^2).*((phi-1).^2)/4.0;
    SE(i) = dot(SEdensity,r.^3);
    Radius(i) = sum((r.^4).*rho)/sum(r.^3.*rho);
    SE(i) = 2*(pi^2)*r(1)*SE(i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CE 
    hk = 1e-4;
    k_end = 100.0/r(end);
    k=linspace(hk,k_end,2.0^10);
    k = k';
    hk = k(2)-k(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4D version
    fft = ((k.^(-1))*(r'.^2).*besselj(1,k*r'))*phi;
    fk=fft.*conj(fft);
    fkmax = max(fk);
    fk = fk/fkmax;
    CEdensity = (fk.*log(fk)).*(k.^(3));
    CE(i)=-2*pi*pi*hk*sum(CEdensity);
    filename = strcat('epsilon',num2str(epsilon,'%.2f'),'kfile','.mat');
%     save(filename,'epsilon','k','fft','fk','CEdensity');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3D version
%     fft = ((k.^(-0.5))*(r'.^1.5).*besselj(0.5,k*r'))*phi;
%     fk=fft.*conj(fft);
%     fkmax = max(fk);
%     fk = fk/fkmax;
%     CEdensity = (fk.*log(fk)).*(k.^2);
%     CE(i)=-4*pi*hk*sum(CEdensity);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result data process
% CE = (CE-mean(CE))/std(CE);
% SE = (SE-mean(SE))/std(SE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('result_4d.mat','Epsilon','Phic','CE','SE','Radius');

