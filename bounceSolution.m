function bounceSolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unchanging parameter 
absTol = 1e-4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% changing parameter to adjust
r_range= 25;
r_start = 1e-5;
hr = 1e-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epsilon iteration
for epsilon = [0.67:0.01:1.0]
    dphi0 = 0.0;
    sigma = (3.0*(epsilon+1)+sqrt(9.0*(1+epsilon)^2-8))/4.0;
    sigma = 0.9999*sigma;
    r = [hr:hr:r_range];
    r = r';
    phis = zeros(length(r));
    phi = zeros(length(r));
    dphi = zeros(length(r));
    phi0_i = 0.0;
    phi0_t = sigma;
    isolution = false;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eom solver
    options = odeset('RelTol',absTol,'MaxStep',1e-2,'Events',@events);
    while(~isolution)
%         if(abs(phi0-(phi0_i+phi0_t)/2.0)<1e-10)
%             isolution = true;
%         end
        phi0 = (phi0_i+phi0_t)/2.0;
        dphi0 = -r_start^3;
        [r_raw,phis_raw,te,ye,ie] = ode15s(@bounceSolver,[r_start r_range],[phi0 dphi0],options); 
        if isempty(ie)
            isolution = true;
        else
        if phis_raw(end,1)<0
            phi0_t = phi0;
        else phi0_i = phi0;
        end     
        end
    end

    phis = interp1(r_raw,phis_raw,r);
    phi = phis(:,1);
    dphi = phis(:,2);
    filename = strcat('epsilon',num2str(epsilon,'%.2f'),'_phic',num2str(phi0,'%.4f'),'.mat');
    save(filename,'epsilon','r','phi','dphi');
    display('solution found');
    phi0/sigma
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% potential function
    function dV = dpotential(phi)
    dV = 0.5*(phi)-1.5*((phi).^2)*(1+epsilon)+((phi).^3);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution adjuster
    function [result,phi,phi_i,phi_t] = solutionAdjuster(ie,phi0,phi0_i,phi0_t)
        result = false;
        phi = phi0;
        phi_i = phi0_i;
        phi_t = phi0_t;
        if isempty(ie)
            result = true;
        elseif find(ie==1)
            phi_t = phi0;
        else phi_i = phi0;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution judger
    function [value,isterminal,direction] = events(r,phis)
        value = [phis(1) phis(2)];     % Detect height = 0
        isterminal = [0 0];   % Stop the integration
        direction = [0 0];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eom
function dy = bounceSolver(r,y)
dy= zeros(2,1);
dy(1) = y(2);
dy(2) = 0.5*(y(1))-1.5*((y(1)).^2)*(1+epsilon)+((y(1)).^3)-3.0*y(2)./r;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end