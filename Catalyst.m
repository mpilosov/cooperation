% This program solves the PDE's used in the 2014 paper entitled
% "Group-Level Events Are Catalysts in the Evolution of Cooperation" 
% Authors: Burt Simon & Michael Pilosov
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

clear all; format shortG; tic;

Xmax = 75; Ymax = 75; n = 75;       % grid parameters. default: n = 75
dx = Xmax/n; dy = Ymax/n;           % step sizes computed

% TIME CONTROL PARAMETERS
T = 8000;                           % time (years)
dt = .04;                           % temporal step size
freq = 2;                           % write to Theta every 'freq' years (positive integer). If set to 'dt', every frame will be written. 
prog = 10;                          % every 'prog' years, data will print to screen.

 
% KIN SELECTION
r = 0.15;            % base case relatedness. Hamilton's rule: rb>c => cooperation 
b = 0.4; c = 0.1;     % benefit & cost parameters in Prisoner's Dilemma 

% BIRTH/DEATH RATE
d = 0.0005;          % Death rate parameter
s = 0.025;           % birth/death rate scaling factor for within-group dynamics

% GROUP SELECTION
rho = 1;            % base case group-level events (benchmark rho = 1)
e0 = 0.025; phi = 1; % extinction rate e(x,y) = rho*e0 * Groups / (-phi*x+y)^2 
lambda = 0.00025;    % fissioning rate f(x,y) = rho*lambda*(x+y)
psi = 0;            % fission density = (1-psi)*h_1 + psi*h_2 
mu = 0.0025;         % per individual migration rate (avg rate a person switches groups per time unit dt) 
init0 = 471.25;
init_v = 5;
init_D=48;
init_C=2;

% OTHER PARAMETERS
% mutC=0; mutD=0;   % mutation rates for coops/defecs (future work)

% current baseline:
% b=0.4 | c=0.1 | s=0.025 | d=0.0005 | e0=0.025 | phi=1 | psi=0 | lambda=0.00025
% mu=0.00025 | rho=1 | r=0.15 | init_C=2 | init_D=48 | init_v=5 | init0=471.25




%-----------------------STRUCTURAL EDITING CAN OCCUR BELOW THIS LINE---------------------

% clearvars -except n dx dy dt freq prog T init_C init_D init_v init0 ...
%     r b c d s rho e0 phi lambda psi mutC mutD mu tim Xmax Ymax % use if nesting following section in a for-loop

alphaC = zeros(n); alphaD = zeros(n);                       % initialize vectors/matrices
growthC = zeros(n); growthD = zeros(n);
theta = zeros(n);
Theta = zeros(n,n,T/freq);
DeltaB = zeros(n,n,T/freq);
Population = zeros(T/freq,5);             % [t*dt Groups Coops Defectors Relatedness]

for i = 1:n                                                 % initialize theta and find
    for j = 1:n                                             % baseline Coop/Def growth rates
        x = (i-.5)*dx; y = (j-.5)*dy; pC = x/(x+y); pD = 1-pC;
        betaC = r*(1+b-c) + (1-r)*(pC*(1+b-c) + pD*(1-c));	% expected cooperator payoff
        betaD = r + (1-r)*(pC*(1+b) + pD);                  % expected defector payoff
        BrateC = x*betaC;                                   % total birth rates
        BrateD = y*betaD;
        %BrateC = x*betaC*(1-mutC) + y*betaD*mutD; % For mutation
        %BrateD = x*betaC*mutC + y*betaD*(1-mutD);
        growthC(i,j) = s * BrateC - x*d*(x+y);              % scaled net growth rates for (x,y)-group
        growthD(i,j) = s * BrateD - y*d*(x+y);              % alpha(i,j)in main loop since it depends on C,D,&G
        theta(i,j) = init0/(2*pi*init_v) * exp( -( (x-init_C)^2 + (y-init_D)^2 ) / (2*init_v));
    end
end
theta(:,n) = 0; theta(n,:) = 0; % set boundary to 0

fprintf('\n r=%.2f | rho=%.2f | mu=%.4f | psi=%.2f | phi=%.2f | init(C,D,v,0)=(%d,%d,%d,%.2f) | (dx,dy,dt)=(%.0f,%.0f,%.2f) \n\n ', ...
    r, rho, mu, psi, phi, init_C, init_D, init_v, init0, dx, dy, dt );
fprintf('      T    Groups    Delta\t  Coops    Delta \t   Defecs    Delta       %%Coops \n')
G1 = 0; C1 = 0; D1 = 0;
for t = 1:ceil(T/freq)*freq/dt  % MAIN LOOP
    if t*dt/prog == floor(t*dt/prog)
        if t*dt/prog==1
            fprintf('\n %8.2f  %8.2f  %+6.2f \t %8.2f  %+6.2f\t%8.2f  %+6.2f   %8.2f\n',...
                t*dt, Groups, Groups-G1, Coops, Coops-C1, Defectors, Defectors-D1, 100*Coops/(Coops+Defectors));  % current stats to screen
        else
            fprintf('\n %8.2f  %8.2f  %+6.2f \t %8.2f  %+6.4f\t%8.2f  %+6.4f   %8.2f\n',...
                t*dt, Groups, Groups-G1, Coops, Coops-C1, Defectors, Defectors-D1, 100*Coops/(Coops+Defectors));  % current stats to screen
        end
        G1 = Groups;  C1 = Coops;  D1 = Defectors;
    end
    Groups = 0; Coops = 0; Defectors = 0; R1 = 0; R2 = 0;
    for i = 1:n                 % integrate to get Groups, Coops, Defecs, Relatedness
        for j = 1:n
            x = (i-.5)*dx; y = (j-.5)*dy;
            Groups = Groups + theta(i,j)*dy*dx;
            Coops = Coops + x*theta(i,j)*dy*dx;
            Defectors = Defectors + y*theta(i,j)*dy*dx;
            R1 = R1 + (x^2/(x+y))* theta(i,j)*dy*dx;
            R2 = R2 + (x*y/(x+y))* theta(i,j)*dy*dx;
        end
    end
    Relate = r + (1-r)*(R1/Coops - R2/Defectors);
    
    for i = 1:n                 % compute alpha's from current state
        for j = 1:n
            x = (i-.5)*dx; y = (j-.5)*dy;
            alphaC(i,j) = growthC(i,j) + mu*(Coops/Groups - x);
            alphaD(i,j) = growthD(i,j) + mu*(Defectors/Groups - y);
        end
    end
    
    % PARTIAL DERIVATIVE TERMS FROM INDIVIDUAL DYNAMICS
    newtheta = zeros(n);
    for i = 1:n-1
        for j = 1:n-1
            A = (dy - abs(alphaD(i,j))*dt) * abs(alphaC(i,j))*dt / (dy*dx);
            B = abs(alphaD(i,j))*dt * abs(alphaC(i,j))*dt / (dy*dx);
            C = (dx - abs(alphaC(i,j))*dt) * abs(alphaD(i,j))*dt / (dy*dx);
            if A<0 || C<0       % error alert: dt too big
                disp([A C t i j]);
                error('Error: dt too large');
            end
            D = 1 - A - B - C;
            
            if alphaD(i,j)<0 && alphaC(i,j)<0
                if i==1         % enforcing boundary conditions
                    newtheta(i,j) = newtheta(i,j) + A*theta(i,j);
                    if j==1
                        newtheta(i,j) = newtheta(i,j) + B*theta(i,j);
                    else
                        newtheta(i,j-1) = newtheta(i,j-1) + B*theta(i,j);
                    end
                else
                    newtheta(i-1,j) = newtheta(i-1,j) + A*theta(i,j);
                    if j==1
                        newtheta(i-1,j) = newtheta(i-1,j) + B*theta(i,j);
                    else
                        newtheta(i-1,j-1) = newtheta(i-1,j-1) + B*theta(i,j);
                    end
                end
                if j==1
                    newtheta(i,j) = newtheta(i,j) + C*theta(i,j);
                else
                    newtheta(i,j-1) = newtheta(i,j-1) + C*theta(i,j);
                end
                newtheta(i,j) = newtheta(i,j) + D*theta(i,j);
            elseif alphaD(i,j)<0 && alphaC(i,j)>=0
                newtheta(i+1,j) = newtheta(i+1,j) + A*theta(i,j);
                if j==1         % enforcing boundary conditions
                    newtheta(i+1,j) = newtheta(i+1,j) + B*theta(i,j);
                    newtheta(i,j) = newtheta(i,j) + C*theta(i,j);
                else
                    newtheta(i+1,j-1) = newtheta(i+1,j-1) + B*theta(i,j);
                    newtheta(i,j-1) = newtheta(i,j-1) + C*theta(i,j);
                end
                newtheta(i,j) = newtheta(i,j) + D*theta(i,j);
            elseif alphaD(i,j)>=0 && alphaC(i,j)<0
                if i==1         % enforcing boundary conditions
                    newtheta(i,j) = newtheta(i,j) + A*theta(i,j);
                    newtheta(i,j+1) = newtheta(i,j+1) + B*theta(i,j);
                else
                    newtheta(i-1,j) = newtheta(i-1,j) + A*theta(i,j);
                    newtheta(i-1,j+1) = newtheta(i-1,j+1) + B*theta(i,j);
                end
                newtheta(i,j+1) = newtheta(i,j+1) + C*theta(i,j);
                newtheta(i,j) = newtheta(i,j) + D*theta(i,j);
            else  % alphaD(i,j)>=0 && alphaC(i,j)>0
                newtheta(i+1,j) = newtheta(i+1,j) + A*theta(i,j);
                newtheta(i+1,j+1) = newtheta(i+1,j+1) + B*theta(i,j);
                newtheta(i,j+1) = newtheta(i,j+1) + C*theta(i,j);
                newtheta(i,j) = newtheta(i,j) + D*theta(i,j);
            end
            
        end
    end
    
    % enforcing boundary conditions - perimeter set to zero
    newtheta(:,n-1) = newtheta(:,n-1) + newtheta(:,n); newtheta(:,n) = 0;
    newtheta(n-1,:) = newtheta(n-1,:) + newtheta(n,:); newtheta(n,:) = 0;
    
    deltaI = newtheta - theta;  % change due to individual dynamics (flux)
    
    e = zeros(n); f = zeros(n);
    for i = 1:n                 % extinction rate and fission rate for each i,j
        for j = 1:n
            x = (i-.5)*dx; y = (j-.5)*dy;
            e(i,j) = rho*e0*Groups/(1+phi*x+y)^2;	% extinction rate
            f(i,j) = rho*lambda*(x+y);              % fissioning rate
        end
    end
    
    deltaE = e .* theta * dt;   % change due to extinction and fissioning
    deltaF = f .* theta * dt;
    
    births = zeros(n);
    F = f.*theta;               % fissioning rate at current state of theta
    XY = [1:n]'*[1:n];          % matrix where XY(i,j)=i*j
    F_int = 2*F./XY;            % product of integrand and its operators, ready to be integrated.
    
    for i = 1:n                 % birth rate at i,j due to fissioning of larger groups
        for j = 1:n
            births(i,j) = (1-psi)*sum(sum(F_int(i:n,j:n)));
        end
    end
    
    if psi~=0                   % assortive splitting
        for i = 1:n
            for j = 1:n
%                 births(i,ceil(j/2)) = births(i,ceil(j/2)) + psi*F(i,j);
%                 births(1,ceil(j/2)) = births(1,ceil(j/2)) + psi*F(i,j);
                births(i,1) = births(i,1) + psi*F(i,j);
                births(1,j) = births(1,j) + psi*F(i,j);
            end
        end
    end
    deltaB = births * dt;       % change due to fissioning of larger groups
    
    tt=t*dt/freq;
    if tt == floor(tt)
        Theta(:,:,tt) = theta;  % write state variables to memory
        DeltaB(:,:,tt) = deltaB;
        Population(tt,:) = [t*dt Groups Coops Defectors Relate];
    end
    
    theta = theta + deltaB + deltaI - deltaE - deltaF; % update theta
    if min(min(theta)) < 0      % print to screen if error occurs. Indicative of 'dt' being too large
        disp(sprintf('At T=%.2f with r=%.2f, rho=%.2f, mu=%.2f, eps=%.2f, some theta(i,j) < 0. min(min(theta)) = %d',...
            t*dt, r, rho, mu, psi, min(min(theta))) ); 
        theta(theta<0)=0;      	% Force values to prevent program from crashing.
    end
      
end
tim=toc; disp(tim);             % display time elapsed

% optional saving protocol - MAKE SURE YOU HAVE A SUBFOLDER 'SimResults'
% str=sprintf('%s\\SimResults\\PDE_phi-%.0f_psi-%.0f_mu-%.0f_r-%.0f_rho-%.0f_dt-%.0f',...
str=sprintf('%s/SimResults/PDE_phi-%.0f_psi-%.0f_mu-%.0f_r-%.0f_rho-%.0f_dt-%.0f',...
    cd,100*phi,100*psi,100000*mu,100*r,100*rho,100*dt);
save(str);
