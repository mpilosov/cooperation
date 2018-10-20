% This code will allow one to make animations of the Theta state variable 
% as seen in figure 2 and in the supplementary information. 
% Author: Michael Pilosov
% This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.

%% INITIALIZE

clear all; close all;
r = .3;
rho = 1;
dt = .04;
mu = 0.0025;
phi = 1;
psi = 0;
eps = 0;

% assumes files are saved to folder 'SimResults' inside current directory
str = sprintf('%s/SimResults/PDE_phi-%.0f_psi-%.0f_mu-%.0f_r-%.0f_rho-%.0f_dt-%.0f',...
    cd,100*phi,100*psi,100000*mu,100*r,100*rho,100*dt);
load(str);

%% QUICKLY DISPLAY RESULTS (useful for testing)

for i=1:size(Theta,3); 
    surf(dx*[1:n],dy*[1:n],Theta(:,:,i)','EdgeColor', 'flat','EdgeLighting','gouraud'); 
    view(0,90); axis([1 Xmax 1 Ymax]); 
    pause(.01); 
end

%% CREATE ANIMATIONS - AVI (default) or GIF

len = size(Theta,3); % default length. 
skip = 1; % skip frames in Theta to make animation play faster
rep = 15; % repeat final frame at the end of animation

% --- format .gif
% delay=.04;
% colors=64;

% --- format .avi
avi_object = VideoWriter(sprintf('%s.avi',str)); % saves with same filename as loaded
avi_object.FrameRate = 30;
open(avi_object);

for i=[1:skip:len len*ones(1,rep)]
    
    str2=sprintf('     r = %.1f     rho = %.1f \n           mu = %.4f \n  psi = %.2f phi = %1.2f \n\n  T = %.0f \n  G = %6.3f \n  C = %6.3f \n  D = %6.3f \n  C/P = %2.3f \n  R = %.3f', r, rho, mu, psi, phi, Population(i), Population(i,2), Population(i,3), Population(i,4), Population(i,3)/( Population(i,3)+Population(i,4) ), Population(i,5) ); 
    surf(dx*[1:n],dy*[1:n],Theta(:,:,i)','EdgeColor', 'flat','EdgeLighting','gouraud'); 
    axis([1 Xmax 1 Ymax]);
    view(0,90); 
    legend(str2);
    pause(.01);
    drawnow;
    frame=getframe;
    
% --- format gif
%     im=frame2im(frame);
%     [imind,map]=rgb2ind(im,colors);
%     if i==1
%         imwrite(imind,map, sprintf('%s.gif',str2),'DelayTime',delay,'LoopCount',inf);
%     else
%         imwrite(imind,map, sprintf('%s.gif',str2), 'DelayTime',delay, 'WriteMode', 'append');
%     end

% --- format avi
    writeVideo(writerObj,frame)
end

close(avi_object); % close out avi object and write to file

disp('done');
close

%% Author's Note

% To recreate figures 3 & 4, one must run multiple tests using a 
% for-loop (lines 50-237) in Synergy.m, then use data stored in the 
% variable 'Population' to determine the time at which XX% cooperators
% occurs. To automate the process, one can put a stop criteria inside
% Synergy.m's for loop, set T sufficiently high, and store the appropriate 
% values of interest (t*dt, variables r, rho, mu, etc.) in another
% variable. This is what was done to create the figures in the paper.
% Direct technical questions to mpilosov@gmail.com.