%% u_plot.m

close all
clc
clear all

% Read the data file
load Turbulent_BL;

%% Data File Content
% Laminar_BL.m contains u, i.e. the x-component of the velocity at five
% locations stored in five variables: U1l, U2l,...,U5l. The number 
% indicates the spatial streamwise position and "l" indicates laminar.
% Each column corresponds to a point in the wall normal direction
% (totally 46 points), and each row a sample (totally 4003). 
% The y-coordinates of the points are stored in the vector t and the times 
% at which the values are sampled are stored in tl. 

% This file shows how the data can be loaded into Matlab, ploted, and how
% various informative values can be calculated. The file is NOT complete 
% for the assignment. Study the commands and what happens, copy, and look
% for other useful commands in order to compute and vissualize according to
% the instructions.

% Reuse the file (copy to a new) for the turbulent flow (Turbulent_BL.mat),
% note though that for this there are data sampled at four streamwise 
% locations stored in U1t, U2t, U3t and U4t, and that the number of points 
% in the wall normal directions differs. Otherwise data can be handled in 
% the same way.
% 

%%
N = length(U1t(:,1));  %Number of samples obtained as the length of a (the first) column
V = 50; % Freestream velocity


%% Plot of instantantaneous u-velocity
% Store the size of the screen
scrsz = get(0,'ScreenSize');
figure('Position',[100 50 scrsz(3)/3 scrsz(4)*.75]) 
subplot(2,1,1)
% Plot the u-velocity signals at the first, second, third, fourth, and fifth
% point at the first spatial position vs time. 
plot(tt,U1t(:,1),'r',tt,U1t(:,2),'g',tt,U1t(:,3),'b',tt,U1t(:,4),'c',tt,U1t(:,5),'m');
legend('y=0.05 \mum','y=0.10\mum','y=0.20\mum','y=0.25\mum')
set(gca,'XLim',[0 tt(length(tt))],'YLim', [0 10])
title('Instantaneous x-velocity @ .3 m from LE')
xlabel('Time [s]');
ylabel('[m/s]');

subplot(2,1,2)
plot(tt,U1t(:,10),'r',tt,U1t(:,20),'g',tt,U1t(:,30),'b',tt,U1t(:,40),'c',tt,U1t(:,46),'m');
legend('y=122 \mum','y=554 \mum','y=1.7 mm','y=12 mm','y=30 mm')
set(gca,'XLim',[0 tt(length(tt))],'YLim', [0 60])
title('Instantaneous x-velocity @ .3 m from LE')
xlabel('Time [s]');
ylabel('[m/s]');

%% Calculation and plot of velocity profiles
for i=1:51
    % Time averages
    U1ta(i) = mean(U1t(:,i));
    U1tm= mean(U1ta);
    
    U2ta(i) = mean(U2t(:,i));        % Time average at point U2
    U2tm= mean(U2ta);
    
    U3ta(i) = mean(U3t(:,i));        % Time average at point U3
    U3tm= mean(U3ta);
    
    U4ta(i) = mean(U4t(:,i));        % Time average at point U4
    U4tm= mean(U4ta);

 
end

% Possible plot of the velocity profiles
figure('Position',[900 200 scrsz(3)/3 scrsz(4)/3]) 
plot(U1ta,yt,'-r')
xlabel('[m/s]')
ylabel('y [m]')
title('Boundary layer profile')

%% Calculation of boundary layer thickness

% Experimental 

blte=[0 0 0 0]
for j=1:1:51
    if(U1ta(j)<0.99*V)
        k=j;
    end
    m=k+1;
    blte(1)=yt(m);
end


for j=1:1:51
    if(U2ta(j)<0.99*V)
        k2=j;
    end
    m2=k2+1;
    blte(2)=yt(m2);
end

for j=1:51
    if(U3ta(j)<0.99*V)
        k3=j;
    end
    m3=k3+1;
    blte(3)=yt(m3);
end

for j=1:51
    if(U4ta(j)<0.99*V)
        k4=j;
    end
    m4=k4+1;
    blte(4)=yt(m4);
end


% Theoretical 
x=[0.10 0.15 0.20 0.25];
mu=1.7894.*10.^-5;
rho=1.225;
UT=[ U1tm U2tm U3tm U4tm];
Re = (rho.*V.*x).*mu.^-1;
blt= (0.3747.*x).*(Re.^-0.2);


 figure('Position',[700 300 scrsz(3)/3 scrsz(4)/3]) 
 plot(x,blt,'r',x,blte,'g') 
 legend('BL thickness,theory','BL thickness,experimental');
 xlabel('Distance from LE [m]');
 ylabel('Boundary layer [ - ]')
 ylim([0 0.009])
 title('Boundary layer thickness - Turbulent Flow')

%% Calculation of skin friction theoritical
cf_theory = 0.0583.*Re.^-0.2;

%% Calculation of skin friction experimental
du_dy=U1tm./yt;
%   %% Calculation of Shear stress
tau=[0 0 0 0];
%for tau
tau(1)=mu*mean(U1tm./yt);
tau(2)=mu*mean(U2tm./yt);
tau(3)=mu*mean(U3tm./yt);
tau(4)=mu*mean(U4tm./yt);


cf_e=[0 0 0 0];
cf_e(1)=tau(1).*(0.5.*rho*V.^2).^-1;
cf_e(2)=tau(2).*(0.5.*rho*V.^2).^-1;
cf_e(3)=tau(3).*(0.5.*rho*V.^2).^-1;
cf_e(4)=tau(4).*(0.5.*rho*V.^2).^-1;


 figure
 plot(x,cf_theory,'r*-',x,cf_e,'g*-') 
 legend('C_{f,theory}','C_{f,experimental}');
 ylim([0 18*10^-3])
 xlabel('Distance from LE [m]');
 ylabel('Skin Friction Coefficient [ C_f]')
 title('Skin Friction Coefficient - Turbulent Flow')