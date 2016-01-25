%------------------------------
% GEOL 5700 
% Model the geotherm using finite difference methods (FTCS)
%------------------------------
% Franklin Hinckley
%------------------------------

%% Cleanup
clearvars
close all
clc

% Define parameters
k = 2.5; % [W/mK]
kappa = 1e-6; % [m^2/sec]
rhoCp = k/kappa;
Qm = 45/1000; % [W/m^2]

% Set boundary conditions
Pyr = 3600*24*365.25; % pne year [sec]
Pday = 86400; % one day [sec]
dT_yr = 15; % variation in temperature per year [deg C]
dT_day = 10; % variation in temperature per day [deg C]
mT_yr = -10; % mean yearly temperature [deg C]

% Set initial conditions
dz = 0.25;
zmax = 10*sqrt(kappa*Pyr/pi); % ten characteristic depths so gradient at depth is just mantle
z = 0:dz:zmax;
Ts = mT_yr;
T = (Qm./k).*z + Ts(1); % start simulation with linear T(z) from mantle

% Define time vector
dt = 0.75*(dz^2/(2*kappa)); % 0.75 gives margin from what is req'd for stability
tProp = 8*(3600*24*365.25);
t = 0:dt:tProp;

Ts = mT_yr + ...
        dT_yr.*sin(((2*pi.*t)/Pyr)) + dT_day.*sin(((2*pi.*t)/Pday));

%% Define analytic bounds
zStar_yr = sqrt(kappa*Pyr/pi);
zStar_day = sqrt(kappa*Pday/pi);
boundLow = T - dT_yr.*exp(-z./zStar_yr) - dT_day.*exp(-z./zStar_day);
boundHigh = T + dT_yr.*exp(-z./zStar_yr) + dT_day.*exp(-z./zStar_day);


%% Set initial temperature profile
T = T + ...
        dT_yr.*exp(-z./zStar_yr).*sin(((2*pi*(0))/Pyr) - (z./zStar_yr)) + ...
        dT_day.*exp(-z./zStar_day).*sin(((2*pi*(0))/Pday) - (z./zStar_day));

%% Main loop
TS = zeros(length(z),length(t));
TS(:,1) = T;
for ii = 1:length(t)
    % Get surface temperature at current time
    TsC = Ts(ii);
    
    % Determine temperature gradient at depth
    dTdzmax = Qm/k; %[C/m]

    % Compute temperature gradient
    dTdz = [diff(T)./dz dTdzmax];

    % Compute heat flux
    q = -k*dTdz;

    % Compute temperature rate
    dqdz = diff(q)/dz; 
    dTdt = (-1/(rhoCp))*dqdz;

    % Update temperature
    T(2:end) = T(2:end) + dTdt*dt;
    T(1) = TsC;
    
    % Save output
    TS(:,ii) = T;
end

%% Plot results
figure
hold on
plot(TS(:,7*(3600*24*365.25)/dt:end),z)
plot(boundLow,z,'-k','Linewidth',2)
plot(boundHigh,z,'-k','Linewidth',2)
plot([0 0],[min(z) max(z)],'--r')
%plot([min(T(:))-2 max(T(:))+2],[z(zInd) z(zInd)],'--k')
hold off
%text(min(boundLow)-1,z(zInd)-0.5,zStr)
xlim([min(boundLow)-2 max(boundHigh)+2])
ylim([0 zmax])
set(gca,'Ydir','reverse')
title('Geotherm','Fontsize',14)
ylabel('Depth [m]','Fontsize',12)
xlabel('Temp [\circC]','Fontsize',12)

% Animation
% figure
% for ii = 1:length(t)
%     plot(TS(:,ii),z,'b')
%     hold on
%     plot(boundLow,z,'-k')
%     plot(boundHigh,z,'-k')
%     plot([0 0],[min(z) max(z)],'r')
%     hold off
%     xlim([min(TS(:))-2 max(TS(:))+2])
%     ylim([0 zmax])
%     set(gca,'Ydir','reverse')
%     title('Geotherm','Fontsize',14)
%     ylabel('Depth [m]','Fontsize',12)
%     xlabel('Temp [\circC]','Fontsize',12)
% %    M(:,ii) = getframe(gcf);
%     pause(0.01)    
% end


