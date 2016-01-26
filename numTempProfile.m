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

%% Cape Thompson data
% Load and parse data
load('cape_thompson.dat');
CT_z = cape_thompson(:,1)';
CT_T = cape_thompson(:,2)';

% Year of data set
t_CT = 1961;

%% Fit step change to Cape Thompson data
% Estimate original surface temperature by extrapolating slope at depth
p_depth = polyfit(CT_z(11:end),CT_T(11:end),1);
CT_TS0 = polyval(p_depth,0);

% Estimate original temperature profile for steady state
CT_T0 = polyval(p_depth,[0 CT_z]);

% Assert current surface temperature is the first data point
%   This assumes that temperature has been constant at current surface
%   temperature for some time, daily and yearly variations will have damped
%   by this depth (27m, 8.5 length scales for zStar_year)
CT_TS = CT_T(1);
CT_T = [CT_T(1) CT_T];

% Plot Cape Thompson data with linear fit
figure
hold on
plot(CT_T,[0 CT_z])
plot(CT_T0,[0 CT_z],'-r')
hold off
xlim([min(CT_T)-1 max(CT_T)+1])
ylim([0 max(CT_z)])
set(gca,'Ydir','reverse')
title('Cape Thompson: Original Geotherm','Fontsize',14)
ylabel('Depth [m]','Fontsize',12)
xlabel('Temp [\circC]','Fontsize',12)
legend('Measured Temperature','Linear Fit')

% Compute magnitude of step
Tstep = (CT_TS - CT_TS0); % [deg C] 

% Define non-dimensional length
L = max(CT_z);
xstar = CT_z/L;

% Define vector of times since step
tvec_step = 1:1:100;
tvec_step = tvec_step*(86400*365.25); % [sec]

% Determine length step and time step
dz = min(diff(CT_z));
dt = 0.1*(dz^2/(2*kappa));

% Initialize temperatures
T = CT_T0;

% Initialize error
minErr = 1e50;

% Loop over vector of times since step
for jj = 1:length(tvec_step)
    % Compute expected temperature profile
    t = 0:dt:tvec_step(jj);
    for ii = 1:length(t)
        % Get surface temperature at current time
        TsC = CT_TS;

        % Determine temperature gradient at depth
        dTdzmax = p_depth(1); %[C/m]

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
    end
    
    % Compute error (sum of squared errors)
    err = sum((T-CT_T).^2);
    
    % Check against minimum error
    if err < minErr
        t_min = tvec_step(jj);
        T_est = T;
        minErr = err;
    end    
end

% Print estimated time since step
fprintf('Time since step: %3.0f years',t_min/(86400*365.25))
fprintf('\n')

% Plot temperature profiles
figure
hold on
plot(CT_T,[0 CT_z])
plot(T_est,[0 CT_z],'r')
hold off
xlim([min(CT_T)-1 max(CT_T)+1])
ylim([0 max(CT_z)])
set(gca,'Ydir','reverse')
title('Cape Thompson: Best Fit','Fontsize',14)
ylabel('Depth [m]','Fontsize',12)
xlabel('Temp [\circC]','Fontsize',12)
legend('Measured Temperature','Estimated Temperature')

