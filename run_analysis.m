%% Pipeline surge caused by downstream valve closure
% Representative industrial cooling-water case. The data are synthetic but
% physically plausible. This script uses core MATLAB and moc_solver.m only.

clear; close all; clc

%% Physical data
g = 9.81;                              % m/s^2

% Water properties at approximately 25 degC
fluid.rho = 997.0;                     % kg/m^3
fluid.cp = 4180.0;                     % J/(kg K)
fluid.mu = 0.890e-3;                   % Pa s
fluid.bulkModulus = 2.20e9;            % Pa
fluid.vapourPressure = 3.17e3;         % Pa absolute
fluid.atmosphericPressure = 101325;     % Pa absolute

% Commercial-steel pipe, approximately DN 300 / NPS 12 Schedule 40
pipe.L = 100.0;                        % m
pipe.OD = 0.3238;                      % m
pipe.thickness = 0.00953;              % m
pipe.D = pipe.OD - 2*pipe.thickness;   % m, internal diameter
pipe.E = 200e9;                        % Pa
pipe.roughness = 4.5e-5;               % m
pipe.area = pi*pipe.D^2/4;             % m^2
pipe.relativeRoughness = pipe.roughness/pipe.D;

% Receiving-side head and steady local-loss coefficients
headAfterValve = 8.0;                  % m
K_exchanger = 40;
K_fittings = 10;
K_valve = 10;                          % fully open valve
K_total = K_exchanger + K_fittings + K_valve;

%% Cooling-water duty and steady operating point
thermalDuty = 4.2e6;                   % W
waterTemperatureRise = 10;             % K
Qdesign = thermalDuty/(fluid.rho*fluid.cp*waterTemperatureRise);

reynoldsNumber = @(Q) fluid.rho.*(Q./pipe.area).*pipe.D./fluid.mu;
frictionFactor = @(Q) 0.25./log10( ...
    pipe.relativeRoughness/3.7 + 5.74./reynoldsNumber(Q).^0.9).^2;
velocityHead = @(Q) (Q./pipe.area).^2/(2*g);

systemHead = @(Q) headAfterValve + ...
    (frictionFactor(Q)*pipe.L/pipe.D + K_total).*velocityHead(Q);

% Representative centrifugal-pump curve
shutoffHead = 23;                      % m
designHead = 17;                       % m at Qdesign
pumpHead = @(Q) shutoffHead - ...
    (shutoffHead-designHead)*(Q/Qdesign).^2;

Q0 = fzero(@(Q) pumpHead(Q)-systemHead(Q),Qdesign);
upstreamHead0 = pumpHead(Q0);
valveInletHead0 = headAfterValve + K_valve*velocityHead(Q0);

Qplot = linspace(0,0.18,100);
figure('Color','w','Position',[120 120 760 430]);
plot(Qplot,pumpHead(Qplot),'LineWidth',1.6); hold on
plot(Qplot,systemHead(Qplot),'LineWidth',1.6);
plot(Q0,upstreamHead0,'ko','MarkerFaceColor','k');
grid on
xlabel('Flow rate, Q [m^3/s]');
ylabel('Head [m]');
title('Steady operating point');
legend('Pump curve','System curve','Operating point','Location','best');

%% Wave speed and transient-model inputs
waveSpeed = sqrt((fluid.bulkModulus/fluid.rho)/(1 + ...
    fluid.bulkModulus*pipe.D/(pipe.E*pipe.thickness)));
waveRoundTripTime = 2*pipe.L/waveSpeed;

% The control valve is represented explicitly at the downstream boundary.
% Therefore, its K value is excluded from the distributed equivalent factor.
fEquivalent = frictionFactor(Q0) + ...
    (K_exchanger + K_fittings)*pipe.D/pipe.L;

operatingPoint = table( ...
    ["Thermal design flow"; "Operating flow"; "Mean velocity"; ...
     "Reynolds number"; "Darcy friction factor"; ...
     "Equivalent distributed factor"; "Upstream head"; ...
     "Valve-inlet head"; "Receiving-side head"; "Wave speed"; ...
     "Wave round-trip time"], ...
    [Qdesign*3600; Q0*3600; Q0/pipe.area; reynoldsNumber(Q0); ...
     frictionFactor(Q0); fEquivalent; upstreamHead0; valveInletHead0; ...
     headAfterValve; waveSpeed; waveRoundTripTime], ...
    ["m^3/h"; "m^3/h"; "m/s"; "-"; "-"; "-"; "m"; "m"; ...
     "m"; "m/s"; "s"], ...
    'VariableNames',{'Quantity','Value','Unit'});
disp(operatingPoint)

%% MOC configuration and Joukowsky verification
gridSize = 100;                        % computational pipe reaches
endTime = 50*waveRoundTripTime;
linearValveLaw = @(s) max(0,1-min(max(s,0),1));

runTransient = @(closureTime) moc_solver(g,pipe,Q0,upstreamHead0, ...
    valveInletHead0,fEquivalent,headAfterValve,waveSpeed,closureTime, ...
    gridSize,endTime,linearValveLaw);

instantaneousResult = runTransient(0);
deltaPJoukowsky = fluid.rho*waveSpeed*Q0/pipe.area;
deltaPMOC = fluid.rho*g*(instantaneousResult.head(end,2) - ...
    instantaneousResult.head(end,1));
verificationError = 100*abs(deltaPMOC-deltaPJoukowsky)/deltaPJoukowsky;

verification = table(deltaPJoukowsky/1e5,deltaPMOC/1e5,verificationError, ...
    'VariableNames',{'JoukowskyRise_bar','MOCFirstStepRise_bar', ...
    'RelativeDifference_percent'});
disp(verification)

%% Linear valve-closing-time study
closureTimes = 0:1:5;                  % s
results = cell(numel(closureTimes),1);

for i = 1:numel(closureTimes)
    results{i} = runTransient(closureTimes(i));

    % The pipe is treated as horizontal with elevation datum z = 0, so the
    % computed piezometric head is converted directly to gauge pressure.
    results{i}.pressureGauge = fluid.rho*g*results{i}.head;
    results{i}.pressureAbsolute = results{i}.pressureGauge + ...
        fluid.atmosphericPressure;
    results{i}.pressureBarG = results{i}.pressureGauge/1e5;
    results{i}.pressureBarA = results{i}.pressureAbsolute/1e5;
end

peakBarA = zeros(numel(results),1);
minimumBarA = zeros(numel(results),1);
vapourMarginBar = zeros(numel(results),1);
peakHoopStressMPa = zeros(numel(results),1);
vapourScreen = strings(numel(results),1);

for i = 1:numel(results)
    peakBarA(i) = max(results{i}.pressureBarA,[],'all');
    minimumBarA(i) = min(results{i}.pressureBarA,[],'all');
    vapourMarginBar(i) = minimumBarA(i) - fluid.vapourPressure/1e5;

    peakGaugePressure = max(results{i}.pressureGauge,[],'all');
    peakHoopStressMPa(i) = peakGaugePressure*pipe.D/ ...
        (2*pipe.thickness)/1e6;

    if vapourMarginBar(i) > 0
        vapourScreen(i) = "ABOVE VAPOUR";
    else
        vapourScreen(i) = "CAVITATION ONSET";
    end
end

closureComparison = table(closureTimes(:),peakBarA,minimumBarA, ...
    vapourMarginBar,peakHoopStressMPa,vapourScreen, ...
    'VariableNames',{'ClosureTime_s','PeakPressure_barA', ...
    'MinimumPressure_barA','VapourMargin_bar', ...
    'PeakHoopStress_MPa','VapourScreen'});
disp(closureComparison)

%% Valve-inlet pressure histories
% Instantaneous closure is reserved for verification; practical timed
% closures are shown here so that the design cases remain readable.
pVapourGaugeBar = ...
    (fluid.vapourPressure-fluid.atmosphericPressure)/1e5;

figure('Color','w','Position',[120 120 1000 450]); hold on
for i = 2:numel(results)
    plot(results{i}.t,results{i}.pressureBarG(end,:),'LineWidth',1.3);
end
yline(pVapourGaugeBar,'--k','Vapour-pressure threshold', ...
    'LabelHorizontalAlignment','right');
grid on
xlabel('Time [s]');
ylabel('Valve-inlet pressure [bar(g)]');
title('Pressure at the control valve with linear closure');
legend(arrayfun(@(t) sprintf('t_{close} = %g s',t),closureTimes(2:end), ...
    'UniformOutput',false),'Location','bestoutside');
exportgraphics(gcf,'valve_pressure_history.jpg','Resolution',200);

%% Spatial pressure envelopes
figure('Color','w','Position',[120 120 1100 430]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on
for i = 2:numel(results)
    plot(results{i}.x,max(results{i}.pressureBarA,[],2),'LineWidth',1.2);
end
grid on
xlabel('Distance from source [m]');
ylabel('Maximum pressure [bar(a)]');
title('Maximum-pressure envelope');

nexttile; hold on
for i = 2:numel(results)
    plot(results{i}.x,min(results{i}.pressureBarA,[],2),'LineWidth',1.2);
end
yline(fluid.vapourPressure/1e5,'--k','Vapour-pressure threshold', ...
    'LabelHorizontalAlignment','right');
grid on
xlabel('Distance from source [m]');
ylabel('Minimum pressure [bar(a)]');
title('Minimum-pressure envelope');
legend(arrayfun(@(t) sprintf('t_{close} = %g s',t),closureTimes(2:end), ...
    'UniformOutput',false),'Location','best');

%% Engineering interpretation
firstAcceptable = find(vapourMarginBar > 0,1,'first');
if isempty(firstAcceptable)
    fprintf('\nNo tested closure time remains above vapour pressure.\n');
else
    fprintf(['\nThe shortest tested closure retaining a positive vapour-' ...
        'pressure margin is %.1f s.\n'],closureTimes(firstAcceptable));
end
fprintf(['Sub-vapour single-phase results indicate possible column ' ...
    'separation; they are not physical pressure predictions.\n']);
