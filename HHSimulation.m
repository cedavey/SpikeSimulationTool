% Modified from teaching material from Workshop 3 Neural Information Processing
% University of Mebourne
% Darren Tan 910828

% Function that generates and plots action potentials using the Hodgkin-Huxley model

clear all
close all

% Initialize channel constants
% All these also can be varied
const.vRest = -65;               % This can be varied later on
const.eNa   = 115 + const.vRest; % [mV]
const.gNa   = 120;
const.eK    = -12 + const.vRest; % [mV]
const.gK    = 36;
const.eLeak = 10.6 + const.vRest;% [mV]
const.gLeak = 0.3;
const.C     = 1.0;               % Capacitance of membrane

duration = 100; % [msec]
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];

% Iapp function
Iapp = @Iapp_func; % applied current injection (can be written as function of t)

% Run ODE
[t, x] = ode45('HHode', tInit, xInit, [], Iapp, const);

% Generate new template
[new_t, new_template] = gen_template(t, x(:,1), duration);

% Plot
plot_simulation(t, x, duration, Iapp, new_t, new_template);


%% Iapp function
function Iapp_out = Iapp_func(t)

% Bell curve function
Iapp_out = 10  * exp(-((t-15)*2).^2) + ...
           15 * exp(-((t-30)*3).^2) + ...
           20 * exp(-((t-70)/2).^2); 


% Constant function
%Iapp_out = 5;

% Sine function
%Iapp_out = 10*sin((t-10)/5);

% Removes negative numbers
%Iapp_out(Iapp_out < 0) = 0; 

end

%% Create new template function
function [new_t, new_template] = gen_template(t, data, duration)
fr = 5000; %sampling rate

% Interpolate data for constant sampling rate
new_t = 0 : 1/fr*1000 : duration;
new_template = interp1(t, data, new_t, 'spline');

% Extract only bits where there is significant change and removes points
% where membrane is resting 
% The > 0.2 and 0.05 difference values currently chosen to fit data
difference = abs(diff(new_template));
start_ap = new_t(find(difference > 0.2, 1, 'first')); % Start of AP [ms]
end_ap   = new_t(find(difference > 0.01, 1, 'last') + 1);   % End of AP [ms]

if (numel(start_ap) == 0) || (numel(end_ap) == 0) 
    % Checks if membrane is just resting the entire duration
    new_t        = 0; 
    new_template = 0;
else
    % Extract the relevant template section
    new_template((new_t < start_ap) | (new_t > end_ap)) = [];
    new_t((new_t < start_ap) | (new_t > end_ap)) = [];
end

end

%% PLOT function
function plot_simulation(t, x, duration, Iapp, new_t, new_template)
% Plot action potential
subplot(2,1,1);
hold on
yyaxis left
plot(t, x(:,1));
axis([0 duration -80 60]);
xlabel('time (ms)');
ylabel('voltage (mV)');

% Plot current applied
Iapp_ = Iapp(t);
if numel(Iapp_) == 1; Iapp_ = Iapp_*ones(size(t)); end % Changes a constant Iapp_ to an array for plotting

% Plot new template
plot(new_t, new_template, '-r');

yyaxis right
p = plot(t, Iapp_, '-');
p.Color(4) = 0.2; % Change transparency
axis([0, duration, -(4*max(abs(Iapp_))), (3*max(abs(Iapp_)))]);
xlabel('time (ms)');
ylabel('Current applied');

legend('AP voltage', 'Extracted template','Applied current');
hold off

% Plot the m, h and n variables
subplot(2,1,2);
hold on
plot(t, x(:,2));
plot(t, x(:,3));
plot(t, x(:,4));
ylim([-0.1 1.1]);
xlabel('Time (ms)');
ylabel('Activation');
title('Channel activations');
legend('m', 'h', 'n');
hold off

end