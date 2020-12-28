% Modified from teaching material from Workshop 3 Neural Information Processing
% University of Mebourne
% Darren Tan 910828

% Function that generates and plots action potentials using the Hodgkin-Huxley model

clear all
close all

% Initialize constants
duration = 100; % duration of simulation in msec
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];

% Iapp function
Iapp     = @Iapp_func; % applied current injection (can be written as function of t)

% Run ODE
[t, x] = ode45('HHode', tInit, xInit, [], Iapp);

% Generate new template
[new_t, new_template] = gen_template(t, x(:,1), duration);

% Plot
plot_simulation(t, x, duration, Iapp, new_t, new_template);


%% PLOT function
function plot_simulation(t, x, duration, Iapp, new_t, new_template)
% Plot action potential
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
plot(new_t, new_template, 'r');

yyaxis right
p = plot(t, Iapp_, '-');
p.Color(4) = 0.2; % Change transparency
axis([0, duration, -(4*max(abs(Iapp_))), (3*max(abs(Iapp_)))]);
xlabel('time (ms)');
ylabel('Current applied');

legend('AP voltage', 'New extracted template','Applied current');

end

%% Iapp function
function Iapp_out = Iapp_func(t)

Iapp_out = 10*exp(-(2*t-50).^2); % Bell curve

%Iapp_out(Iapp_out < 0) = 0; 

end

%% Create new template function
function [new_t, new_template] = gen_template(t, data, duration)
sr = 5000; %sampling rate

% Interpolate data for constant sampling rate
new_t = 0 : 1/sr*1000 : duration;
new_template = interp1(t, data, new_t, 'spline');

% Extract only bits where there is an action potential and removes points
% where membrane is resting (values currently hardcoded and needs to be updated for each
% new spike generation)
start_ap = 23.2; % Start of AP [ms]
end_ap   = 51;   % End of AP [ms]

new_template((new_t < start_ap) | (new_t > end_ap)) = [];
new_t((new_t < start_ap) | (new_t > end_ap)) = [];

end


