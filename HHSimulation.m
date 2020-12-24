% Modified from teaching material from Workshop 3 Neural Information Processing
% University of Mebourne
% Darren Tan 910828

% Function that generates and plots action potentials using the Hodgkin-Huxley model

close all

% Initialize constants
duration = 100; % duration of simulation in msec
Iapp     = @(t) 10*square(t/10 + 100) + 10; % applied current injection (can be written as function of t)
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];

% Run ODE
[t, x] = ode45('HHode', tInit, xInit, [], Iapp);

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

yyaxis right
plot(t, Iapp_, '-');
axis([0, duration, -(4*max(abs(Iapp_))), (3*max(abs(Iapp_)))]);
xlabel('time (ms)');
ylabel('Current applied');

legend('AP voltage', 'Applied current');
