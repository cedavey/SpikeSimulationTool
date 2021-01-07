% Modified from teaching material from Workshop 3 Neural Information Processing
% University of Mebourne
% Darren Tan 910828

% Function that generates and plots action potentials using the Hodgkin-Huxley model

clear all
close all

tic;
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
    
    refract_time = refract_period(tInit, xInit, const, duration);

    % Plot
    plot_simulation(t, x, duration, Iapp, new_t, new_template);
    
    % Calculate the refractory period
%     [t2, x2] = ode45('HHode', tInit, xInit, [], Iapp(2), const);
%     plot_simulation(t2, x2, duration, Iapp(2), [], []);
    
simulationTime = toc;
fprintf('Time elapsed = %f\n', simulationTime);

%% Iapp function
function Iapp_out = Iapp_func(t)

% Bell curve function
Iapp_out = 10 * exp(-((t-15)*2).^2);
           %10 * exp(-((t-28)*2).^2);
           %20 * exp(-((t-70)/2).^2); 


% % Constant function
% Iapp_out = 0;
% for i = 15:50
%     Iapp_out = Iapp_out + 10 * exp(-((t - i)*2).^2);
% end

% Sine function
%Iapp_out = 10*sin((t-10)/5);

% Removes negative numbers
%Iapp_out(Iapp_out < 0) = 0; 

end

%% Create new template function
function [new_t, new_template] = gen_template(t, data, duration)
fr = 5000; %sampling rate

% Interpolate data for constant sampling rate
new_t = 0 : 1/fr*1000 : duration;   % Time vector
new_template = interp1(t, data, new_t, 'spline');   % Value vector

% Differentiate the new_template to find start and end index
slope = diff(new_template)./diff(new_t);

% Do not extract template if no action potential is fired
if ~any(slope > 20) 
    new_t = 0; 
    new_template = 0; 
    return 
end

% If there is a large change in the first 10 ms, ignore these assuming the
% system is going back into equilibrium
start_ap = find(slope(10:end) >= 1, 1, 'first') + 10;   % Start AP index
end_ap = find(slope >=0.05 | slope <= -0.03, 1, 'last');    % End AP index

% Extracting only the relevant template
new_t = new_t(start_ap - 5:end_ap);
new_template = new_template(start_ap - 5:end_ap);

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

%% Calculate refractory period function

function refract_time = refract_period(tInit, xInit, const, duration)

refract = 0;
i = 16; % This corresponds to +1 of the initial input spike at t = 15
while refract == 0
    Iapp = @(t) 10 * exp(-((t-15)*2).^2) + 10 * exp(-((t-i)*2).^2);  % Need to figure a way to pass the Iapp equation into ode45 (with the t variable) and vary the equation whilst at it to iterate
    [t, x] = ode45('HHode', tInit, xInit, [], Iapp, const);
    
    ap_max1 = find(x == max(x([1:320],1))); % Will need to change the 310 to a variable that can vary (since it will iterate)
    
    ap_max2 = find(x == max(x([320:end] , 1))); % Need to change the 310 to a variable as well
    
    % Calculates if the difference between the two action potentials is
    % less than 3, then it will consider it an action potential and change
    % refract to 1
    if abs(x(ap_max1,1) - x(ap_max2,1)) < 3
        refract = 1;
        refract_time = i - 15;   % This is in ms
    end
    i = i + 1;
    
    % This will make sure this function does not exceed duration of
    % simulation
    if i > duration
        error('Unable to find refractory time due to lack of second spike');
        refract_time = NaN;
    end
end

end    
    
    
