% Modified from teaching material from Workshop 3 Neural Information Processing
% University of Mebourne
% Darren Tan 910828

% Function that generates and plots action potentials using the Hodgkin-Huxley model

clear all
close all

% tic;
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
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);

    % Generate new template
    [templates.t, templates.d] = gen_template(t, x(:,1), duration);
    
    templates.refract_time = refract_period(tInit, xInit, const, duration);
    
    templates.transition = gen_trans(tInit, xInit, const, templates, duration);

    % Plot
    plot_simulation(t, x, duration, Iapp, templates);
    
    % Save the template struct as .mat file
    save_templates(templates);
    
% simulationTime = toc;
% fprintf('Time elapsed = %f\n', simulationTime);

%% Iapp function
function Iapp_out = Iapp_func(t)

% Bell curve function
Iapp_out = 10 * exp(-((t-30)*2).^2) + ...
           10 * exp(-((t-46)*2).^2);
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

% Spike template
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
function plot_simulation(t, x, duration, Iapp, template)
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
plot(template.t, template.d, '-r');

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
initial_ap = 30;
i = initial_ap + 1; % This corresponds to +1 of the initial input spike at t = 15
while refract == 0
    Iapp = @(t) 10 * exp(-((t-initial_ap)*2).^2) + 10 * exp(-((t-i)*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);
    
    % Use findpeaks to find the first 2 peaks which will be the 2 aps
    pks = findpeaks(x(:,1));
    pks = pks(pks>=0);
    if length(pks) >= 2 && diff(abs(pks([1,2]))) < 5
        refract = 1;
        refract_time = i - initial_ap;
    end
    if i > duration
        error('Unable to find refractory time due to lack of second spike');
        refract_time = NaN;
    end
    i = i + 1;
end
end

%% Generating templates for in between spikes for efficiency
% This function will generate the entire template with transition periods.
% This includes the two spikes that occur

% HOWEVER one issue that is not currently addressed (on the very very rare
% occasion) is if whilst generating the action potential, if somehow a value
% is == 0, then this will cause an irregularity in the template
function trans_temp = gen_trans(tInit, xInit, const, templates, duration)

% Initialising values
temp1 = [];
% trans_temp = {};
initial_ap = 30;    % This is the first input spike
% inx = 0;

% Iterate through this for the number of ms of templates you want
% It is currently at 10 templates from right after the abs refractory
% period
for i = templates.refract_time:templates.refract_time + 9
    Iapp = @(t) 10 * exp(-((t-initial_ap) * 2).^2) + 10 * exp(-((t-(initial_ap + i))*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);
    
    [~, temp2] = gen_template(t, x(:,1),duration);
    
%     temp2 = {temp2' - temp2(1,1)};
%     trans_temp = [trans_temp temp2];
    
    temp2 = temp2';
    
    % This little bit will concatenate different size matricies and fill in
    % the gaps with 0's
    [i1, j1] = ndgrid(1:size(temp1, 1), 1:size(temp1, 2));
    [i2, j2] = ndgrid(1:size(temp2, 1), (1:size(temp2, 2)) + size(temp1, 2));
    temp1 = accumarray([i1(:), j1(:); i2(:), j2(:)], [temp1(:); temp2(:)]);
end

% Find all the 0s (which is 99% of the time going to be just at the end of
% temp1 and change then to the initial value and then normalise the whole
% array
temp1(find(temp1 == 0)) = temp1(1,1);
trans_temp = temp1 - temp1(1,1);
end

%% Save function
% Saves the struct template into a .mat file for SST to open
function save_templates(templates)

templates.d = (templates.d + abs(templates.d(1)))';
templates.t = templates.t';

[file, path] = uiputfile(['templates' filesep 'template_test2.mat'], 'Save file name'); % Will want to change the name
if file
    file_name = [path filesep file];
    save(file_name, 'templates');
else
    fprintf('User didnt''t chose a file location. Template was not saved\n');
end
end