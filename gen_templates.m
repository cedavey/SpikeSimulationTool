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

duration = 200; % [msec]
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];
sampling_rate = 5000;
templates.initial_ap = 30;

% Calculate the refractory period
[templates.abs_refract_time, templates.abs_refract_index] = refract_period(tInit, xInit, const, duration, templates, sampling_rate);

% Calculate when second input will generate same output as previous output
% and generate transition templates
[templates.end_time, templates.end_index, templates.transition] = gen_trans(tInit, xInit, const, duration, sampling_rate, templates);

% Input function
Iapp = @(t) 10*exp(-((t - templates.initial_ap)*2).^2);

% Runs ODE
[t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);

% Interpolate data
[interp_t, interp_d] = interpolate(t, x(:,1), duration, sampling_rate);

% Generate single spike template
[templates.t, templates.d] = gen_template(interp_t, interp_d, 1, templates);

% Adjust the transition template lengths
templates = adj_templates(templates);

% Plot the simulation
plot_simulation(t, x, duration, Iapp, templates);

% Save the templates struct as a .mat file
save_templates(templates);

%% Interpolate data function

% This function will interpolate the data to the correct sampling rate
function [int_t, int_d] = interpolate(t, d, duration, sampling_rate)

int_t = 0 : 1/sampling_rate*1000 : duration;
int_d = interp1(t, d, int_t, 'spline');

end

%% Calculate the absolute refractory period

function [refract_time, refract_index] = refract_period(tInit, xInit, const, duration, templates, sampling_rate)

refract = 0;
i = templates.initial_ap + 1; % 1+ to the initial_ap that was defined

% Will loop until refract variable becomes true
while refract == 0
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i) * 2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);
    [int_t, int_d] = interpolate(t, x(:,1), duration, sampling_rate);
    
    % Use find peaks to find the first 2 peaks which will be the 2 APs
    pks = findpeaks(int_d);
    pks = pks(pks >= 0);
    if length(pks) >= 2 && diff(abs(pks([1,2]))) < 5    % If more than 2 peaks above 0 and the diff between the two is less than 5
        refract = 1;    % Change refract to be true to break while loop
        refract_time = i - templates.initial_ap;  % This refractory time is the time between the TWO INPUT SPIKES. AT this time is when the first following spike can occur
        refract_index = find(int_t >= i, 1, 'first') - find(int_t >= templates.initial_ap, 1, 'first');   % This refractory index is the values within the time matrix between the TWO INPUT SPIKES
    end
    if i > duration
        error('Unable to find refractory time due to lack of second spike');
        refract_time = NaN;
        refract_index = NaN;
    end
    i = i + 1;   
end
end

%% Generating the transition templates and calculating the end times/index

function [end_time, end_index, trans_temp] = gen_trans(tInit, xInit, const, duration, sampling_rate, templates)

trans_temp = {};
count = 0;

% Loops through the entirety of the duration if required in steps of 0.2ms
for i = templates.abs_refract_time + templates.initial_ap:1/sampling_rate * 1000:duration
    count = count + 1;
    
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i)*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, const);
    [int_t, int_d] = interpolate(t, x(:,1), duration, sampling_rate);
    [~, temp_d] = gen_template(int_t, int_d, 2);
    
    % Places data into cell arrays
    temp_d = {temp_d - temp_d(1,1)};
    trans_temp = [trans_temp temp_d];
    
    % Checks if the last 3 templates (in intervals of 5) have the same max 
    % value. If they do then it breaks out of loop and records value to adjust templates
    if (i > templates.abs_refract_time + templates.initial_ap + 3) ...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - 5}) + 0.05) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - 5}) - 0.05)) ...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - 10}) + 0.05) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - 10}) - 0.05))...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - 15}) + 0.05) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - 15}) - 0.05))
        trans_temp = {trans_temp{:, 1:count - 15}};
        end_time = i - 3;   % End time is absolute time (includes initial_ap time)
        end_index = find(int_t == i) - 15;  % Relative to INTERPOLATED index (includes initial_ap)
        break
    end
    
end
end

%% Template generator function

function [new_t, new_d] = gen_template(t, d, num_inputs, templates)

% Differentiate new_d to help find start and end index
slope = diff(d)./diff(t);

% If no AP is fired, do not generate template
if ~any(slope > 20)
    new_t = 0;
    new_d = 0;
    return
end

start_ap = 145; % Start AP index. This index's as 7 in the interpolated time data

% Will adjust end_ap depending on the number of inputs
if num_inputs == 1
    end_ap = templates.end_index;
    %end_ap = find(slope >= 0.05 | slope <= -0.05, 1, 'last');
else
    end_ap = length(t); % If there are 2 inputs, then it will adjust template later
    
end

% Adjusts the template
new_t = t(start_ap : end_ap)';
new_d = (d(start_ap : end_ap))';

% Find the initial_ap index


end

%% Adjust templates.transition templates

% Currently only works if there are MAXIMUM of 2 inputs
function templates = adj_templates(templates)

start_ap = 145;     % This is just a variable from the gen_templates
templates.initial_ap_index = find(templates.t == templates.initial_ap) - 1;  % - 1 because of weird matlab indexing

% This is the number of positions between input spike and end index which
% will correspond to second input spike and end index for the second input
% spike
% Assumes that the following spike dynamics are same as the first
index_bw_ap = templates.end_index - start_ap;

% Find the maximum value of the templates.transition cell array for
% normalisation
max_val = max(cellfun(@(x) max(x), templates.transition));

% Runs through the transition templates and shortens them to the correct
% length
for i = 1:size(templates.transition,2)
    templates.transition{i} = templates.transition{i}(1:templates.initial_ap_index ...
        + templates.abs_refract_index + i + index_bw_ap);
    
    % Normalises the templates amplitude so that the max peak is = 1
    temp = templates.transition{i};
    temp = temp./max_val;
    templates.transition{i} = temp;
end
end

%% Plot function

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

%% Save function
% Saves the struct template into a .mat file for SST to open
function save_templates(templates)

% This will then inverse the arrays to work with main code
templates.d = templates.d + abs(templates.d(1));

[file, path] = uiputfile(['templates' filesep 'templates_test3.mat'], 'Save file name'); % Will want to change the name
if file
    file_name = [path filesep file];
    save(file_name, 'templates');
else
    fprintf('User didnt''t chose a file location. Template was not saved\n');
end
end