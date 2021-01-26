function [templates] = HHSim(parameters)
% 
% duration = 100; % [msec]
% tInit    = [0 duration];
% xInit    = [-65; 0.052; 0.059; 0.317];
% 
% % Iapp function
% Iapp = @Iapp_func; % applied current injection (can be written as function of t)
% 
% % Run ODE
% [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
% 
% % Generate new template
% [templates.t, templates.d] = gen_template(t, x(:,1), duration);
% 
% templates.refract_time = refract_period(tInit, xInit, parameters, duration);
% 
% templates.transition = gen_trans(tInit, xInit, parameters, templates, duration);
% 
% templates.HHparameters = parameters;
% 
% end
% 
% %% Iapp function
% function Iapp_out = Iapp_func(t)
% 
% % Bell curve function
% Iapp_out = 15 * exp(-((t-30)*2).^2);
% 
% end
% 
% 
% %% Create new template function
% function [new_t, new_template] = gen_template(t, data, duration)
% fr = 5000; %sampling rate
% 
% % Interpolate data for constant sampling rate
% new_t = 0 : 1/fr*1000 : duration;   % Time vector
% new_template = interp1(t, data, new_t, 'spline');   % Value vector
% 
% % Spike template
% % Differentiate the new_template to find start and end index
% slope = diff(new_template)./diff(new_t);
% 
% % Do not extract template if no action potential is fired
% if ~any(slope > 20) 
%     new_t = 0; 
%     new_template = 0; 
%     return 
% end
% 
% % If there is a large change in the first 10 ms, ignore these assuming the
% % system is going back into equilibrium
% start_ap = find(slope(10:end) >= 1, 1, 'first') + 10;   % Start AP index
% end_ap = find(slope >=0.05 | slope <= -0.03, 1, 'last');    % End AP index
% 
% % Extracting only the relevant template
% new_t = new_t(start_ap - 5:end_ap);
% new_template = new_template(start_ap - 5:end_ap);
% 
% % Change to column vector and normalize to 0
% new_t = new_t';
% new_template = (new_template + abs(new_template(1)))';
% 
% end
% 
% 
% %% Calculate refractory period function
% function refract_time = refract_period(tInit, xInit, parameters, duration)
% 
% refract = 0;
% initial_ap = 30;
% i = initial_ap + 1; % This corresponds to +1 of the initial input spike at t = 15
% while refract == 0
%     Iapp = @(t) 10 * exp(-((t-initial_ap)*2).^2) + 10 * exp(-((t-i)*2).^2);
%     [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
%     
%     % Use findpeaks to find the first 2 peaks which will be the 2 aps
%     pks = findpeaks(x(:,1));
%     pks = pks(pks>=0);
%     if length(pks) >= 2 && diff(abs(pks([1,2]))) < 5
%         refract = 1;
%         refract_time = i - initial_ap;
%     end
%     if i > duration
%         error('Unable to find refractory time due to lack of second spike');
%         refract_time = NaN;
%     end
%     i = i + 1;
% end
% end
% 
% %% Generating templates for in between spikes for efficiency
% % This function will generate the entire template with transition periods.
% % This includes the two spikes that occur
% 
% % HOWEVER one issue that is not currently addressed (on the very very rare
% % occasion) is if whilst generating the action potential, if somehow a value
% % is == 0, then this will cause an irregularity in the template
% function trans_temp = gen_trans(tInit, xInit, parameters, templates, duration)
% 
% % Initialising values
% trans_temp = {};    % Used for if want to store data in cell arrays
% initial_ap = 30;    % This is the first input spike
% 
% % Iterate through this for the number of ms of templates you want
% % It is currently at 10 templates from right after the abs refractory
% % period
% for i = templates.refract_time:templates.refract_time + 9
%     Iapp = @(t) 10 * exp(-((t-initial_ap) * 2).^2) + 10 * exp(-((t-(initial_ap + i))*2).^2);
%     [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
%     
%     [~, temp2] = gen_template(t, x(:,1),duration);
%     
%     % This uses cell arrays instead
%     temp2 = {temp2' - temp2(1,1)};
%     trans_temp = [trans_temp temp2];
%     
% end
% end



%%

duration = 200; % [msec]
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];
% sampling_rate = 5000;   % Add this into parameters
templates.initial_ap = 30;

% Calculate the refractory period
[templates.abs_refract_time, templates.abs_refract_index] = refract_period(tInit, xInit, duration, templates, parameters);

% Calculate when second input will generate same output as previous output
% and generate transition templates
[templates.end_time, templates.end_index, templates.transition] = gen_trans(tInit, xInit, duration, templates, parameters);

% Input function
Iapp = @(t) 10*exp(-((t - templates.initial_ap)*2).^2);

% Runs ODE
[t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);

% Interpolate data
[interp_t, interp_d] = interpolate(t, x(:,1), duration, parameters);

% Generate single spike template
[templates.t, templates.d] = gen_template(interp_t, interp_d, 1, templates);

% Adjust the transition template lengths
templates = adj_templates(templates);

end

%% Interpolate data function

% This function will interpolate the data to the correct sampling rate
function [int_t, int_d] = interpolate(t, d, duration, parameters)

int_t = 0 : 1/parameters.sampling_rate*1000 : duration;
int_d = interp1(t, d, int_t, 'spline');

end

%% Calculate the absolute refractory period

function [refract_time, refract_index] = refract_period(tInit, xInit, duration, templates, parameters)

refract = 0;
i = templates.initial_ap + 1; % 1+ to the initial_ap that was defined

% Will loop until refract variable becomes true
while refract == 0
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i) * 2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    [int_t, int_d] = interpolate(t, x(:,1), duration, parameters);
    
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

function [end_time, end_index, trans_temp] = gen_trans(tInit, xInit, duration, templates, parameters)

trans_temp = {};
count = 0;

% Loops through the entirety of the duration if required in steps of 0.2ms
for i = templates.abs_refract_time + templates.initial_ap:1/parameters.sampling_rate * 1000:duration
    count = count + 1;
    
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i)*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    [int_t, int_d] = interpolate(t, x(:,1), duration, parameters);
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
