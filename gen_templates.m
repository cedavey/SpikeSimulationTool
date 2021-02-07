function [templates] = gen_templates(parameters)

duration = 200; % [msec]
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];
% sampling_rate = 5000;   % Add this into parameters
templates.initial_ap = 30;
templates.initial_ap_index = find(0:1/parameters.sampling_rate*1000:duration >= templates.initial_ap, 1, 'first');

% Progress bar
w = waitbar(0, 'Generating templates...');

% Store parameters into the templates struct for reference
templates.parametersUsed = store(parameters);

% Calculate the refractory period
[templates.abs_refract_time, templates.abs_refract_index, templates.rel_refract_time, templates.rel_refract_index] = refract_period(tInit, xInit, duration, templates, parameters, w);

try w = waitbar(2/8, w); catch, delete(w); error('Manually stopped'); end

% Calculate when second input will generate same output as previous output
% and generate transition templates
templates.transition = gen_trans(tInit, xInit, duration, templates, parameters);

try w = waitbar(3/8, w); catch, delete(w); error('Manually stopped'); end

% Input function
Iapp = @(t) 10*exp(-((t - templates.initial_ap)*2).^2);

try w = waitbar(4/8, w); catch, delete(w); error('Manually stopped'); end

% Runs ODE
[t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);

try w = waitbar(5/8, w); catch, delete(w); error('Manually stopped'); end

% Interpolate data
[interp_t, interp_d] = interpolate(t, x(:,1), duration, parameters);

try w = waitbar(6/8, w); catch, delete(w); error('Manually stopped'); end

% Generate single spike template
[templates.t, templates.d] = gen_template(interp_t, interp_d, 1, templates);
% Change to column vector and normalize to 0 and normalize max
templates.d = (templates.d + abs(templates.d(1)));
templates.d = templates.d / max(templates.d);

try w = waitbar(7/8, w); catch, delete(w); error('Manually stopped'); end

% Adjust the transition template lengths
templates = adj_templates(templates);

% Insert the sampling rate
templates.sampling_rate = parameters.sampling_rate;

try w = waitbar(8/8, w); catch, delete(w); error('Manually stopped'); end

% Close progress bar
try delete(w); catch E, fprintf(2,'\t%s\n',E.message); end

end

%% Interpolate data function

% This function will interpolate the data to the correct sampling rate
function [int_t, int_d] = interpolate(t, d, duration, parameters)

int_t = 0 : 1/parameters.sampling_rate*1000 : duration;
int_d = interp1(t, d, int_t, 'spline');

end

%% Calculate the absolute refractory period

function [refract_time, refract_index, rel_refract_time, rel_refract_index] = refract_period(tInit, xInit, duration, templates, parameters, w)

max_sampling_rate = 100000;
buffer = 0.0001;

refract = 0;
i = templates.initial_ap; % 1+ to the initial_ap that was defined

% Will loop until refract variable becomes true
while refract == 0
    i = i + 1/max_sampling_rate*1000; % This time should be 1/sampling_rate
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i) * 2).^2);
    [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
%     [int_t, int_d] = interpolate(t, x(:,1), duration, templates.sampling_rate);
    
    % Use find peaks to find the first 2 peaks which will be the 2 APs
    pks = findpeaks(x(:,1));
    pks = pks(pks >= 0);
    if length(pks) >= 2 && diff(abs(pks([1,2]))) < 5    % If more than 2 peaks above 0 and the diff between the two is less than 5
        refract = 1;    % Change refract to be true to break while loop
        max_sr_refract_time = i - templates.initial_ap;
        refract_index = find(0:1/parameters.sampling_rate*1000:duration >= i - templates.initial_ap, 1, 'first');  % This refractory time is the time between the TWO INPUT SPIKES. AT this time is when the first following spike can occur
        refract_time = 1/parameters.sampling_rate*1000 * (refract_index - 1);   % This refractory index is the values within the time matrix between the TWO INPUT SPIKES
    end
    if i > 60
        error('Unable to find refractory time due to lack of second spike');
        refract_time = NaN;
        refract_index = NaN;
    end
end

try w = waitbar(1/8, w); catch, delete(w); error('Manually stopped'); end

% Finding the relative refractory period at max sampling rate to keep
% consistency
count = 0;
rel_temp = [];

for i = max_sr_refract_time + templates.initial_ap:1/max_sampling_rate * 1000:duration
    count = count + 1;
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i)*2).^2);
    [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    d = x(:,1);
    rel_d = {d - d(1)};
    rel_temp = [rel_temp rel_d];
    
    if (i > max_sr_refract_time + templates.initial_ap + 5)...
            && ((max(rel_temp{:, count}) <= max(rel_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(rel_temp{:, count}) >= max(rel_temp{:, count - max_sampling_rate/1000}) - buffer))...
            && ((max(rel_temp{:, count}) <= max(rel_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(rel_temp{:, count}) >= max(rel_temp{:, count - (max_sampling_rate/1000) * 2}) - buffer))...
            && ((max(rel_temp{:, count}) <= max(rel_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(rel_temp{:, count}) >= max(rel_temp{:, count - (max_sampling_rate/1000) * 3}) - buffer))...
            && ((max(rel_temp{:, count}) <= max(rel_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(rel_temp{:, count}) >= max(rel_temp{:, count - (max_sampling_rate/1000) * 4}) - buffer))...
            && ((max(rel_temp{:, count}) <= max(rel_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(rel_temp{:, count}) >= max(rel_temp{:, count - (max_sampling_rate/1000) * 5}) - buffer))
        rel_refract_index = find(0:1/parameters.sampling_rate * 1000:duration >= (i - 5 - templates.initial_ap), 1, 'first');   % Takes first value that is larger than the one found with max_sampling_rate. This value is FROM initial spike at index = 151
        rel_refract_time = 1/parameters.sampling_rate * 1000 * (rel_refract_index - 1);    % This value is FROM the initial spike at t = 30
        break
    end
end
end

%% Generating the transition templates and calculating the end times/index

function trans_temp = gen_trans(tInit, xInit, duration, templates, parameters)

trans_temp = {};

% Loops through the entirety of the duration if required in steps of 0.2ms
for i = templates.abs_refract_time + templates.initial_ap:1/parameters.sampling_rate * 1000:templates.rel_refract_time+templates.initial_ap
    
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i)*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    [int_t, int_d] = interpolate(t, x(:,1), duration, parameters);
    [~, temp_d] = gen_template(int_t, int_d, 2, templates);
    
    % Places data into cell arrays
    temp_d = {temp_d - temp_d(1,1)};
    trans_temp = [trans_temp temp_d];
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

%start_ap = find(0:1/templates.sampling_rate*1000:200 == templates.initial_ap) - 5; % Start AP index. This index's as 7 in the interpolated time data

% Will adjust end_ap depending on the number of inputs
if num_inputs == 1
    end_ap = templates.rel_refract_index + templates.initial_ap_index;
    %end_ap = find(slope >= 0.05 | slope <= -0.05, 1, 'last');
else
    end_ap = length(t); % If there are 2 inputs, then it will adjust template later
    
end

% Adjusts the template
new_t = t(templates.initial_ap_index - 10: end_ap)';
new_d = (d(templates.initial_ap_index - 10: end_ap))';

end

%% Adjust templates.transition templates

% Currently only works if there are MAXIMUM of 2 inputs
function templates = adj_templates(templates)

% still use templates.t, maybe just adjust this variable to be outside the
% struct
templates.initial_ap_index = find(templates.t <= templates.initial_ap, 1, 'last');  % - 1 because of weird matlab indexing

% This is the number of positions between input spike and end index which
% will correspond to second input spike and end index for the second input
% spike
% Assumes that the following spike dynamics are same as the first
% index_bw_ap = templates.rel_refract_index - start_ap;

% Find the maximum value of the templates.transition cell array for
% normalisation
max_val = max(cellfun(@(x) max(x), templates.transition));

% Runs through the transition templates and shortens them to the correct
% length
for i = 1:size(templates.transition,2)
    templates.transition{i} = templates.transition{i}(1:templates.initial_ap_index ...
        + templates.abs_refract_index + i + templates.rel_refract_index);
    
    % Normalises the templates amplitude so that the max peak is = 1
    temp = templates.transition{i};
    temp = temp./max_val;
    templates.transition{i} = temp;
end
end

%% Store the parameters

function parametersUsed = store(parameters)

parametersUsed{1} = {['Sampling rate = ' num2str(parameters.sampling_rate)]};
parametersUsed{2} = {['eNa = ' num2str(parameters.eNa - parameters.vRest)] ['gNa = ' num2str(parameters.gNa)]};
parametersUsed{3} = {['eK = ' num2str(parameters.eK - parameters.vRest)] ['gK = ' num2str(parameters.gK)]};
parametersUsed{4} = {['eLeak = ' num2str(parameters.eLeak - parameters.vRest)] ['gLeak = ' num2str(parameters.gLeak)]};
parametersUsed{5} = {['C = ' num2str(parameters.C)]};
parametersUsed{6} = {['vRest = ' num2str(parameters.vRest)]};

end
