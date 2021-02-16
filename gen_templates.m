% gen_templates script which is to be used in tandem with the SpikeSimulationTool to generate axon templates in accordance to the
% Hodgkin Huxley neural model

% Note: HH = Hodgkin Huxley

% University of Melbourne
% Department of Biomedical Engineering
% Summer Research Experience Program - Spike Simulation Suite
% Christopher Ong & Chi Yung Darren Tan | 13 February 2021

%% gen_template function

function [templates] = gen_templates(parameters)

duration = 200;                         % Duration of the simulation [msec]
tInit    = [0 duration];                % Start and end time
xInit    = [-65; 0.052; 0.059; 0.317];  % Initial HH parameters
templates.initial_ap = 30;              % Fixed time for initial current input
templates.initial_ap_index = find(0:1/parameters.sampling_rate*1000:duration >= templates.initial_ap, 1, 'first');  % Initial current input index

% Progress bar
w = waitbar(0, 'Generating templates...');

% Store parameters into the templates struct for reference
templates.HHparameters = parameters;

% Calculate the refractory period
[templates.refract_time, templates.refract_index, templates.trans_time, templates.trans_index] = refract_period(tInit, xInit, duration, templates, parameters, w);

% Update progress bar
try w = waitbar(2/8, w); catch, delete(w); error('Manually stopped'); end

% Calculate when second input will generate same output as previous output
% and generate transition templates
templates.transition = gen_trans(tInit, xInit, duration, templates, parameters);

% Update progress bar
try w = waitbar(3/8, w); catch, delete(w); error('Manually stopped'); end

% Input function
Iapp = @(t) 10*exp(-((t - templates.initial_ap)*2).^2);

% Update progress bar
try w = waitbar(4/8, w); catch, delete(w); error('Manually stopped'); end

% Runs ODE
[t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);

% Update progress bar
try w = waitbar(5/8, w); catch, delete(w); error('Manually stopped'); end

% Interpolate data
[interp_t, interp_d] = interpolate(t, x(:,1), duration, parameters);

% Update progress bar
try w = waitbar(6/8, w); catch, delete(w); error('Manually stopped'); end

% Generate single spike template
[templates.t, templates.d] = adj_template_index(interp_t, interp_d, 1, templates);

% Update progress bar
try w = waitbar(7/8, w); catch, delete(w); error('Manually stopped'); end

% Adjust the transition template lengths
templates = adj_trans_templates(templates);

% 1st derivative of all templates to change intracellular spike to
% extracellular spike shape
templates = intra2extra(templates, parameters);

% Update progress bar
try w = waitbar(8/8, w); catch, delete(w); error('Manually stopped'); end

% Close progress bar
try delete(w); catch E, fprintf(2,'\t%s\n',E.message); end

end

%% Interpolate data function

% This function will interpolate the data to the correct and a constant sampling rate
function [int_t, int_d] = interpolate(t, d, duration, parameters)

int_t = 0 : 1/parameters.sampling_rate*1000 : duration; % Interpolated time matrix
int_d = interp1(t, d, int_t, 'spline');                 % Interpolated data matrix

end

%% Calculate the absolute refractory period

% Function will calculate the refractory period times
% To ensure for consistency, these values are calculated with reference to
% the maximum sampling rate of 100kHz which will be set as the maximum
% sampling rate that can be inputted in parameters struct.
function [refract_time, refract_index, trans_time, trans_index] = refract_period(tInit, xInit, duration, templates, parameters, w)

max_sampling_rate = 100000; % Maximum sampling rate that will ever be used
buffer = 0.0001;            % Error buffer when calculating transition templates

refract = 0;                % Conditional variable
i = templates.initial_ap;

% Will loop until refract variable becomes true
while refract == 0
    i = i + 1/max_sampling_rate*1000;   % Iterating variable
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i) * 2).^2);   % Current input for HH
    [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);                  % HH ode
    
    % Use find peaks function to find the time when two peaks will be
    % generated above 0. This will represent that two action potentials
    % have been fired
    pks = findpeaks(x(:,1));
    pks = pks(pks >= 0);
    if length(pks) >= 2 && diff(abs(pks([1,2]))) < 5                            % If more than 2 peaks above 0 and the diff between the peak value is less than 5
        refract = 1;                                                            % Change refract to be true to break while loop
        max_sr_refract_time = i - templates.initial_ap;                         % Maximum sampling rate refract time
        refract_index = find(0:1/parameters.sampling_rate*1000:duration >= i - templates.initial_ap, 1, 'first');  % This refractory time is the time between the TWO INPUT SPIKES. AT this time is when the first following spike can occur
        refract_time = 1/parameters.sampling_rate*1000 * (refract_index - 1);   % This refractory index is the values within the time matrix between the TWO INPUT SPIKES
    end
    if i > 60   % When iterating variable exceeds this value, second spike will not occur due to non-compliance with HH model
        error('Unable to find refractory time due to lack of second spike');
        refract_time = NaN;
        refract_index = NaN;
    end
end

% Update progress bar
try w = waitbar(1/8, w); catch, delete(w); error('Manually stopped'); end

% Finding the relative refractory period at max sampling rate to keep
% consistency
count = 0;  % Count variable for template cells
trans_temp = [];

for i = max_sr_refract_time + templates.initial_ap:1/max_sampling_rate * 1000:duration
    count = count + 1;
    Iapp = @(t) 10 * exp(-((t - templates.initial_ap) * 2).^2) + 10 * exp(-((t - i)*2).^2); % Current input for HH
    [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);              % HH ode
    trans_d = x(:,1);                                 % Voltage output 
    trans_temp = [trans_temp {trans_d}];          % Storing templates in an matrix
    
    % Tests for when 5 consecutive integers for i produces a max value that
    % is within the buffer range. When there are 5 consecutive values that
    % are within this range, it is considered that any further templates
    % will generate the same result, hence loop will be stopped for
    % efficiency purposes
    if (i > max_sr_refract_time + templates.initial_ap + 5)...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - max_sampling_rate/1000}) - buffer))...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - (max_sampling_rate/1000) * 2}) - buffer))...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - (max_sampling_rate/1000) * 3}) - buffer))...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - (max_sampling_rate/1000) * 4}) - buffer))...
            && ((max(trans_temp{:, count}) <= max(trans_temp{:, count - max_sampling_rate/1000}) + buffer) && (max(trans_temp{:, count}) >= max(trans_temp{:, count - (max_sampling_rate/1000) * 5}) - buffer))
        trans_index = find(0:1/parameters.sampling_rate * 1000:duration >= (i - 5 - templates.initial_ap), 1, 'first');   % Takes the first index value which corresponds to 1 from the 5 consecutive integers
        trans_time = 1/parameters.sampling_rate * 1000 * (trans_index - 1);    % Respective time value to the index
        break
    end
end
end

%% Generating the transition templates

% Function will generate transition templates in accordance to the
% index/time that was calculated in refract_period function
function trans_temp = gen_trans(tInit, xInit, duration, templates, parameters)

trans_temp = {};    % Store transition templates in cell array

% Loops through the entirety of the duration if required in steps of
% 1/sampling_rate
for i = templates.refract_time + templates.initial_ap:1/parameters.sampling_rate * 1000:templates.trans_time+templates.initial_ap
    
    Iapp = @(t) 10 * exp(-((t - 30) * 2).^2) + 10 * exp(-((t - i)*2).^2);       % Current input for HH
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);  % HH ode
    [int_t, int_d] = interpolate(t, x(:,1), duration, parameters);              % Interpolate data to inputted sampling rate
    [~, temp_d] = adj_template_index(int_t, int_d, 2, templates);               % Correct the starting index of templates
    
    trans_temp = [trans_temp {temp_d}];   % Store transition template as cell arrays
end
end

%% Template generator function

% This function will adjust the start and end of templates.
% For the master templates (single spike), this function will taper both
% the start and end to the correct length.
% For transition templates, this function will taper the start of the
% template. End will be adjusted later
function [new_t, new_d] = adj_template_index(t, d, num_inputs, templates)

% Will adjust end_ap depending on the number of inputs
% Input == 1 is for master templates
% Anything else is for transition templates, however this code is
% specifically done for MAXIMUM of 2 inputs
if num_inputs == 1
    end_ap = templates.trans_index + templates.initial_ap_index;
else
    end_ap = length(t); % If there are 2 inputs, then it will adjust template later
end

% Adjusts the template
new_t = t(templates.initial_ap_index - 10: end_ap)';
new_d = (d(templates.initial_ap_index - 10: end_ap))';

end

%% Adjust templates.transition templates

% This function will adjust the transition templates. It will taper the end
% of the templates making use of the previously calculated values
% (Currently only works if there are MAXIMUM of 2 inputs)
function templates = adj_trans_templates(templates)

% (Still need to try find a way to not use .t)
templates.initial_ap_index = find(templates.t <= templates.initial_ap, 1, 'last');  % Adjusts initial AP index to the new start index from adj_template_index function

max_val = max(cellfun(@(x) max(x), templates.transition));  % Find the maximum value of the transition cell array for normalisation

% Iterates through the transition templates and tapers the end of the
% template
for i = 1:size(templates.transition,2)
    templates.transition{i} = templates.transition{i}(1:templates.initial_ap_index ...
        + templates.refract_index + i + templates.trans_index);
end
end

%% Converts the intracellular spike templates into extracellular spike templates

% This function converts the templates from intracellular to extracellular.
% This is done by taking the 1st derivative of the intracellular template

% It will also normalise the extracellular template to work with the SpikeSimulationTool
% and the SpikeExtractionTool. Normalise to the minimum point of the extracellular template

function templates = intra2extra(templates, parameters)

% Master template
templates.d = (diff(templates.d)./diff(templates.t));   % 1st derivative to convert from intracellular to extracellular
templates.d = -templates.d./min(templates.d);           % Normalising template to the minimum point of extracellular template
templates.t = templates.t(1:end-1);                     % Adjusting time matrix to be same length as data matrix

% Transition templates
for i = 1:size(templates.transition,2)
    templates.transition{i} = (diff(templates.transition{i})./(1/parameters.sampling_rate));    % 1st derivative to convert from intracellular to extracellular
    templates.transition{i} = -templates.transition{i}./min(templates.transition{i});           % Normalising template to minimum point
end
end
