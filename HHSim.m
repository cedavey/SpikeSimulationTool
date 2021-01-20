function [templates] = HHSim(parameters)

duration = 100; % [msec]
tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];

% Iapp function
Iapp = @Iapp_func; % applied current injection (can be written as function of t)

% Run ODE
[t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);

% Generate new template
[templates.t, templates.d] = gen_template(t, x(:,1), duration);

templates.refract_time = refract_period(tInit, xInit, parameters, duration);

templates.transition = gen_trans(tInit, xInit, parameters, templates, duration);

end

%% Iapp function
function Iapp_out = Iapp_func(t)

% Bell curve function
Iapp_out = 10 * exp(-((t-30)*2).^2);
           %10 * exp(-((t-46)*2).^2);
           %20 * exp(-((t-70)/2).^2); 

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

% Change to column vector and normalize to 0
new_t = new_t';
new_template = (new_template + abs(new_template(1)))';

end


%% Calculate refractory period function
function refract_time = refract_period(tInit, xInit, parameters, duration)

refract = 0;
initial_ap = 30;
i = initial_ap + 1; % This corresponds to +1 of the initial input spike at t = 15
while refract == 0
    Iapp = @(t) 10 * exp(-((t-initial_ap)*2).^2) + 10 * exp(-((t-i)*2).^2);
    [~, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    
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
function trans_temp = gen_trans(tInit, xInit, parameters, templates, duration)

% Initialising values
trans_temp = {};    % Used for if want to store data in cell arrays
initial_ap = 30;    % This is the first input spike

% Iterate through this for the number of ms of templates you want
% It is currently at 10 templates from right after the abs refractory
% period
for i = templates.refract_time:templates.refract_time + 9
    Iapp = @(t) 10 * exp(-((t-initial_ap) * 2).^2) + 10 * exp(-((t-(initial_ap + i))*2).^2);
    [t, x] = ode45('gen_templates_HHode', tInit, xInit, [], Iapp, parameters);
    
    [~, temp2] = gen_template(t, x(:,1),duration);
    
    % This uses cell arrays instead
    temp2 = {temp2' - temp2(1,1)};
    trans_temp = [trans_temp temp2];
    
end
end
