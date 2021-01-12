function [new_template] = HHSim(duration, transitions)

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

tInit    = [0 duration];
xInit    = [-65; 0.052; 0.059; 0.317];

% Run ODE
[t, x] = ode45('HHode', tInit, xInit, [], transitions, const);

% Generate new template
[~, new_template] = gen_template(t, x(:,1), duration);
new_template = new_template + 65; % Adjust RMP to 0
new_template = new_template / max(abs(new_template)); % Normalize

end


%% Create new template function
function [new_t, new_template] = gen_template(t, data, duration)
fr = 5000; %sampling rate

% Interpolate data for constant sampling rate
new_t = 0 : 1/fr*1000 : duration;   % Time vector
new_template = interp1(t, data, new_t, 'spline');   % Value vector

end
