% SPIKESIMULATIONTOOL Script to simulate a spike train from single channel extracellular recording from various cells and drift*.
%
% *Drift is the growth of spike amplitude and noise level at a different
%  rate.
%
% University of Melbourne
% Department of Biomedical Engineering
% Artemio Soto-Breceda | 6/August/2019

%% Copyright 2019 Artemio Soto-Breceda
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%% Options
Naxons      = 5;                    % Number of different templates
SNR         = 20;                   % Initial signal to noise ratio (it will change with drift)
growth      = [(1.1 + (2 - 1.1) * rand) (0.8 + (1.2 - 0.8) * rand)]; % [1.9 1.1];         % Growth of: [<spamp> <noise>]
total_time  = 2000;                 % Seconds
fs          = 5000;                 % Sampling rate (current template file has this sampling rate, so it should stay like this unless the templates are fixed)
sr          = randi(10,1,Naxons)/2; % Spike rate
overlap     = false;                % If true, it allows spikes of diff axons to overlap
rpt_temp    = true;                 % True if more than one axon can have the same template
pw_locs     = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs); % Locations of the change of drift. Only controls the growth of the spike amplitudes, not noise
pw_grow     = 1.2 + (1.5 - 1.2) * rand; % The new drift for each location. Only controls the growth of the spike amplitudes, not noise
has_drift   = true; 
has_noise   = true;
pre_noise   = true;                 % Append period of just noise at the start
do_filter   = true;                 % Bandpass filter the signal
passband    = [40 1200];%[80, 600]; % Passband
PLOT        = false;

%% Events
evnts.inflammation_onset   = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs);  % High frequency at time
evnts.inflammation_tau     = 5e-3*fs;  % Time constant for increased spike rate to decay to spontaneous activity
evnts.inflammation_axons   = floor(0 + ((Naxons/2 - 0) * rand)); % Number of inflamed axons (increase the spike rate).

evnts.amplitude_nat_onset  = 500*fs;   % Change of amplitude in just some of axons
evnts.amplitude_nat_axons  = 0;        % Change of amplitude in just a couple of axons

evnts.amplitude_dist_onset = round((total_time/4 + ((total_time*3/4) - (total_time/4)) * rand) * fs);  % Change of amplitude in all axons
evnts.amplitude_dist_value = 0.5 + (1.5 - 0.5) * rand;      % Value of the new amplitude multiplier
evnts.amplitude_dist_prob  = 0.2; % Probability of having a change in the amplitude

evnts.prob_start  = floor(0 + ((Naxons/2 - 0) * rand)); % (Recruited) Number of axons that don't start at the beginning. They will randomly start somewhere along the recording.
evnts.prob_end    = floor(0 + ((Naxons/2 - 0) * rand)); % (Dismissed) Number of axons that don't last the whole recording. They will randomly end somewhere along the recording.

%% Run
% Load the templates matrix
load(['templates' filesep 'templates4']);

% Normalize templates amplitude, max peak = 1
for i =1:size(d,2)
   dd = d(:,i);
   max_d = max(d);
   dd = dd./max_d(i);
   d(:,i) = dd; %#ok<SAGROW>
end

dt = 1/fs;
try % Generate a train of extracellular spikes. There is no noise
   [v, vv, report] = gen_train(d, Naxons, fs, total_time/dt, 'SpikeRate',...
      sr, 'Overlap', overlap, 'Recruited', evnts.prob_start,...
      'Dismissed', evnts.prob_end, 'Events', evnts, 'RepeatTemplates',...
      rpt_temp);
catch E
   if strcmp('Manually stopped', E.message)
      fprintf(2,'\tManually stopped\n');
      return;
   else
      rethrow(E);
   end
end

% Add drift
if has_drift
   % Add the noise and drift
   v = add_drift(v, 'SNR', SNR, 'Noise', has_noise, 'Growth', growth,...
      'PrecedingNoise', pre_noise, 'Linear', false,...
      'PwLocs', pw_locs, 'PwGrowth', pw_grow);
elseif has_noise
   % Add only noise
   v = add_drift(v, 'SNR', SNR, 'Noise', has_noise, 'Growth', [1 1],...
      'PrecedingNoise', false);
end

if do_filter && has_noise
   % Lowpass filter the signal
   v = bandpass(v, passband, fs);
end

%% Format recording for SpikeExtractionTool
vsim = struct;
vsim.type = 'voltage';
vsim.name = 'vsim';
vsim.dt = dt;
vsim.params = [];
vsim.data = v;
vsim.time = (dt:dt:length(vsim.data)*dt);
vsim.axons = vv;
vsim.report = report;
% Fix the dimensions if they are wrong
if size(vsim.data,1) < size(vsim.data,2), vsim.data = vsim.data'; end
if size(vsim.time,1) < size(vsim.time,2), vsim.time = vsim.time'; end


if PLOT
%%
   figure;
   plot(vsim.time, vsim.data);
   xlabel('Time (s)');
   ylabel('Amplitude');
   
   figure;
   vv = vsim.axons;
   max_vv = max(vv);
   [~, srt] = sort(max_vv, 'descend');
   try
      if ~has_drift || ~pre_noise
         plot(vsim.time, vv(:,srt));
      else
         plot(vsim.time(101:end), vv(:,srt));
      end
   catch
      plot(vsim.time(101:end), vv(:,srt));
   end
   xlabel('Time (s)');
   ylabel('Amplitude');
   naxons = size(vsim.axons, 2);
   % Plot templates
   figure;
   for i = 1:naxons
      nr = ceil(sqrt(naxons));
      nc = ceil(naxons/3);
      if naxons == 3
         nc = 1; nr = 3;
      end
      subplot(nc, nr, i);
      plot(vv( max(vsim.report.locs{i}(1) - 100, 1) : min(vsim.report.locs{i}(1) + 100, size(vv,1)), i));
   end
   clear vv
end

% Print report
try
   if isfield('inf_time', report)
      fprintf('\tInflammation: %.02f s| Number of inflamed axons: %d\n', report.inf_time * dt, numel(report.inflamed));
   else
      fprintf('\tNumber of inflamed axons: %d\n', 0);
   end
   fprintf('\tAmplitude change time: %.02f s\n', report.opts.Events.amplitude_dist_onset * dt);
   fprintf('\tRecruited axons: %s\n', num2str(find(report.recruit > min(report.recruit))'));
   fprintf('\tDismissed axons: %s\n', num2str(find(report.dismiss < max(report.dismiss))'));
catch E
   fprintf('\tCouldn''t print the report. Unexpected event: %s\n', E.message);
end
%% Save
% Check size of the variable to save. If it is larger than 2Gb, remove the
% per axon field: vsim.axons
s = whos('vsim');
if s.bytes > 2e9
   vsim.axons = [];
   fprintf('\tThe file is too large, only the final recording will be saved, not the per-axon information.\n');
end

% [file,path] = uiputfile(['simulations' filesep 'sim.mat'],'Save file name');
% if file
%    file_name = [path filesep file];
%    save(file_name, 'vsim');
% else
%    fprintf('\tUser didn''t chose a file location. The simulation wasn''t saved.\n');
% end
sufix = 0;
file = 'sim';
valid_file = file;
while exist([valid_file,'.mat'], 'file')
   sufix = sufix + 1;
   valid_file = [file num2str(sufix)];
end
file = valid_file;
save(file, 'vsim');
