% GEN_TRAIN Simulates a spike train from extracellular recording
%
% Syntax:
%     v = gen_train(templates, Naxons, fs, duration, <options>);
%     [v, vv] = gen_train(templates, Naxons, fs, duration, <options>);
%     v = gen_train(templates, Naxons, fs, duration, 'SpikeRate', sr)
%     v = gen_train(templates, Naxons, fs, duration, 'Overlap', true);
%
% Inputs:
%     templates   -  matrix of spike templates where each column represents
%                    a spike template.
%     Naxons      -  Scalar indicating the number of different axons. Each
%                    axon is assigned a template randomly and a random
%                    average amplitude.
%     fs          -  Sampling rate
%     duration    -  Duration of the simulation in samples.
%
%     <options>
%     'SpikeRate' -  (Optional) Average spike rate of all axons if scalar,
%                    spike rate of each axon if a vector size Naxons x 1.
%                    Default value is 100 spikes per second for all axons.
%     'Overlap'   -  (Optional) If true (default), spikes can occur
%                    simultaneously.
%     'Recruited' -  (Optional) Number  of axons that are recruited along
%                    the recording. All the others start at time = 0.
%                    Default is 0.
%     'Dismissed' -  (Optional) Number of axons that are lost or stop
%                    spiking before the end of the recording. All others
%                    continue until time = end. Default is 0.
%     'Events'    -  (Optional, Struct) Contains fields that reflect in
%                    changes of spike rate or amplitude of some or all
%                    axons
%
% Outputs:
%     v           -  Simulation of an extracellular spiking recording.
%     vv          -  Simulation before adding the axons. Matrix of
%                    dimension Naxons x duration in samples.
%
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

%%
function [v, vv, report] = gen_train(templates, Naxons, fs, duration, varargin)
   % Default inputs
   opts.SpikeRate = 100;
   opts.Overlap   = true;
   opts.Recruited = 0;
   opts.Dismissed = 0;
   
   % Options
	% Validate inputs and assign optional values
   if nargin > 4 && mod(nargin,2)
      fprintf(2,'\tOptions must be pairs of Name(string) and Value\n');
   else
      for i = 1:2:nargin-4
         opts.(varargin{i}) = varargin{i + 1};
      end
   end
   
   % If SpikeRate is scalar, produce a vector of Naxons x 1 with random
   % numbers in the range [0.5*SpikeRate , 1.5*SpikeRate]
   if numel(opts.SpikeRate) == 1
      r = 0.5 + (1.5 - 0.5) .* rand(Naxons,1);
      opts.SpikeRate = r .* opts.SpikeRate;
   end
   
   % Amplitudes of each axon. Amplitude variability sampled from the range
   % [sqrt(2)/2, sqrt(2)]. Range of amplitudes taken from:
   % Rossant et al, 'Spike sorting for large, dense electrode arrays',
   % Nature Neuroscience, 2016 
   % sqrt(2)/2 + (sqrt(2) - sqrt(2)/2)
   r = rand(Naxons,1) .* (1 + poissrnd(1, [Naxons, 1])); % Random from Poisson distribution without 0
   amplitudes = ones(size(opts.SpikeRate));
   amplitudes = r' .* amplitudes;
   amplitudes(amplitudes < 0.5) = 0.5;
   
   % Starting time of each axon
   st_time = ones(size(opts.SpikeRate)); % Sample 1
   axs = randperm(Naxons, opts.Recruited); % Axons to be recruited along the recording
   st_time(axs) = randi(round(2*duration/3), [opts.Recruited, 1]);
   
   % Ending time of each axon
   end_time = duration * ones(size(opts.SpikeRate)); % Sample 1
   axs = randperm(Naxons, opts.Dismissed); % Axons to be dismissed along the recording
   end_time(axs) = randi(round([(duration/3) ,duration]), [opts.Dismissed, 1]);
      
   % Run the simulation
   [vv, rr] = run_simulation(Naxons, templates, fs, duration ,opts ,amplitudes, st_time, end_time);
   
   % Sum all axons
   v = sum(vv,2);
   
   % Normalize to maximum value
   v = v / max(v);
   
   % Return information about the simulation
   report = struct;
   report.opts    = opts;
   report.recruit = st_time;
   report.dismiss = end_time;
   % Copy the report from 'rr' ('run_simulation' function) to 'report'
   for fn = fieldnames(rr)'
      report.(fn{1}) = rr.(fn{1});
   end   
end

% Run the simulation. Returns spike trains per axon in vv.
function [vv, report] = run_simulation(Naxons, templates, fs, duration ,opts ,amplitudes ,st_time , end_time)
   dt            = 1/fs;
   T             = dt:dt:duration*dt;
   vv            = zeros(duration, Naxons);
   spks          = zeros(duration, Naxons);
   locs          = cell(Naxons, 1);
   max_spike_num = ceil(max(opts.SpikeRate(:)*duration*dt)); % Maximum number of spikes
%    rest          = round(100e-3/dt);% Refractory period in seconds
   
   % Output variable
   report = struct;
   % Check for special events
   if isfield(opts,'Events')
      % Inflammation
      if opts.Events.inflammation_axons > 0
         inflamed = randperm(Naxons, opts.Events.inflammation_axons);
         inf_time = opts.Events.inflammation_onset;
         tau = opts.Events.inflammation_tau;
         report.inf_time = inf_time;
         report.inflamed = inflamed;
      else
         inflamed = 0;
         report.inflamed = inflamed;
      end
      
      % Sudden change of amplitude for some axons (natural)
      if opts.Events.amplitude_nat_axons > 0
         amped = randperm(Naxons, opts.Events.amplitude_nat_axons);
         amp_time = opts.Events.amplitude_nat_onset;
         report.amped = amped;
         report.amp_time = amp_time;
      else
         amped = 0;
         report.amped = amped;
      end
      
      % Sudden change of all amplitudes (disturbance)
      amp_disturbance_onset = opts.Events.amplitude_dist_onset;
      amp_disturbance_value = opts.Events.amplitude_dist_value;
      amp_disturbance_probability = opts.Events.amplitude_dist_prob;
   end
   
   if opts.RepeatTemplates
      % Random list of templates (can repeat)
      templates_ = randi(length(templates), 1, Naxons);
   else
      % Random list of templates (non repeating)
      templates_ = randperm(length(templates), Naxons);
   end
   % Progress bar
   w = waitbar(0, 'Generating simulation...');
   for i = 1:Naxons
      % Check that start time is smaller than end time
      if end_time(i) < st_time(i)
         % if end is smaller, invert the values
         et = end_time(i);
         end_time(i) = st_time(i);
         st_time(i) = et;
      end
     
      
      currentTemplate = templates_(i); % Randomly pick 1 of the templates to assign to this axon.
      duration_of_spike = templates(currentTemplate).refract_index + length(templates(currentTemplate).transition); % Duration(idx) until first spike no longer transitions
      isi = random('Exponential', fs/opts.SpikeRate(i), [3*max_spike_num 1]);
      isi = round(isi);
      % Remove isi that are closer than the duration of a spike or
      % refractory period
      isi(isi <= templates(currentTemplate).refract_index) = []; % isi(isi < (size(templates,1) + rest)) = ceil(size(templates,1) + rest);
      
      
      % If it doesn't get affected by inflammation, it's firing rate
      % remains constant. If it does, we will remove the isi's who's cumsum
      % go beyond the inflamation time.
      % Check if this axon gets an inflammation event
      if ~isempty(find(inflamed == i,1))
         % If it gets inflamed, change the firing rate at 
         % time = evnts.inflammation_onset
         inf_sample = find(cumsum(isi) > inf_time, 1, 'first');
         sr = 9 * exp(-(1:numel(isi)-inf_sample)/(tau)) + 1; % Get an exponential from 10 to 1 with time constant tau
         isi(inf_sample + 1:end) = isi(inf_sample + 1:end)./sr';
         % Remove isi that are closer than the duration of a spike
         isi(isi <= templates(currentTemplate).refract_index) = [];
         isi = ceil(isi);
      end
      
      % Calculate times of all spikes
      sptimes = cumsum(isi);
      
      % Remove spikes that exceed the duration of the recording or the end
      % time of the particular axon
      sptimes(sptimes > end_time(i)) = [];
      % Remove spikes that occur before the axon is recruited
      sptimes(sptimes < st_time(i)) = [];
      
      % Check for overlapping spikes on different recordings. If overlap is
      % false, then don't allow overlapping.
      if ~opts.Overlap
      % Only if opts.overlap is false
         if i > 1
            for ii = 1:length(allsptimes)
               idx = (sptimes >= allsptimes(ii) - duration_of_spike)...
                      & (sptimes <  allsptimes(ii) + duration_of_spike);
               % Remove the opts.overlapped
               sptimes(idx) = [];
               % Update progress
               try w = waitbar((i-1)/Naxons + ii/(length(allsptimes)*Naxons), w); catch E, delete(w); error('Manually stopped'); end
            end
         else
            allsptimes = [];
         end
         allsptimes = [allsptimes; sptimes];
      end
      
      % Seperates the transition and non_transitioning spikes
      [non_transition, transition, transition_cells] = separate_transition_spikes(sptimes, duration_of_spike);
      
      
      % Create a recording of zeroes
      v_non_transition = zeros(duration, 1);
      v_transition = zeros(duration, 1);
      spks(:,i) = zeros(duration, 1); % Used for reports and not template generation
      
      % Assign binary spikes to the vector, the amplitude of the spikes is
      % weighted, instead of being just 1 or 0.
      v_non_transition(non_transition) = amplitudes(i);
      spks(sptimes,i) = 1; % Used for reports and not template generation
      % Vary the amplitude
      rand_amp = 0.99 + (1.01 - 0.99) .* rand(size(non_transition)); % small variation in amplitude
      v_non_transition(non_transition) = v_non_transition(non_transition).*rand_amp;
      
      % If the amplitude of current axon changes suddenly, scale all the
      % spikes after such time.
      % Ref: quirk2001, tsubokawa1996
      log_amp = 0;
      if ~isempty(find(amped == i,1))
         end_amp = (0.5 * rand(1,1) - 0.25); % Change in amplitude limited to 0.15 and -0.15
         log_amp =  (end_amp./(1 + exp(-10 * dt * ([1 : length(v_non_transition)] - amp_time)))); % logistic function
         amp = 1 + log_amp';
         v_non_transition = v_non_transition .* amp;
      end
      
      % Propagate the spike shape along the spikes vector
      v_non_transition = conv(v_non_transition,templates(currentTemplate).d);
      v_non_transition = v_non_transition(1 : duration, 1); % Remove trailing bits of convolution
      
      % Only run gen_transitions if there are transitions
      if ~isempty(transition)
          v_transition = (amplitudes(i) + log_amp') .* gen_transitions(transition_cells, templates(currentTemplate).transition, duration, templates(currentTemplate).refract_index);
      end
      
      % Assign the temporal variable v_ to the matrix of axons
      vv(:,i) = v_non_transition + v_transition;
      
      % Save the sike locations (times) in report.locs
      locs{i} = sptimes;
      
      % Update progress
      try w = waitbar(i/Naxons, w); catch, delete(w); error('Manually stopped'); end
   end
   
   % If the amplitude of all spikes change due to a disturbance, simply
   % multiply all the values of v_ after the change in time.
   if (amp_disturbance_onset > 0) && (rand < amp_disturbance_probability)
      vv(amp_disturbance_onset : end, : ) = vv(amp_disturbance_onset : end, : ) * amp_disturbance_value;
   else
       amp_disturbance_onset = NaN;
       report.opts.Events.amplitude_dist_onset = NaN; % If disturbance does not occur, set the disturbance onset time to NaN
   end
   
   report.locs = locs;
   report.spks = spks;
   report.templatesUsed = templates(templates_); % Record the templates used
   
   if ~exist('amp_time', 'var'); amp_time = NaN; end
   tempFamGroupings = categoriseTempFamGroups(templates_, Naxons, st_time, end_time, amp_disturbance_value, amp_disturbance_onset, amped, amp_time, locs);
   report.tempFamGroupings = tempFamGroupings';
      
   % Close progress bar
   try delete(w); catch E, fprintf(2,'\t%s\n',E.message); end

end

% Seperates the spikes which are close together and have transitions between them and
% those far apart without transitions.
function [non_transition, transition, transition_cells] = separate_transition_spikes(sptimes, duration_of_spike)
% Creates new isi since sptimes has been modified
isi = [sptimes(1); diff(sptimes)];

% Determine the times for spikes that do not transition
non_transition    = sptimes(1); % Allocates first spike into non_transition first regardless of whether its actually correct (this is fixed later down in the code)
for i = 1 : length(isi)-1
   if isi(i) <= duration_of_spike  || isi(i+1) <= duration_of_spike
       % When the isi for the next 2 spikes is <= rest, adds the 
       % current isi to the current cumulative sum
       non_transition(end,1) = sptimes(i+1); 
   elseif isi(i) > duration_of_spike
       % When the isi for the spikes is > rest, it creates a new index to
       % indicate the next normal (without transition) spike
       non_transition(end+1,1) = sptimes(i+1);
   else
       error('Error: There is a missing value in sptimes_normal');
   end
   
end
% Start and end exceptions that the for loop cant handle properly
% Since the first spike doesnt have a spike before it, it can be
% non-transitioning when the isi after it is > rest 
if ( isi(2) > duration_of_spike ) && ( non_transition(1) ~= sptimes(1) )
    non_transition = [sptimes(1); non_transition]; 
end
% Clears the last sptimes_normal if the last isi is less than the rest time
% since the for loop above can't check the preallocated last value (may
% need to change this when considering total duration of simulation)
if isi(end) <= duration_of_spike 
    non_transition(end) = []; 
end

% Determines the times of spikes with transitions (mutually exclusive to
% times without transitions)
transition = sptimes(~ismember(sptimes, non_transition));

% Find groups of transitioning spikes and split into cells
group_start_n_end = [0; find(diff(transition) > duration_of_spike); length(transition)];
for j = 1 : length(group_start_n_end) - 1
  transition_cells{j} = transition(group_start_n_end(j)+1 : group_start_n_end(j+1));
end

end

function tempFamGroupings = categoriseTempFamGroups(templates_, Naxons, st_time, end_time, amp_disturbance_value, amp_disturbance_onset, amped, amp_time, locs)

% Categorizes axons into families (same family if they have the same
% template)
temp_fam_groupings(:,2) = templates_';%temp_fam_groupings = zeros(size(templates_', 1), 2);
temp_fam_groupings(:,1) = (1:Naxons)';
%    n = 1;
%    for i = unique(templates_)
%        temp_fam_groupings(templates_'==i, 2) = n; % Changes the template num to start from 1,2,3... to follow how the extractor orders it
%        n=n+1;
%    end
temp_fam_groupings = num2cell(temp_fam_groupings); % Changes from matrix to cell

% add the onset time of External Disturbance since a change in amplitude greater
% than 10% will be grouped to a different family NEEDS TO BE CONFIRMED
change_amp_diff_fam = 0.1; % ADJUST THIS VALUE TO MAKE IT MORE OR LESS SENSITIVE TO CHANGING FAMILIES
family_st_end = num2cell([st_time' end_time'], 2); % add the start and end times of each axon
if ~isnan(amp_disturbance_onset) && (amp_disturbance_value > 1+change_amp_diff_fam || amp_disturbance_value < 1-change_amp_diff_fam)
    for ii = 1 : Naxons
        if amp_disturbance_onset > st_time(ii) && amp_disturbance_onset < end_time(ii) % Checks if the disturbance occurs while the axon is active
            family_st_end{ii} = [family_st_end{ii}(1), amp_disturbance_onset, family_st_end{ii}(2)];
        end
    end
end

% add the onset time for Natural Disturbance
if ~isequal(amped, 0) && (end_amp > change_amp_diff_fam || end_amp < -change_amp_diff_fam)
    for ii = amped
        if ~isnan(amp_time) && amp_time > st_time(ii) && amp_time < end_time(ii) % Checks if the disturbance occurs while the axon is active
            family_st_end{ii}(end+1) = amp_time; % Adds the natural disturbance onset to the end, this will be fixed in the next line
            family_st_end{ii} = sort(family_st_end{ii}, 'ascend'); % Sorts it so that family times will occur in the right order
        end
    end
end

% Allocates the sptimes of each naxon to a family and the family num
for ii = 1 : Naxons
    for iii = 1 : length(family_st_end{ii})-1
        temp_fam_groupings{ii, 3}{iii,1} = locs{ii}(and(locs{ii} >= family_st_end{ii}(iii),locs{ii} < family_st_end{ii}(iii+1)));
        temp_fam_groupings{ii, 3} = temp_fam_groupings{ii, 3}(~cellfun('isempty', temp_fam_groupings{ii, 3})); % Deletes family group with no actual sp in it (rare occurs but does)
        
        % Family num
        family_offset = 0;
        if ii > 1 % pass the first axon
            repeat_shape_idx = find(temp_fam_groupings{ii,2} == [temp_fam_groupings{1:ii-1, 2}]);
            if any(repeat_shape_idx)
                for n = (repeat_shape_idx)
                    family_offset = length(temp_fam_groupings{n,3}) + family_offset;
                end
            end
        end
        temp_fam_groupings{ii, 4} = (1:(length(family_st_end{ii})-1)) + family_offset;
    end
end

% Assigns the family and template groupings to a field in the report
tempFamGroupings = cell2struct(temp_fam_groupings, {'Naxon', 'template', 'family_sptimes', 'family_num'}, 2);
end