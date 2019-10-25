% ADD_DRIFT adds drift and noise to an extracellular recording
%
% Syntax: v = add_drift(recording, <options>);
%
% Inputs:
%     v_in     - Input recording
%     options  - [optional] String defined parameters:
%                 'Growth':   1x2 vector. Defines the total growth of spike
%                             amplitude and noise respectively. If the
%                             input is a 1x1 vector, both will grow
%                             equally. Default is [2 1.2], i.e. resulting
%                             spikes amplitude will be 100% larger at the
%                             end of the signal, and resulting noise will
%                             20% larger at the end.
%                 'SNR':      Signal to noise ratio. Default is 20
%                 'Noise':    true or false. If true, there is added noise.
%                             Default is true.
%                 'PrecedingNoise': If true, concatenates 100 samples of
%                             noise at the start. Default is true.
%                 'Th':       Only affect signal within the noise amplitude
%                             Default is false
%                 'Linear':   true or false. If true, the noise and spike
%                             drift is full linear, if false, the drift
%                             becomes piecewise linear. Default is true.
%                 'Pieces':   How mant times will the growth change if
%                             Linear is false, i.e. if piecewise linear.
%                             Default is 5 Pieces and there is a 25% 
%                             chance for each of the pieces to occur.
%                 'PwLocs':   Locations of the change of drift. If
%                             piecewise location drift is enabled, i.e.
%                             Linear is false, then the locations of the
%                             changes in drift are decided randomly.
%                             However, if PwLocs is given as an input,
%                             instead of randomly chosing where the drift
%                             changes, it is determined by PwLocs. It is a
%                             vector of times in samples.
%                 'PwGrowth': The new drift for each PwLocs. If PwLocs is
%                             input, there has to be PwGrowth.
%
% Outputs:
%
%     v        - The output recording with added drift.
%
% University of Melbourne
% Department of Biomedical Engineering
% Created by Artemio Soto-Breceda | 1/August/2019
% Last edit | 7/August/2019

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
function v = add_drift(v_in, varargin)
   % Default options
   opts.Linear = true;
   opts.Pieces = 5;
   opts.Growth = [2 1.2];
   opts.Noise = true;
   opts.SNR = 20;
   opts.PrecedingNoise = true;
   opts.Th = false;
   opts.PwLocs = 0;
   opts.PwGrowth = 0;
   
   % Validate inputs and assign optional values
   if ~mod(nargin,2)
      fprintf(2,'\tOptions must be pairs of Name(string) and Value\n');
   else
      for i = 1:2:nargin-1
         opts.(varargin{i}) = varargin{i + 1};
      end
      
      if length(opts.Growth) == 1
         opts.Growth = [opts.Growth(1) opts.Growth(1)];
      end
      
      if opts.PwLocs ~= 0 && opts.PwGrowth == 0
         error('If ''PwLocs'' is an input, ''PwGrowth'' must be provided as well.');
      end
      
      if opts.PwLocs > numel(v_in)
         error('''PwLocs'' can''t be larger than the duration of the recording.');
      end
   end
   
   v = get_drift(v_in, opts);
   % Intercede with zeros
   if opts.PrecedingNoise
      if isrow(v)
         v = [zeros(1,100) v];
      else
         v = [zeros(100,1); v];
      end
   end
   
   if opts.Noise
      % Calculate the variance of the signal
      [a, b] = findpeaks( v_in, 'MinPeakHeight', 0.01, 'MinPeakWidth', 3);
      rms_v = (sum(a.^2))/numel(a);
      % Change the variance of the noise to comply with the SNR
      % var_n = (rms_v^2)/opts.SNR;
      var_n = rms_v/(10^(opts.SNR/10));
      % Generate noise with variance var_n
      noise = get_noise(length(v), opts, var_n);           
      
      vv = v + noise; % Add it to signal
      if opts.Th
         idx = find(abs(vv) > 1.1*abs(noise));
         vv(idx) = v(idx);
      end
      v = vv;
   end
      
end

function v = get_drift(v, opts)
   g = opts.Growth(1);
   L = length(v);
   
   if size(v,1) > size(v,2), v = v'; end
   
   % If linear growth
   if opts.Linear
      drift_v = linspace(1, g, L); % Drift of signal
   else
      drift_v = linspace(1, g, L); % Drift of signal
      if opts.PwLocs == 0
         [piecewise_locations, piecewise_growth] = get_piecewise(opts.Pieces, L);

         for i = 1:numel(piecewise_locations)
            new_drift = drift_v(piecewise_locations(i) : end);
            drift_v(piecewise_locations(i) : end) = linspace(new_drift(1), new_drift(1) + piecewise_growth(i), length(new_drift));
         end
      else
         for i = 1:numel(opts.PwLocs)
            new_drift = drift_v(opts.PwLocs(i) : end);
            drift_v(opts.PwLocs(i) : end) = linspace(new_drift(1), opts.PwGrowth(i), length(new_drift));
         end
      end
   end
   v = v .* drift_v; % Drifted signal
end

function n = get_noise(L, opts, variance)
   g = opts.Growth(2);
   
   % If linear growth
   if opts.Linear
      drift_n = linspace(1, g, L); % Drift of noise
   else
      drift_n = linspace(1, g, L); % Drift of signal
      [piecewise_locations, piecewise_growth] = get_piecewise(opts.Pieces, L);
      
      for i = 1:numel(piecewise_locations)
         new_drift = drift_n(piecewise_locations(i) : end);
         drift_n(piecewise_locations(i) : end) = linspace(new_drift(1), new_drift(1) + piecewise_growth(i), length(new_drift));
      end
   end
   
   % Example: Generate values from a normal distribution with mean 1
   % and standard deviation 2.
   % r = 1 + 2.*randn(100,1);
   n = sqrt(variance).*randn(1,L); % Normal (white gaussian) distribution
   n = n .* drift_n; % Drifted noise
end

function [pl, pg] = get_piecewise(Npieces, L)
   pg = rand(Npieces, 1) - 0.5; % Random growth between -0.5 and 0.5
   pg(rand(Npieces, 1) < 0.75) = []; % 75% probability of each to be removed
   
   % Find random Piece Locations
   pl = randi(round(2*L/3), [numel(pg), 1]);
   % Sort the locations ascendently 
   pl = sort(pl);
end