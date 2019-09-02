% SIMULATION_TO_WAVE_CLUS Formats the SST simulation into Wave_Clus format
%
% 	Syntax:
%     success = simulation_to_wave_clus();
%     success = simulation_to_wave_clus(vsim);
%     success = simulation_to_wave_clus(vsim, vout);
%     %(Run as a script to chose from a browse window)
%
%  Input:
%     vsim  -  Optional (struct), this is the output of a SST simulation.
%              If not provided, then the function prompts a file selection
%              for the user to chose the .mat variable to load from the
%              hard drive.
%     vout  -  Optional (string), Path and Filename of the converted
%              simulation.
%
%  Output:
%     success - (boolean) False if something failed.
% 
% Artemio Soto-Breceda | 20/August/2019
function varargout = simulation_to_wave_clus(varargin)
   success = false;
   if nargin >= 1
      
      vsim = varargin{1};
      if nargin >= 2
         save_file = vout;
         varargout = {success};
      else
         % Chose a name and location for the new file
         [file,path] = uiputfile(['simulations' filesep 'sim_voltage.mat'],'Save file name');
         if file
            save_file = [path filesep file];
         else
            fprintf('\tUser didn''t chose a file location. The simulation wasn''t saved.\n');
            return;
         end
      end
   else
      
      [file,path] = uigetfile(['..' filesep 'simulations'],'Load simulation');      
      if file         
         file_name = [path filesep file];
         d = load(file_name);
         fn = fieldnames(d);
         vsim = d.(fn{1});
         % Set an aribtrary name for the new file
         file_name = file_name(1 : strfind(file_name, '.mat') - 1);
         save_file = [file_name, 'waveclus'];
      else
         fprintf('\tUser didn''t chose a file location. The simulation wasn''t loaded.\n');
         return;
      end      
      
   end
   try
      data = vsim.data; % data vector
   catch
      error('Wrong file chosen, try again');
   end
   sr = 1/vsim.dt;% Sampling frequency
   sort_fmax = min(3000, floor(sr/2)); % Maximum low-pass cut-off frequency (half of the sampling freq)
   
   % Save file
   save(save_file, 'data', 'sr', 'sort_fmax');
   % Success
   success = true;
   
   if nargin >= 1, varargout = {success}; end
end