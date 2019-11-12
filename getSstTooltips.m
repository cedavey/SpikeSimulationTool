% GETTOOLTIPS Retrieves the list of all tooltips requested by objectName
%
% Syntax: str = getTooltips(objectName);
%
% Inputs:
%     objectName - Cell vector: String or strings containing the names of
%                  the objects for which to retrieve the tooltip.
%     
% Outputs:
%     str - Cell vector: String/s with the tooltip.
%
% Artemio - 12/November/2019
function str = getSstTooltips(objectName)
   % Create the list of all available tooltip strings
   tooltips = struct('nAxonsTextBox',        'Number of different axons to simulate',...
      'totalTimeTextBox',                    'Total duration of the recording in seconds',...
      'samplingRateTextBox',                 'Samples per second',...
      'maxSpikeRateTextBox',                 sprintf('Maximum possible spike rate. The spike rate is randomized\nfor each axon as a number between 1 and Max'),...
      'overlapCheckBox',                     'If ticked, more than one spike can occur at the same time. If untick, templates will never overlap',...
      'repeatTemplateCheckBox',              'If ticked, more than one spike can have the same template (spike shape)',...
      'snrTextBox',                          'Initial signal to noise ratio (SNR)',...
      'hasNoiseCheckBox',                    'If ticked, the recording will include noise',...
      'hasDriftCheckBox',                    'If ticked, the recording will have drift, i.e. the overall amplitude (and SNR) will vary with time',...
      'preNoiseCheckBox',                    'If ticked, there will be 100 samples of pure noise appended at the beginning of the signal',...
      'doFilterCheckBox',                    'If ticked, the signal will be filtered with a bandpass filter',...
      'lowBandTextBox',                      'Low band cut-off frequency for the bandpass filter',...
      'highBandTextBox',                     'High band cut-off frequency for the bandpass filter',...
      'plotCheckBox',                        'If ticked, the simulation will plot the results in new figures',...   
      'runButton',                           'Run Simulation');
   
   % Create output variable as an empty cell
   str = cell(1,length(objectName));
   
   for i = 1:length(objectName)
      try
         % Replace any space in the object name for an underscore
         objectName{i}(strfind(objectName{i}, ' ')) = '_';
         % Assign each tooltip string to the i-th cell. There will be i
         % output string cells for i input object names
         str{i} = tooltips.(objectName{i});
      catch
         % Most likely it didn't find the desired tooltip string
         % Return empty string.
         str{i} = '';
         % And print a message
         msgStr = cprintf('\tThe tooltip ''%s'' does not exist\n', objectName{i});
         cprintf('Text', msgStr);
      end
   end

end