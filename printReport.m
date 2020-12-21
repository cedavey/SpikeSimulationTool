% 17th December 2020
% Chi Yung Darren Tan 910828

% This function replaces the print report within SpikeSimulationTool.m

% Function is called within SST.mIapp and will print the report within
% app.reportTextArea in the GUI

function report_op = printReport(report, dt)

try
   if isfield(report, 'inf_time')
      % Finding which fields have different values
      i = find(report.recruit > min(report.recruit));
      j = find(report.dismiss < max(report.dismiss));
      
      %Initialisation
      recruited_s = [];
      dismissed_s = [];
      
      % If find reports vector length 1, simply print that value
      if length(i) == 1
          recruited_s = num2str(i');
          
      % If length of find is larger than 1, iterate through its length to
      % print correct statement
      elseif length(i) > 1
          recruited_s = [num2str(i(1))];
          for x = 2:length(i)
              recruited_s = [recruited_s ', ' num2str(i(x))];
          end
      end
      % Same as above but for dismissed
      if length(j) == 1
          dismissed_s = num2str(j');
      elseif length(j) > 1
          dismissed_s = [num2str(j(1))];
          for x = 2:length(j)
              dismissed_s = [dismissed_s ', ' num2str(j(x))];
          end
      end
      report_op = sprintf('Start time of inflammation: %.02f s | Number of inflamed axons: %d\nTime when amplitude changes: %.02f s\nRecruited axon(s): %s\nDismissed axon(s): %s\n', ...
          report.inf_time * dt, numel(report.inflamed), report.opts.Events.amplitude_dist_onset * dt, recruited_s, dismissed_s);
   else
      report_op = sprintf('Number of inflamed axons: %d\nTime when amplitude changes: %.02f s\nRecruited axons: %s\nDismissed axons: %s\n', ...
          0, report.opts.Events.amplitude_dist_onset * dt, recruited_s, dismissed_s);
   end
catch E
   report_op = sprintf('Couldn''t print the report. Unexpected event: %s\n', E.message);
end
end


% try
%    if isfield(report, 'inf_time')
%       fprintf('\tInflammation: %.02f s| Number of inflamed axons: %d\n', report.inf_time * dt, numel(report.inflamed));
%    else
%       fprintf('\tNumber of inflamed axons: %d\n', 0);
%    end
%    fprintf('\tAmplitude change time: %.02f s\n', report.opts.Events.amplitude_dist_onset * dt);
%    fprintf('\tRecruited axons: %s\n', num2str(find(report.recruit > min(report.recruit))'));
%    fprintf('\tDismissed axons: %s\n', num2str(find(report.dismiss < max(report.dismiss))'));
% catch E
%    fprintf('\tCouldn''t print the report. Unexpected event: %s\n', E.message);
% end