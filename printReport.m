% 17th December 2020
% Chi Yung Darren Tan 910828

% This function replaces the print report within SpikeSimulationTool.m

% Function is called within SST.mIapp and will print the report within
% app.reportTextArea in the GUI

function report_op = printReport(report, dt)

try
   if isfield(report, 'inf_time')
      report_op = sprintf('\tInflammation: %.02f s | Number of inflamed axons: %d\n\tAmplitude change time: %.02f s\n\tRecruited axons: %s\n\tDismissed axons: %s\n', ...
          report.inf_time * dt, numel(report.inflamed), report.opts.Events.amplitude_dist_onset * dt, num2str(find(report.recruit > min(report.recruit))'), ...
          num2str(find(report.dismiss < max(report.dismiss))'));
   else
      report_op = sprintf('\tNumber of inflamed axons: %d\n\tAmplitude change time: %.02f s\n\tRecruited axons: %s\n\tDismissed axons: %s\n', ...
          0, report.opts.Events.amplitude_dist_onset * dt, num2str(find(report.recruit > min(report.recruit))'), num2str(find(report.dismiss < max(report.dismiss))'));
   end
catch E
   report_op = sprintf('\tCouldn''t print the report. Unexpected event: %s\n', E.message);
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