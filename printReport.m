% 17th December 2020

% This function replaces the print report within SpikeSimulationTool.m

% Function is called within SST.mIapp and will print the report within
% app.reportTextArea in the GUI

function report_op = printReport(report, dt)

% Determines which axons are recruited and dismissed for the report
axons_recruited = sprintf('%d, ', find(report.recruit > min(report.recruit)));
axons_recruited = axons_recruited(1 : end-2); % removes final ', '

axons_dismissed = sprintf('%d, ', find(report.dismiss < max(report.dismiss)));
axons_dismissed = axons_dismissed(1 : end-2); % removes final ', '

% Print report depending on whether there is inflammation
try
   if isfield(report, 'inf_time')
      report_op = sprintf('Start time of inflammation: %.02f s | Number of inflamed axons: %d\nAmplitude change time: %.02f s\nRecruited axon(s): %s\nDismissed axon(s): %s\n', ...
                          report.inf_time * dt,                         ... % Start time of inflammation [s]
                          numel(report.inflamed),                       ... % Num of inflammed axons
                          report.opts.Events.amplitude_dist_onset * dt, ... % Amplitude change time [s]
                          axons_recruited,                              ... % Axons Recruited
                          axons_dismissed                               ... & Axons Dismissed
                          );
   else
      report_op = sprintf('Number of inflamed axons: %d\nAmplitude change time: %.02f s\nRecruited axon(s): %s\nDismissed axon(s): %s\n', ...
                          0,                                            ... % Num of inflammed axons
                          report.opts.Events.amplitude_dist_onset * dt, ... % Amplitude change time [s]
                          axons_recruited,                              ... % Axons Recruited
                          axons_dismissed                               ... & Axons Dismissed
                          );
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