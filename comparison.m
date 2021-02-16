% Script to compare extractor tool spike times with the simulated times

close all
clear

% Open simulation
[file, path] = uigetfile('*.mat');
fullName = fullfile(path, file);
load(fullName);

for i = 1:size(simulation.report.spks, 2)
    axon_simulated(:,i) = simulation.report.spks(:,i);
end

% Open extracted spikes
[file, path] = uigetfile('*.mat');
fullName = fullfile(path, file);
load(fullName);

for i = 1:size(simulation.report.spks, 2)
    axon_extracted{:,i} = (vsim_Rescale_APs_spikes.APstimes{1,1}{i,1})';
end