% Script to compare extractor tool spike times with the simulated times

%% INSTRUCTIONS FOR USE
% 1) Choose either manual or automatic loading of files by commenting out
%    the unused one (line 19-33)
%
% 2) Adjust the positive and negative tolerance to suit your data(line 39-40). 
%   -If unsure, run the script once first and the code will tell you whether 
%     to increase/decrease if you were wrong
%   - OR run script and manually check the plots for the optimal tolerance
%     by checking the distance between each simulated and extracted sp location
%   WARNING: If a particular sim naxon suddenly has 'no matches' when you've
%   bumped up tolerance when it had some matches before, DONT PANIK, its
%   because the matching sp are like split 50-50 with another naxon and its
%   just flipping to another naxon when more sp are picked up
%
% 3) OPTIONAL enable Allow Overlap mode which allows simulated sp which have 
%    detected >1 extraced sp to only save the closest one while others are reported 
%    in the command line. Try not to do it when you first run and only as a
%    last resort when tolerance doesnt work
%
% 4) Run script
%
% 5) Pray everyting goes right
%

close all
clear
%% LOAD DATA (COMMENT OUT THE MANUAL OR AUTOMATIC SECTION AS NEEDED)
% (MANUAL)Open simulated and extracted file manually using file UIs
%[sim_data, extrac_data] = open_sim_fileUI;

% (AUTO)Open simulated and extracted file automatically (Change file directories when needed)
file_sim    = 'sim_C&D4.mat';%'simulated1.mat';%'sim_artemio8.mat';
path_sim    = 'C:\Users\chris\Desktop\sim test';
file_extrac = 'data_C&D4.mat';%'extracted1.mat'; %'data_Art8.mat';
path_extrac = 'C:\Users\chris\Desktop\sim data';
sim_data = load(fullfile(path_sim, file_sim));
fieldname = fieldnames(sim_data);
sim_data = sim_data.(fieldname{1});
extrac_data = load(fullfile(path_extrac, file_extrac));
fieldname = fieldnames(extrac_data);
extrac_data = extrac_data.(fieldname{1});

%% ADJUST TOLERANCE HERE
% The extracted sp must be between the positive and negative tolerance
% values of ONE simulated sp for it to be considered match
% Can be adjusted smaller to be more sensitive
pos_tolerance = 90; % Right of simulated sp (main one to adjust since extracted times(based on peak) are usually after simulated times(before peak))
neg_tolerance = 0; % Left of simulated sp (less important but still can be changed if extracted times somehow occurs before sim times)
fprintf('pos_tolerance = %d ; neg_tolerance = %d\n', pos_tolerance, neg_tolerance);

%% CHANGE OVERLAP MODE HERE
% This will allow simulated sp which have detected >1 extraced sp to save the
% closest one while others are reported in the command line. ONLY USE IF
% TOLERANCE CAN'T BE TUNED BETTER
allow_overlap = 1; % 0=disabled   1=enabled
if allow_overlap == 1; fprintf('<strong>OVERLAP MODE ENABLED</strong> simulated sp with >1 matching extracted sp will be allowed and reported\n'); end

%% Set variables
end_loc       = length(sim_data.data); % End loc of the simulation
extracted_dt  = extrac_data.dt; % dt used to calculate sampling rate of extracted
simulated_loc = sim_data.report.locs; % simulated sp locations 
extracted_loc = time2idx_units(extrac_data.APstimes, extracted_dt); % extracted sp locations 

%% Plot simulated and extracted spike locations
plot_sp_loc(simulated_loc, extracted_loc, end_loc);
title(['RAW DATA COMPARISON' newline 'Top half is simulated and bottom half is extracted' newline 'Extracted is sorted into tiers based on shape and family']);

%% Compare
% Iterate through all extracted sp to find matching simulated sp
% Loops through each AP shape
for AP_num = 1:length(extracted_loc)
    % Loops through each family in AP shape
    for family_num = 1:length(extracted_loc{AP_num})
        % Loops through each simulated naxon
        matched_loc_temp=[];
        for sim_axon_num = 1:length(simulated_loc)
            family_matched_axon_flag = 0;
            
            % Loops through each indivdual simulated sp loc
            for sim_sploc = simulated_loc{sim_axon_num}'
                % Determines the idx of matching sp which are within the
                % tolerance range
                idx = (extracted_loc{AP_num}{family_num} >= sim_sploc - neg_tolerance & ...
                       extracted_loc{AP_num}{family_num} <= sim_sploc + pos_tolerance   ...
                       );
                
                % Used for overlap mode
                if allow_overlap == 1 && sum(idx) > 1
                    ignored_overlap_sp = find(idx == 1, sum(idx)-1, 'last'); % records the sp which have been ignored
                    idx = find(idx == 1, 1, 'first'); % Adjust idx so that only the first closest sp is matched
                    ignored_str = [];
                    for ignore = ignored_overlap_sp
                        ignored_str = [ignored_str sprintf('extracted_loc{%d}{%d}(%d), ', AP_num, family_num, ignore)];
                    end
                    fprintf(['<strong>Overlap ignored sim naxon %d: </strong>' ignored_str newline], sim_axon_num);
                end
                   
                % Records which simulated when there is a match
                if ~isempty(extracted_loc{AP_num}{family_num}(idx))
                    try
                        matched_loc_temp(end+1,1) = extracted_loc{AP_num}{family_num}(idx);
                        matched_loc_temp(end,2) = sim_axon_num; % Records which simulated naxon it belongs to
                    catch
                        if sum(idx) > 1
                            error('A sim sp has detected >1 matching extracted sp within tolerance. @ extracted_loc{%d}{%d}(%d). Try decreasing tolerance or enable OVERLAP MODE', AP_num, family_num, find(idx==1, 1, 'first'));
                        else
                            error('I dun know whats wrong ¯\_( ͡ಥ ͜ʖ ͡ಥ)_/¯');
                        end
                    end
                end
            end
            if ~isempty(matched_loc_temp)
                % Records the extracted sp that match with a simulated sp
                matched_loc{AP_num}{family_num, 1} = matched_loc_temp;
                matched_loc{AP_num}{family_num, 2} = unique(matched_loc_temp(:,2))'; % Records which simulated axon the family most likely belongs based on the the majority within the family
                
                family_matched_axon_flag = 1; 
            end
        end
        if family_matched_axon_flag == 0
            fprintf('<strong>Extracted AP %d Family %d doesnt have matching naxon.</strong> Try increasing tolerance to make sure all AP shapes are detected if possible\n', AP_num, family_num);
        end
    end
end

try 
    % Plot compared data
    plot_sp_loc(simulated_loc, matched_loc, end_loc);
    title('EXTRACTED SPIKES THAT MATCH A SIM SPIKE');
catch
    error('No matching spikes found (╯°□°）╯︵ ┻━┻. Try increasing the tolerance variable to increase detection range');
end


%% Compare report
total_num_simulated_sp = 0;
for naxon = 1:length(simulated_loc)
    total_num_simulated_sp = total_num_simulated_sp + length(simulated_loc{naxon});
end

total_num_extracted_sp = 0;
for AP_num = 1:length(extracted_loc)
    for family_num = 1:size(extracted_loc{AP_num}, 1)
        total_num_extracted_sp = total_num_extracted_sp + size(extracted_loc{AP_num}{family_num, 1}, 1);        
    end
end

total_num_matched_sp = 0;
for AP_num = 1:length(matched_loc)
    for family_num = 1:size(matched_loc{AP_num}, 1)
        total_num_matched_sp = total_num_matched_sp + size(matched_loc{AP_num}{family_num, 1}, 1);        
    end
end

% Checks if total number of extracted > total number of sim. This does
% cause some statisitcal ratios to be a bit weird but no biggie 
if total_num_extracted_sp > total_num_simulated_sp; fprintf('<strong>WARNING:</strong> total number of extracted > total number of sim\n'); end

match_accuracy = total_num_matched_sp / total_num_simulated_sp * 100;
fprintf('<strong>%.1f%%</strong> of spikes were correctly identified and matched to the original simulation.\n', match_accuracy);

incorrectly_identified_sp = (total_num_extracted_sp - total_num_matched_sp) / total_num_extracted_sp * 100;
fprintf('<strong>%.1f%%</strong> of spikes were out of tolerance range / ignored(if overlap mode) based on the total extracted.\n', incorrectly_identified_sp);

% Prints a table which identifies the extracted spikes belonging to a
% certain simulated axon
row_border = '-----------------------------------------------\n';
fprintf([' SIM  |               EXTRAC'                       newline ...
         'Naxon |  AP      Family    Certainty    %%naxon'    newline ...
                              row_border                              ...
         ]);
for naxon = 1:length(simulated_loc)
    family_matched_axon_flag = 0; % Set a marker/flag which indicates whether the naxon succesfully matched a family
    
    for AP_num = 1:length(matched_loc)
        for family_num = 1:size(matched_loc{AP_num}, 1)
            for associated_axon = matched_loc{AP_num}{family_num, 2}
                if naxon == associated_axon
                    % certainty is the number of identified spikes in a family which actually belong to the simulated axon group (since the
                    % the associated axon is determined by the majority of spikes within  the family, some may be incorrectly assigned to the wrong naxon)
                    % certainty < 1.00 means that at least one of the sp in the family is matched to a different naxon
                    certainty = sum(matched_loc{AP_num}{family_num, 1}(:,2) == naxon) / size(matched_loc{AP_num}{family_num}, 1);
                    % percent_of_naxon is the fraction of sp within the family over the total simulated sp in the naxon
                    % it is multiplied by certainty since not all sp within the family may come from the same naxon
                    percent_of_naxon = (size(matched_loc{AP_num}{family_num, 1}, 1) / length(simulated_loc{naxon})) * certainty;
                    
                    % Print the row
                    fprintf(['  %d   |   %d         %d        %0.2f        %0.2f' newline], ...
                             naxon,           ...
                             AP_num,          ...
                             family_num,      ...
                             certainty,       ...
                             percent_of_naxon ...
                             );
                    
                    family_matched_axon_flag = 1; % Set a marker that a family was found to match the naxon
                end
            end %if naxon == matched_loc{AP_num}{family_num, 2}
        end %for family_num = 1:size(matched_loc{AP_num}, 1)
    end %for AP_num = 1:length(matched_loc)
    
    if family_matched_axon_flag == 0
        % Print the row if there are no matching families to the simulated naxon
        fprintf(['  %d   |              no match' newline], naxon);
        fprintf([row_border]);
    else
        fprintf([row_border]);
    end

end %for naxon = 1:length(simulated_loc)

%%
function plot_sp_loc(simulated_loc, extacORmatch_loc, end_loc)

figure;
hold on;

% Plot simlated sp locations
for naxon = 1:length(simulated_loc)
    naxon_tier = ((naxon-1) / length(simulated_loc)) * 0.1; % Factor used to visually separate each naxon on the plot
    stem(simulated_loc{naxon}, ones(size(simulated_loc{naxon})) - naxon_tier);
end

% Plot extracted sp locations
for AP_num = 1:length(extacORmatch_loc)
    for family_num = 1:size(extacORmatch_loc{AP_num}, 1)
        % These tier values are used to visually categorize and differentiate the AP
        % shape and family from each other by putting them in different
        % layers/tiers on the plot
        AP_tier     = (AP_num - 1) * 0.1; % Each new AP shape occupies a new tier with height of 0.1
        family_tier = (family_num / size(extacORmatch_loc{AP_num}, 1)) * 0.1; % Factor that plots the each family on a separate tier level on the plot to differentiate the families
        
        stem(extacORmatch_loc{AP_num}{family_num}(:,1), 0.5*ones(size(extacORmatch_loc{AP_num}{family_num}(:,1))) - family_tier - AP_tier);
    end
end

% Plot a line indicating the top and bottom half of the figure
plot([0 end_loc], [0.6 0.6], 'g');

end

function [sim_data, extrac_data]= open_sim_fileUI
% Open using simulation
[file, path] = uigetfile('*.mat', 'Select simulation');
fullName = fullfile(path, file);
sim_data = load(fullName);
fieldname = fieldnames(sim_data);
sim_data = sim_data.(fieldname{1});

% Open extracted spikes
[file, path] = uigetfile('*.mat', 'Select extracted spikes');
fullName = fullfile(path, file);
extrac_data = load(fullName);
fieldname = fieldnames(extrac_data);
extrac_data = extrac_data.(fieldname{1});

end

function extracted_loc = time2idx_units(extracted_loc, extracted_dt)

% Convert the extracted spike locs which are originally in units of time[s]
% to index units
for AP_num = 1:length(extracted_loc)
    for family_num = 1:length(extracted_loc{AP_num})
        temp = round(extracted_loc{AP_num}{family_num} / extracted_dt); % Convert from time to idx by dividing by dt
        
        if size(temp,2) > 1; temp = temp'; end % Sometimes the input is row and sometimes its column vector, this will turn it into column vector only
        
        extracted_loc{AP_num}{family_num} = temp;
    end
end

end

% function associated_axon_list = associate_family2axon(associated_axon_of_each_sp)
% % Sorts the into most associated axons to least associated axon in the
% % family
% 
% associated_axons = unique(associated_axon_of_each_sp');
% 
% if length(associated_axons) == 1
%     associated_axon_list = associated_axons;
% else
%     [n, axons] = hist(associated_axon_of_each_sp, associated_axons);
% 
%     [~, idx] = sort(n, 'descend');
% 
%     associated_axon_list = axons(idx);
% end
% 
% end

% function simulated_loc_out = row2col(simulated_loc_in)
% 
% % Converts the simlated sp loc data from rows to columns so they match the
% % extracted data
% for axon_num = 1:length(simulated_loc_in)
%     simulated_loc_out{axon_num} = simulated_loc_in{axon_num}';
% end
% 
% end