% Script to compare extractor tool spike times with the simulated times

%% INSTRUCTIONS FOR USE
% 1) Choose either manual or automatic loading of files by commenting out
%    the unused one (line 19-33)
%
% 2) Adjust the positive and negative tolerance to suit your data(line 30-31). 
%   -If unsure, run the script once first and the code will tell you whether 
%     to increase/decrease if you were wrong
%   - OR run script and manually check the plots for the optimal tolerance
%     by checking the distance between each simulated and extracted sp location
%   WARNING: If a particular sim naxon suddenly has 'no matches' when you've
%   bumped up tolerance when it had some matches before, DONT PANIK, its
%   because the matching sp are like split 50-50 with another naxon and its
%   just flipping to another naxon when more sp are picked up
%
% 3) OPTIONAL enable Allow Overlap mode (line 38) which allows simulated sp which have 
%    detected >1 extraced sp to only save the closest one while the ignored ones are reported 
%    in the command line. Try not to do it when you first run and only as a
%    last resort when you see that the extracted sp are too close together
% 3.5) Set which plots and files you want to save.
% 4) Run script
%
% 5) Pray everyting goes right

close all
%clear

%% ADJUST TOLERANCE HERE
% The extracted sp must be between the positive and negative tolerance
% values of ONE simulated sp for it to be considered match
% Can be adjusted smaller to be more sensitive
pos_tolerance = 40; % Right of simulated sp (main one to adjust since extracted times(based on peak) are usually after simulated times(before peak))
neg_tolerance = 20; % Left of simulated sp (less important but still can be changed if extracted times somehow occurs before sim times)
fprintf('pos_tolerance = %d ; neg_tolerance = %d\n', pos_tolerance, neg_tolerance);

%% CHANGE OVERLAP MODE HERE
% This will allow simulated sp which have detected >1 extraced sp to save the
% closest one while others are reported in the command line. ONLY USE IF
% TOLERANCE CAN'T BE TUNED BETTER
allow_overlap = 1; % 0=disabled   1=enabled
if allow_overlap == 1; fprintf('<strong>OVERLAP MODE ENABLED</strong> simulated sp with >1 matching extracted sp will be allowed and reported\n'); end

%% CHANGE WHICH PLOTS WILL SHOW AND which files types to save
plot_raw_sptimes     = 0;
plot_matched_sptimes = 0;
plot_amp_noise       = 1;
plot_certainty_tbl   = 1;
save2excel = 0;
savereport = 0;

%% FILE DIRECTORIES / LOAD DATA (COMMENT OUT THE MANUAL OR AUTOMATIC SECTION AS NEEDED)
% (AUTO)Open simulated and extracted file automatically (Change file directories when needed)
% Simulation file
filename_sim    = 'simtest.mat';%'sim_artemio8.mat';%'simulated1.mat';
filepath_sim    = 'C:\Users\chris\Desktop\testfiles\sim';
% Extrac file
filename_extrac = 'extractest.mat';%'data_Art8.mat';%'extracted1.mat';
filepath_extrac = 'C:\Users\chris\Desktop\testfiles\extract';
% Report folder location (for saving report struct .mat)
filepath_report = 'C:\Users\chris\Desktop\testfiles\reports';
% Excel sheet file
filename_excel  = 'C:\Users\chris\Desktop\testfiles\compare_report.xlsx';


% Extracts data
sim_data = load(fullfile(filepath_sim, filename_sim));
fieldname = fieldnames(sim_data);
sim_data = sim_data.(fieldname{1});
extrac_data = load(fullfile(filepath_extrac, filename_extrac));
fieldname = fieldnames(extrac_data);
extrac_data = extrac_data.(fieldname{1});

% (MANUAL)Open simulated and extracted file manually using file UIs
%[sim_data, extrac_data] = open_sim_fileUI;

%% Set variables
end_loc       = length(sim_data.data); % End loc of the simulation
extracted_dt  = extrac_data.dt; % dt used to calculate sampling rate of extracted
simulated_loc = sim_data.report.locs; % simulated sp locations 
extracted_loc = time2idx_units(extrac_data.APstimes, extracted_dt); % extracted sp locations 
sim_          = sim_data.data; % Actual simulation data (all naxons with noise, drift etc)
sim_axons     = sim_data.axons; % Individual axon simulation (without noise)
sim_tempFamGroupings = sim_data.report.tempFamGroupings; % Groupings of sim axons based on family and shape

%% Plot simulated and extracted spike locations
if plot_raw_sptimes
    plot_sp_loc(simulated_loc, extracted_loc, end_loc, sim_);
    title(['RAW DATA COMPARISON' newline 'Top half is simulated and bottom half is extracted' newline 'Extracted is sorted into tiers based on shape and family']);
end

%% Compare
% Iterate through all extracted sp to find matching simulated sp
% Loops through each AP shape
for extrac_AP_num = 1:length(extracted_loc)
    % Loops through each family in AP shape
    for extrac_family_num = 1:length(extracted_loc{extrac_AP_num})
        % Loops through each simulated naxon
        matched_loc_temp=[];
        for sim_axon_group = sim_tempFamGroupings
            family_matched_axon_flag = 0;
            
            % Loops through each indivdual simulated sp loc
            for sim_family_num = 1:length(sim_axon_group.family_sptimes)
                for sim_sploc = sim_axon_group.family_sptimes{sim_family_num}'
                    % Determines the idx of matching sp which are within the
                    % tolerance range
                    idx = (extracted_loc{extrac_AP_num}{extrac_family_num} >= sim_sploc - neg_tolerance & ...
                           extracted_loc{extrac_AP_num}{extrac_family_num} <= sim_sploc + pos_tolerance   ...
                           );

                    % Used for overlap mode
                    if allow_overlap == 1 && sum(idx) > 1
                        ignored_overlap_sp = find(idx == 1, sum(idx)-1, 'last'); % records the sp which have been ignored
                        idx = find(idx == 1, 1, 'first'); % Adjust idx so that only the first closest sp is matched
                        ignored_str = [];
                        for ignore = ignored_overlap_sp
                            ignored_str = [ignored_str sprintf('extracted_loc{%d}{%d}(%d), ', extrac_AP_num, extrac_family_num, ignore)];
                        end
                        fprintf(['<strong>OVERLAP MODE ignored: </strong>' ignored_str newline]);
                    end


                    if ~isempty(extracted_loc{extrac_AP_num}{extrac_family_num}(idx))
                        try
                            % Checks if it the extracted sp was assigned to
                            % more than one sim sp (this doenst break the code
                            % but might mess up statistics)
                            if ~isempty(matched_loc_temp) && ismember(extracted_loc{extrac_AP_num}{extrac_family_num}(idx), matched_loc_temp(:,1))
                                fprintf('<strong>WARNING: </strong> extracted_loc{%d}{%d}(%d) assinged to >1 sim naxon. Try decreasing tolerance.\n', extrac_AP_num, extrac_family_num, find(idx==1));
                            end

                            % Records which simulated when there is a match into a
                            % temporary variable which is grouped after the for loop
                            matched_loc_temp(end+1,1) = extracted_loc{extrac_AP_num}{extrac_family_num}(idx);
                            matched_loc_temp(end,2) = sim_axon_group.Naxon; % Records which simulated naxon it belongs to
                            matched_loc_temp(end,3) = sim_axon_group.template; % Records matching sim shape/template num
                            matched_loc_temp(end,4) = sim_axon_group.family_num(sim_family_num); 

                        catch E
                            if sum(idx) > 1
                                error('A sim sp has detected >1 matching extracted sp within tolerance. @ extracted_loc{%d}{%d}(%d). Try decreasing tolerance or enable OVERLAP MODE', extrac_AP_num, extrac_family_num, find(idx==1, 1, 'first'));
                            else
                                rethrow(E)
                            end
                        end
                    end
                end
            end
            if ~isempty(matched_loc_temp)
                % Records the extracted sp that match with a simulated sp
                % into a larger group variable
                matched_loc{extrac_AP_num}{extrac_family_num, 1} = matched_loc_temp;
                matched_loc{extrac_AP_num}{extrac_family_num, 2} = unique(matched_loc_temp(:,2))'; % Records which simulated axon the family most likely belongs based on the the majority within the family
                matched_loc{extrac_AP_num}{extrac_family_num, 3} = unique(matched_loc_temp(:,3)); %Records template/shape
                family_matched_axon_flag = 1; 
            end
        end
        if family_matched_axon_flag == 0
            fprintf('<strong>Extracted AP %d Family %d doesnt have matching naxon.</strong> Try increasing tolerance to make sure all AP shapes are detected if possible\n', extrac_AP_num, extrac_family_num);
        end
    end
end

if plot_matched_sptimes
    try 
        % Plot compared data
        plot_sp_loc(simulated_loc, matched_loc, end_loc, sim_);
        title('EXTRACTED SPIKES THAT MATCH A SIM SPIKE');
    catch
        error('No matching spikes found (╯°□°）╯. Try increasing the tolerance variable to increase detection range');
    end
end

%% Compare report
report = compare_report(simulated_loc, extracted_loc, matched_loc, sim_tempFamGroupings, filename_sim, filename_extrac);

%% Amplitude & noise level plot
if plot_amp_noise
    amp_noise_plot(sim_axons, report, simulated_loc, sim_);
end

%% Confusion matrix
% Reorder matrix so that the diagonal contains the match
% if 0
%     [~, highest_certainty_idx] = sort(max(report.certainty, [], 2), 'descend');
%     report.certainty = report.certainty(:, highest_certainty_idx);
%     report.percent_family = report.percent_family(:, highest_certainty_idx);
%     report.col_headings = report.col_headings(:, highest_certainty_idx);
% end

% Makes column headings
col_headings = cell(1, size(report.col_headings, 2));
for certainty_match = 1:size(report.col_headings, 2)
    col_headings{certainty_match} = sprintf('match S%d F%d', report.col_headings(1,certainty_match), report.col_headings(2,certainty_match));
end

% Makes row headings
row_headings = cell(1, size(report.row_headings, 1));
for certainty_match = 1:size(report.row_headings, 1)
    row_headings{certainty_match} = sprintf('sim N%d S%d F%d', report.row_headings(certainty_match,1), report.row_headings(certainty_match,2), report.row_headings(certainty_match,3));
end

certainty_table = array2table(report.certainty, 'RowNames', row_headings, 'VariableNames', col_headings);
percent_family_table = array2table(report.percent_family, 'RowNames', row_headings, 'VariableNames', col_headings);

if plot_certainty_tbl
    ui_handle = uifigure('Position', [50 50 1000 600]);
    
    % Table for certainty
    ui_certainty_table = uitable(ui_handle,'Data', certainty_table, 'RowStriping', 0, 'Position', [10 300 970 250]);
    table_certainty_title = 'Certainty of each matched family to a sim family';
    uilabel(ui_handle, 'Text', table_certainty_title, 'Position', [10 570 970 20], 'FontSize', 16);
    
    % Table for percentage family
    ui_family_per_table = uitable(ui_handle,'Data', percent_family_table, 'RowStriping', 0, 'Position', [10 10 970 250]);
    table_family_per_title = 'Percentage of sim family taken up by that matched family';
    uilabel(ui_handle, 'Text', table_family_per_title, 'Position', [10 270 970 20], 'FontSize', 16);
    
    % FOMRATTING TABLES
    % Coloring the rows according to their naxon for easier readibility
    % (odd naxons colored white, odd naxons white)
    [even_row_idx] = find(~mod(report.row_headings(:,1), 2));
    odd_naxon_style = uistyle('BackgroundColor', [0.97 0.97 0.97]);
    addStyle(ui_certainty_table, odd_naxon_style, 'row', even_row_idx);
    addStyle(ui_family_per_table, odd_naxon_style, 'row', even_row_idx);
    
    for certainty_match = 0.2 : 0.2 : 0.8
        highlight_matches_style = uistyle('BackgroundColor', [1-certainty_match, 1, 1-certainty_match]);
        
        % Highlight matching families with certainties given by
        [match_row_idx, match_col_idx] = find(report.certainty >= certainty_match);
        addStyle(ui_certainty_table, highlight_matches_style, 'cell', [match_row_idx, match_col_idx]);
        
        % Highlight matching families based on percent family
        [match_row_idx_, match_col_idx_] = find(report.percent_family >= certainty_match);
        addStyle(ui_family_per_table, highlight_matches_style, 'cell', [match_row_idx_, match_col_idx_]);
    end
    
end

%% Save to Excel file
%enable_excel = input('Save report to Excel sheet? (Y/blank) :', 's');

if save2excel% && strcmpi(enable_excel,'y')
    writecell({'Sim file:'                   , filename_sim                      ,'' , 'Extrac_file:' , filename_extrac; ...
               report.match_accuracy         , 'total matched / total sim'   ,'' , datetime('now'), ''         ; ...
               report.correctly_identified_sp, 'total matched / total extrac','' , ''             , ''           ...
               }, ...
               filename_excel, 'Sheet', filename_sim, 'WriteMode', 'append');
    writetable(certainty_table,      filename_excel, 'Sheet', filename_sim, 'WriteMode', 'append', 'WriteRowNames', true, 'WriteVariableNames', true);
    writetable(percent_family_table, filename_excel, 'Sheet', filename_sim, 'WriteMode', 'append', 'WriteRowNames', true, 'WriteVariableNames', true);
end

%% Save report to .mat
if savereport
    filename_report = ['report_' filename_sim(1:end-4) '_' filename_extrac];
    save(fullfile(filepath_report, filename_report), 'report', '-mat');
end

%% Compare report function
function report = compare_report(simulated_loc, extracted_loc, matched_loc, sim_shape_family, filename_sim, filename_extrac)
total_num_simulated_sp = 0;
for naxon = 1:length(simulated_loc)
    total_num_simulated_sp = total_num_simulated_sp + length(simulated_loc{naxon});
end

total_num_extracted_sp = 0;
for extrac_AP_num = 1:length(extracted_loc)
    for extrac_family_num = 1:size(extracted_loc{extrac_AP_num}, 1)
        total_num_extracted_sp = total_num_extracted_sp + size(extracted_loc{extrac_AP_num}{extrac_family_num, 1}, 1);        
    end
end

total_num_matched_sp = 0;
for extrac_AP_num = 1:length(matched_loc)
    for extrac_family_num = 1:size(matched_loc{extrac_AP_num}, 1)
        total_num_matched_sp = total_num_matched_sp + size(matched_loc{extrac_AP_num}{extrac_family_num, 1}, 1);        
    end
end

% Checks if total number of extracted > total number of sim. This does
% cause some statisitcal ratios to be a bit weird but no biggie 
if total_num_extracted_sp > total_num_simulated_sp; fprintf('<strong>WARNING:</strong> total number of extracted > total number of sim\n'); end

match_accuracy = total_num_matched_sp / total_num_simulated_sp * 100;
fprintf('<strong>%.1f%%</strong> (%d/%d) total matched sp / total simulated sp.\n', match_accuracy, total_num_matched_sp, total_num_simulated_sp);

correctly_identified_sp = total_num_matched_sp / total_num_extracted_sp * 100;
fprintf('<strong>%.1f%%</strong> (%d/%d) total matched sp / total extracted sp.\n', correctly_identified_sp, total_num_matched_sp, total_num_extracted_sp);

% Calculates certainty and 
row = 0; 
col = 0; 
certainty = [];
percent_family = [];
for match_AP_num = 1:length(matched_loc)
    for match_fam_num = 1:size(matched_loc{match_AP_num}, 1)
        col = col + 1;
        col_headings(1, col) = match_AP_num;
        col_headings(2, col) = match_fam_num;
        for sim_axon_num = 1:length(sim_shape_family)
            n=1;
            for sim_fam_num = sim_shape_family(sim_axon_num).family_num
                row = row + 1;
                row_headings(row, 1) = sim_axon_num;
                row_headings(row, 2) = sim_shape_family(sim_axon_num).template;
                row_headings(row, 3) = sim_fam_num;
                
                matched_axon = matched_loc{match_AP_num}{match_fam_num, 1}(:,2) == sim_axon_num;
                matched_family = matched_loc{match_AP_num}{match_fam_num, 1}(:,4) == sim_fam_num;
                % number of sp in a matched family that ACTUALLY match to a
                % specific simulated family
                num_matched = sum(matched_axon & matched_family);
                
                certainty(row, col) = num_matched / size(matched_loc{match_AP_num}{match_fam_num}, 1);
                
                percent_family(row, col) = num_matched / size(sim_shape_family(sim_axon_num).family_sptimes{n}, 1);
                
                percent_naxon_temp(row, col) = num_matched / length(simulated_loc{sim_axon_num});
                
                n=n+1;
            end
        end
        row = 0;
    end
    
end

percent_naxon = zeros(1,length(simulated_loc));
for i = 1:length(simulated_loc)
    percent_naxon(i) = sum(percent_naxon_temp(row_headings(:,1)==i, :), 'all');
end

report.row_headings   = row_headings;
report.col_headings   = col_headings;
report.certainty      = certainty;
report.percent_family = percent_family;
report.percent_naxon  = percent_naxon;
report.match_accuracy = match_accuracy;
report.correctly_identified_sp = correctly_identified_sp;
report.filename_sim    = filename_sim;
report.filename_extrac = filename_extrac;

% for naxon = 1:length(sim_shape_family)
%     for sim_family_num = sim_shape_family(naxon).family_num
%         family_matched_axon_flag = 0; % Set a marker/flag which indicates whether the naxon succesfully matched a family
% 
%         for extrac_AP_num = 1:length(matched_loc)
%             for extrac_family_num = 1:size(matched_loc{extrac_AP_num}, 1)
%                 for associated_axon = matched_loc{extrac_AP_num}{extrac_family_num, 2}
%                     if naxon == associated_axon
%                         % num of xtracted sp in the family matched to the same naxon
%                         extracted_matched = sum(matched_loc{extrac_AP_num}{extrac_family_num, 1}(:,2) == naxon);
%                         % extracted_matched / num of sp in family
%                         percent_family = extracted_matched / size(matched_loc{extrac_AP_num}{extrac_family_num}, 1);
%                         % extracted_matched / num of sp in sim naxon
%                         percent_naxon = extracted_matched / length(simulated_loc{naxon});
% 
%                         % Print the row onto cmd window
%     %                     fprintf(['  %d   |   %d         %d        %0.2f        %0.2f' newline], ...
%     %                              naxon,           ...
%     %                              extrac_AP_num,          ...
%     %                              extrac_family_num,      ...
%     %                              percent_family,  ...
%     %                              percent_naxon    ...
%     %                              );
%                         % Record into table variable
%                         report.table(row, 1) = naxon;
%                         report.table(row, 2) = extrac_AP_num;
%                         report.table(row, 3) = extrac_family_num;
%                         report.table(row, 4) = percent_family;
%                         report.table(row, 5) = percent_naxon;
% 
%                         report.percent_naxon(naxon) = report.percent_naxon(naxon) + percent_naxon; 
% 
%                         row = row + 1;
%                         family_matched_axon_flag = 1; % Set a marker that a family was found to match the naxon
%                     end
%                 end %if naxon == matched_loc{AP_num}{family_num, 2}
%             end %for family_num = 1:size(matched_loc{AP_num}, 1)
%         end %for AP_num = 1:length(matched_loc)
% 
%         if ~family_matched_axon_flag
%             % Print the row if there are no matching families to the simulated naxon
%             fprintf(['  %d   |              no match' newline], naxon);
%             fprintf(row_border);
%         else
%             fprintf(row_border);
%         end
%     end
% end %for naxon = 1:length(simulated_loc)

end

%% Plot sim&extrac and sim&matched stem plots
function plot_sp_loc(simulated_loc, extacORmatch_loc, end_loc, sim_)

fig_handle = figure;
clf(fig_handle);
hold on;

% Plot simlated sp locations
for naxon = 1:length(simulated_loc)
    naxon_tier = ((naxon-1) / length(simulated_loc)) * 0.1; % Factor used to visually separate each naxon on the plot
    stem(simulated_loc{naxon}, ones(size(simulated_loc{naxon})) - naxon_tier);
end

% Plot actual simulation
sim_plot = plot(sim_/max(sim_));
sim_plot.Color(4) = 0.3;

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

%% Open the files using file UI
function [sim_data, extrac_data]= open_sim_fileUI

try
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

catch
    error('Manually closed. User did not select file.');
end

end

%% 
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

%% Plot amplitude of naxons and noise levels function
function amp_noise_plot(sim_axons, report, simulated_loc, sim_)
figure;
hold on;
max_peak_amps = max(sim_axons);
for i = 1 : size(sim_axons, 2)
    [peaks, locs] = findpeaks(sim_axons(:,i), 'MinPeakHeight', 0.1, 'MinPeakProminence', max_peak_amps(i)/3);
    %pks{i} = [locs peaks];
    plot(locs, peaks/max(max_peak_amps), 'o', 'MarkerSize', 1.5);
    title('Amplitude of each naxon');
    %percent_naxon = NaN;%sum(report.percent_family(report.row_headings(:,1)==i, :), 'all');
    naxon_legend{i} = ['% Naxon ' num2str(i) ': ' num2str(report.percent_naxon(i))]; 
end
axis([0 size(sim_axons, 1) 0 1]);
% Plot the noise standard deviation levels [HAS ITS FLAWS WIP]
offset = 0; 
if simulated_loc{1,1}(1) <= 1100; offset = simulated_loc{1,1}(1)+100; end % Offset incase there is a spike within the first 1000 idxs
noise_std = std(sim_((1:1000)+offset)) / max(max_peak_amps);
plot([0 length(sim_)], [noise_std noise_std], 'r-'); % 1 std
plot([0 length(sim_)], 2*[noise_std noise_std], 'g-'); % 2 std
plot([0 length(sim_)], 3*[noise_std noise_std], 'b-'); % 3 std
legend(naxon_legend);

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