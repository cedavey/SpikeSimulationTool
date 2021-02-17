function v_transition = gen_transitions(transition_cells, transitions, duration, abs_refrac)
% Function that generates the transitions train

%Initialize
v_transition = zeros(duration, 1); 

% Place transition template onto start of transitiong group time
for i = 1 : length(transition_cells)
    
    % isi between first two spikes in each group minus the absolute
    % refractory period used to select appropriate transition index
    transition_idx = transition_cells{i}(2) - transition_cells{i}(1) - abs_refrac;
    
    % Always places the first two transitioning spike regardless of the number
    % of spikes in the group (3 or more transitions handled later)
    % Start and end index of each group of transitioning spikes
    st_trans_gr  = transition_cells{i}(1);
    end_trans_gr = st_trans_gr + length(transitions{transition_idx}) - 1;
    % Vary the amplitude while inserting transition into spike train
    rand_amp = 0.99 + (1.01 - 0.99) .* rand; % small variation in amplitude
    v_transition(st_trans_gr : end_trans_gr, 1) = transitions{transition_idx}*rand_amp;
    
    % Interpolates 3 or more transitioning spikes since there is no
    % premade template for them
    if length(transition_cells{i}) >= 3
        for n = 3 : (length(transition_cells{i}(3 : end)) + 2) % Iterate from 3rd transitioning spike onwards
            
            t_between_trans_sp = transition_cells{i}(n) - transition_cells{i}(n-1); % Time between the spikes
            transition_idx = t_between_trans_sp - abs_refrac; % Minus abs_refrace to get the idx of the trans template being used
            
            % Extracts only the transitioning spike (2nd one in the template instead of the full double spike)
            transition_template = transitions{transition_idx}(t_between_trans_sp : end);
            
            st_trans_gr  = transition_cells{i}(n);
            end_trans_gr = st_trans_gr + length(transition_template) - 1;
            
            % Vary the amplitude while inserting transition into spike train
            rand_amp = 0.99 + (1.01 - 0.99) .* rand; % small variation in amplitude
            v_transition(st_trans_gr : end_trans_gr, 1) = transition_template*rand_amp;
            
            % Smoothen (interpolate) the transition since the connecting
            % parts may not perfectly line up [THIS PART IS A BIT OOFT]
            % This part may need some improvement but it works for now
             
            
            % [Imagine its like lines on the graph and the interpolated region is between t2 and t4]   
            %                      /
            %                    /
            %                  /
            % ______         
            %  t1 t2 [interp]  t3 t4
            
            t1 = st_trans_gr - 2; % 2 time steps before region of interpolation
            t2 = st_trans_gr - 1; % 1 time step before region of interpolation
            
            t3 = find(transition_template > v_transition(st_trans_gr-1), 1, 'first') + st_trans_gr; % 1 time step after region of interpolation where the data point at t2 is larger than the data point at t3
            t4 = t3 + 1; % 2 time steps after region of interpolation
            
            t = [t1; t2; t3; t4]; % Combine time steps into time matrix for interpolation
            
            d = v_transition([t1, t2, t3, t4]); % Extract end bits (data points) where interpolation will join the ends
                        
            v_transition(t1:t4) = interp1(t, d, t1:t4, 'makima'); % Replace with the smoothened (interpolated) version
            
        end
    end
end

% Makes sure it does not go over duration
v_transition(duration+1 : end) = [];

end

