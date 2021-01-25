function v_transition = gen_transitions(transition_cells, transitions, duration)
% Function that generates the transitions train

%Initialize
v_transition = zeros(duration, 1); 

% Place transition template onto start of transitiong group time
for i = length(transition_cells)
    % Start and end index of each group of transitioning spikes
    st_trans_gr  = transition_cells{i}(1);
    end_trans_gr = st_trans_gr + length(transitions{1}) - 1;
    
    % Vary the amplitude while inserting transition into spike train
    rand_amp = 0.99 + (1.01 - 0.99) .* rand; % small variation in amplitude
    v_transition(st_trans_gr : end_trans_gr, 1) = transitions{1}*rand_amp;
    
    %
end

% Makes sure it does not go over duration
v_transition(duration+1 : end) = [];

end

