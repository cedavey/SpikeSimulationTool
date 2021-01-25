function [sptimes, non_transition, transition] = separate_transition_spikes(isi, duration_of_spike)
% Seperates the spikes which are close together and have transitions between them and
% those far apart without transitions.

% Calculate times of all spikes
sptimes = cumsum(isi);

% Determine the times for spikes that do not transition
non_transition    = sptimes(1);
for i = 1 : length(isi)-1
   if isi(i) <= duration_of_spike  || isi(i+1) <= duration_of_spike
       % When the isi for the next 2 spikes is <= rest, adds the 
       % current isi to the current cumulative sum
       non_transition(end,1) = sptimes(i+1); 
   elseif isi(i) > duration_of_spike
       % When the isi for the spikes is > rest, it creates a new index to
       % indicate the next normal (without transition) spike
       non_transition(end+1,1) = sptimes(i+1);
   else
       error('Error: There is a missing value in sptimes_normal');
   end
   
end
% Start and end exceptions that the for loop cant handle properly
% Since the first spike doesnt have a spike before it, it can be
% non-transitioning when the isi after it is > rest 
if ( isi(2) > duration_of_spike ) && ( non_transition(1) ~= sptimes(1) )
    non_transition = [sptimes(1); non_transition]; 
end
% Clears the last sptimes_normal if the last isi is less than the rest time
% since the for loop above can't check the preallocated last value (may
% need to change this when considering total duration of simulation)
if isi(end) <= duration_of_spike 
    non_transition(end) = []; 
end

% Determines the times of spikes with transitions (mutually exclusive to
% times without transitions)
transition = sptimes(~ismember(sptimes, non_transition));

end