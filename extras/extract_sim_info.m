duration = vsim.time(end);

recruit = vsim.report.recruit * vsim.dt;
recruit(recruit < 0.0003) = 0;

dis = vsim.report.dismiss * vsim.dt;

inflamed = zeros(size(recruit));
infl = vsim.report.inflamed;
if infl ~= 0
   inf_time = vsim.report.inf_time * vsim.dt;
   inflamed(infl) = inf_time;
end

Nspikes = sum(vsim.report.spks);
sr = vsim.report.opts.SpikeRate;

amped = vsim.report.amped;

disp('--------------------------------');
disp(['  Duration:' num2str(duration)]);
disp(['Spike rate: ' num2str(sr)]);
disp(['  # Spikes: ' num2str(Nspikes)]);
disp([' Recruited: ' num2str(recruit)]);
disp([' Dismissed: ' num2str(dis)]);
disp(['  Inflamed: ' num2str(inflamed)]);
disp(['     Amped: ' num2str(amped)]);
disp('--------------------------------');