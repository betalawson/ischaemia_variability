function wave_disorder = calcWaveDisorder(simdata)
% This function calculates the "wave disorder", the inconsistency in the
% wavefront's shape, for the input simulation data

% Wavespeed is approximated in a macroscopic sense, detecting at what
% moments in time at least <activation_proportion>% of the sites at certain 
% horizontal locations have been activated. These locations are defined in
% terms of a fraction of the total domain - start_frac and end_frac
activation_proportion = 0.2;
start_frac = 0.35;
end_frac = 0.9;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a list of all sites that actually recorded an activation
activated = (simdata.N_activations > 0);

% Read out the x co-ordinates and activation times of all activated sites
x = double( simdata.coords(activated,1) );
t_obs = simdata.activation_times(activated,1);

% Read out domain length
Lx = max(x);

% Starting position for wavespeed measurement is 35% of the way along the
% domain, ending point is 90%. Find closest x locations to these locations
% and use those
dist_x = abs( x - start_frac*Lx );
[~, I] = sort(dist_x, 'ascend');
start_x = x(I(1));
start_sites = (x == start_x);

dist_x = abs( x - end_frac*Lx );
[~, I] = sort(dist_x, 'ascend');
end_x = x(I(1));
end_sites = (x == end_x);

% Now find moments in time when <activation_proportion>% of these become
% active
start_times = simdata.activation_times(start_sites,1);
start_times = sort(start_times, 'ascend');
quantile_site = ceil( activation_proportion * length(start_times) );
start_t = start_times( quantile_site );

end_times = simdata.activation_times(end_sites,1);
end_times = sort(end_times, 'ascend');
quantile_site = ceil( activation_proportion * length(end_times) );
end_t = end_times( quantile_site );


% Use these to calculate the macroscopic velocity
v_macro = (end_x - start_x) / (end_t - start_t);

% Grab out all locations and times from sites within the range used to
% measure the macroscopic velocity
x_use = x( (x >= start_x) & (x <= end_x) );
t_use = t_obs( (x >= start_x) & (x <= end_x) );

% Calculate the total discrepancy in time from the times predicted by this
% optimal time
wave_disorder = rms( x_use / v_macro - t_use );

end

