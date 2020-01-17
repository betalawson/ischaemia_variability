function [WL, WL_mean, WL_std] = calcWavelengthMap( simdata )
% This function takes the input simulation and creates a map of the
% wavelength at each of the activated sites, using an estimated conduction
% velocity and the APD values recorded for those activations

% Temporal threshold that repolarisation and activation events must fall
% within to represent the front/back of the same wave
threshold = 0.2;     % (in ms)

% Read out the activation times from the provided simulation data
AT = simdata.activation_times(:,1);

% Work with only the active sites
act_sites = AT > 0;
AT = AT(act_sites);

% Read out the index to map back to original
index = find(act_sites);

% Read out the co-ordinates
cx = double( simdata.coords(act_sites,1) );
cy = double( simdata.coords(act_sites,2) );

% Now read out APDs and calculate repolarisation times
APD = simdata.APDs(act_sites,1);
RT = AT + APD;

% Only calculate for those sites 30% or further into the domain
calc_sites = simdata.coords(act_sites,1) > 0.3 * max(simdata.coords(act_sites,1) );
calc_indices = find(calc_sites);

% Loop over activated sites
WL = NaN( length(simdata.activation_times(:,1) ), 1 );
for m = 1:length(calc_indices)
    
    k = calc_indices(m);
    
    % Find the sites with a difference between their repolarisation time
    % and this site's activation time below a threshold
    mask = ( abs( AT(k) - RT ) < threshold );
    
    % Take the minimum distance of the sites falling within the threshold
    %dists =  sqrt( (cx(k) - cx(mask) ).^2 + ( cy(k) - cy(mask) ).^2 );
    dists = abs( cx(k) - cx(mask) ) + abs( cy(k) - cy(mask) );
    
    % Only calculate wavelength when 
    if ~isempty(dists)
        WL( index(k) ) = min(dists);
    end
        
end

% Output the mean and standard deviation of wavelength values
WL_mean = mean( WL, 'omitnan' );
WL_std = std( WL, 'omitnan' );