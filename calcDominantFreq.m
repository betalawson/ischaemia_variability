function [DF, F_std] = calcDominantFreq(simdata)
% This function takes an input 2D simulation and calculates the dominant
% frequency of re-entry. If the majority of sites display no secondary
% activations, the dominant frequency will be recorded as NaN

bin_width = 20;

% First, convert all -1 values in the activation times matrix to be NaN.
% This way no comparisons can be made after the final activation for a site
act_times = simdata.activation_times;
act_times( act_times == -1 ) = NaN;

% Now, subtract these times to get the periods of activation
periods = diff(act_times,1,2);

% Find the most common period recorded, and invert it to get the dominant
% frequency
centres = bin_width/2:bin_width:1000-bin_width/2;
counts = hist( periods(:), centres);
[~, loc] = max(counts);
DF = 1 / centres(loc) * 1000;

% Find the standard deviation of frequency
F = 1 ./ periods(:) * 1000;
F_std = std(F, 'omitnan');