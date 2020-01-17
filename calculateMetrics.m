function metrics = calculateMetrics( simdata )
% This function calculates a set of metrics from the input simulation data.
% The metrics are:
%    [ waveDisorder, blockRatio, blockCount, APDmean, APDstd, WLmean, WLstd, DF, F_std, multiActRatio, multiActCount ]
%
% If an empty value for 'simdata' is provided, then a set of NaN values is 
% output as the metrics (but with the same number of values)

if ~isempty(simdata)

    % Wave disorder
    wave_disorder = calcWaveDisorder( simdata );

    % Sensitivity to Block
    [block_sensitivity, block_count] = calculateNeighbourMetrics( simdata );
    
    % APD metrics - first activation event
    APD_mean = mean( simdata.APDs( simdata.APDs(:,1) > 0 ), 1 );
    APD_std = std( simdata.APDs( simdata.APDs(:,1) > 0 ), 1 );
    
    % Wavelength metrics
    [~, WL_mean, WL_std] = calcWavelengthMap( simdata );
    
    % Dominant Frequency
    [DF, F_std] = calcDominantFreq( simdata );
    
    % Activated sites
    multi_act_ratio = sum( simdata.N_activations > 1 ) / sum( simdata.N_activations > 0);
    multi_act_count = sum( simdata.N_activations > 1 ) / length( simdata.N_activations);
            
    % Store each of these metrics
    metrics = [wave_disorder, block_sensitivity, block_count, APD_mean, APD_std, WL_mean, WL_std, DF, F_std, multi_act_ratio, multi_act_count];
    
else
    
    % Null data
    metrics = nan(1,11);
    
end