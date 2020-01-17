function main_effects = calcMainEffects(emulator)
% This function calculates the main effects for the different parameters,
% as specified by the 'ranges' input. Main effects are calculated by
% evaluating the GP's predictions at a very large number of points, thus
% averaging over all other parameters to get a function defined as the
% 'mean' value of the output(s) for each value of the individual parameter

the 'main' effect of each
% individual parameter

nan_handling = 'zero_nans';

% Define ranges for the parameters

% Full range
ranges = [ 0.00, 0.00, 0.00;       1.25, 1.25, 1.25];
% Short range 
%ranges = [ 0.525, 0.525, 0.525;   0.725, 0.725, 0.725];

% Number of samples to use in Monte Carlo integration
N_samples = 50000;
% Number of points used to define the main effect function
N_curvepoints = 100;

% Read out number of parameters
N_params = size(ranges,2);

% Initialise some matrices for memory benefits
predict_vals = zeros(N_samples, N_params);
main_effects.param_vals = zeros(N_params, N_curvepoints);
main_effects.output_vals = cell(N_params,1);


% Loop over parameters
for k = 1:N_params
    
    fprintf('Calculating mean effects for parameter %g, out of %g \n', k, N_params);
    
    % Create a vector of the other parameters
    others = 1:N_params;
    others(k) = [];
    
    % Loop over values for the parameter in question
    main_effects.param_vals(k,:) = linspace( ranges(1,k), ranges(2,k), N_curvepoints );
    for j = 1:N_curvepoints
        
        % Read out value of the parameter here
        val = main_effects.param_vals(k,j);
       
        % Sample other points randomly using an LHS design
        LHS = lhsdesign(N_samples, N_params-1,'criterion','none');
        LHS_vals = ranges(1,others) + LHS .* ( ranges(2,others) - ranges(1,others) );
        
        % Create a matrix of points to predict at
        predict_vals(:,others) = LHS_vals;
        predict_vals(:,k) = val;
        
        % Predict at all of these points using the partitioned emulator
        % prediction routine
        [GP_preds, ~] = PartitionedEmulatorPrediction(predict_vals, emulator);
        
        % Average over all the calculated values to get the main effect
        % value here. Exact behaviour depends on which NaN value handling
        % option was selected
        switch nan_handling
            case 'remove_nans'
                GP_preds_noNaNs = GP_preds; GP_preds_noNaNs(isnan(GP_preds_noNaNs)) = 0;
                main_effects.output_vals{k}(j,:) = sum(GP_preds_noNaNs) ./ sum(~isnan(GP_preds));
            case 'zero_nans'
                GP_preds( isnan(GP_preds) ) = 0;
                main_effects.output_vals{k}(j,:) = sum(GP_preds) / N_samples;
        end
        
    end

end

% Save the data
save('main_effect_data.mat','main_effects');