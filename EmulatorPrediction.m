function pred_Y = EmulatorPrediction(pred_X, emulator)
% This function takes the supplied emulator, and uses it to make
% predictions at the specified points in the parameter space X_pred. This
% is a simplified version that works with emulators constructed using
% MATLAB's Gaussian process fitting.

% Read out the number of outputs
d = emulator.d;

%%% USE THE TRAINED GAUSSIAN PROCESSES TO EMULATE FEATURES
  
% Each feature is emulated separately. If the variable is not emulated for
% this class, 
for i = 1:d
    if ~strcmp(emulator.Y_GPs{i},'dummy')
        pred_Y(:,i) = predict(emulator.Y_GPs{i}, pred_X);
    else
        pred_Y(:,i) = nan( size(pred_X,1), 1);
    end
end

        
end