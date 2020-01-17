function [pred_Y, pred_C] = PartitionedEmulatorPrediction(pred_X, PartitionedEmulator)

% This function takes a list of locations in the parameter space where
% predictions would like to be made, along with a partitioned emulator that
% is used to generate these predictions


% First, use the emulator's classifier to decide what each location in the
% parameter space should be classified as, and hence which emulator will be
% used to make predictions
if strcmp(PartitionedEmulator.classification_function,'predict')
    pred_C = predict(PartitionedEmulator.classifierModel, pred_X);
else
    [pred_C, pred_X] = feval(PartitionedEmulator.classification_function, pred_X);
end

% Now, loop over each classification, and use the associated emulator to
% make predictions for those points
for i = 0:max(pred_C)
    
    % Pick out the points which correspond to this classification
    this_class = (pred_C == i);
    
    % Only emulate if some data falls in this class
    if sum(this_class) > 0
        
        % Use the emulator associated with this classification to predict
        pred_Y(this_class,:) = EmulatorPrediction(pred_X(this_class,:), PartitionedEmulator.emulatorsByClass{i+1});
        
    end
    
end

