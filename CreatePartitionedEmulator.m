function [PartitionedEmulator] = CreatePartitionedEmulator(train_data, test_data, classification_function)

% This function creates an emulator by first partitioning the parameter
% space using an SVM classifier. Each classification is then emulated
% separately by its own Gaussian process (GP). The output is a single
% struct, PartitionedEmulator, which features a classifier model, GP
% emulators for each classification, and information regarding accuracy
% of both the classification, and the emulation itself (against test data)


%%%%% SETUP

% Load in the options
[code_options, SVM_options, GP_options] = DEFINEOptions;

% Extract contents of provided data
train_X = train_data.X;
train_Y = train_data.Y;
train_C = train_data.C;

test_X = test_data.X;
test_Y = test_data.Y;
test_C = test_data.C;

% Read out data dimensions
[N_outputs] = size(train_Y, 2);

%%%%% BUILD CLASSIFIER IF NEEDED
if strcmp(classification_function,'predict')

    fprintf('Constructing SVM classifier... ');
    
    % Build an SVM classifier from the training data. Actually, we try a
    % number of repeats (due to stochastic nature of Bayesian optimisation) and
    % take the best from those found
    [SVM_learner_template, SVM_hyperparameter_options] = prepareSVM(SVM_options);
    for i = 1:SVM_options.hyperparam_retries
        
        % Build an SVM classifier with Bayesian optimisation of hyperparameters
        SVM{i} = fitcecoc(train_X, train_C, 'Learners', SVM_learner_template, 'Coding', 'onevsone', 'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', SVM_hyperparameter_options);
        
        % Apply this classifier to the test data, and store its performance
        pred_C_test = predict(SVM{i}, test_X);
        performance_ratings(i) = sum( test_C == pred_C_test ) / length(test_C);
        
    end
    
    % Now, select the best found SVM as the one to use
    [best_SVM_performance, loc] = max(performance_ratings);
    PartitionedEmulator.classifierModel = SVM{loc};
    
    % Store classification performance to be parsed out of the function
    PartitionedEmulator.classification_performance = best_SVM_performance;
    
    fprintf('Done! \n');
    
    % Work out what the predicted classes are for the training data
    pred_C_train = predict(PartitionedEmulator.classifierModel, train_X);
    
else
    
    % If not predicting, then classifications are guaranteed to be right
    pred_C_test = test_C;
    
end


% Each emulator is created using separate training data, so we loop over
% the emulators to be built separately
for i = 0:max(train_C)
    
    fprintf('Class %g: ',i);
    
    % Pick out the training data from the full set
    if code_options.remove_misclassified_from_training && strcmp(classification_function,'predict')
        train_X_by_class{i+1} = train_X(train_C == i & pred_C_train == i,:);
        train_Y_by_class{i+1} = train_Y(train_C == i & pred_C_train == i,:);
    else
        train_X_by_class{i+1} = train_X(train_C == i,:);
        train_Y_by_class{i+1} = train_Y(train_C == i,:);
    end
           
    % Construct an emulator for this classification
    PartitionedEmulator.emulatorsByClass{i+1} = buildEmulator(train_X_by_class{i+1}, train_Y_by_class{i+1}, GP_options);
    
    % Now, gather the test data that will use this emulator
    if ~strcmp(classification_function,'predict')
        pred_C_test = test_C;
    end
        
    test_X_this_class = test_X(pred_C_test == i, :);
    
    % Use this emulator to predict values for those points
    pred_Y_this_class = EmulatorPrediction(test_X_this_class, PartitionedEmulator.emulatorsByClass{i+1});
    
    % Store these back in the full list of predicted biomarkers
    pred_Y(pred_C_test == i,:) = pred_Y_this_class;
    
end

fprintf('Emulator training complete! \n');

%%%%% PERFORMANCE

fprintf('Calculating performance... ');

% To calculate emulator performance, we take a simple combined mean square 
% error, but weighted by the range of values of each biomarker
PartitionedEmulator.MSE = 0;
PartitionedEmulator.noMissclassMSE = 0;
for i = 1:N_outputs
    
    % Find the range of this biomarker
    Y_range = max(test_Y(:,i)) - min(test_Y(:,i));
    
    % Add up the MSE, both for the case where misclassified data are
    % included or not
    PartitionedEmulator.noMissclassMSE = PartitionedEmulator.noMissclassMSE + sum((test_Y(pred_C_test == test_C,i) - pred_Y(pred_C_test == test_C,i) ).^2) / Y_range / sum(pred_C_test == test_C);
    PartitionedEmulator.MSE = PartitionedEmulator.MSE + sum((test_Y(:,i) - pred_Y(:,i) ).^2) / Y_range / size(test_Y,1);

end

% Append the classification function to the emulator
PartitionedEmulator.classification_function = classification_function;

fprintf('Done! \n');