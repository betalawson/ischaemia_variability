function [SVM_learner_template, SVM_hyperparameter_options] = prepareSVM(SVM_options)

% This function simply sets up the SVM-building code with the parameters
% that MATLAB requires. It may be edited by the user (refer to templateSVM
% and fitcecoc in MATLAB help)
SVM_hyperparameter_options.Verbose = 0;
SVM_hyperparameter_options.ShowPlots = 0;
SVM_hyperparameter_options.AcquisitionFunctionName = 'expected-improvement-plus';
SVM_hyperparameter_options.MaxObjectiveEvaluations = SVM_options.max_func_evals;
SVM_learner_template = templateSVM('Standardize',1,'KernelFunction',SVM_options.kernel_func);

end

