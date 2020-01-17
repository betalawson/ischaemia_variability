function [code_options, SVM_options, GP_options, isomap_options] = DEFINEOptions

% This is the function used by the user to change options of the relevant
% codes

% Specify whether or not the points which were not correctly classified
% contribute to the training data for emulators or not
code_options.remove_misclassified_from_training = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the kernel function to use for SVM
SVM_options.kernel_func = 'gaussian';

% Specify the number of function evaluations to use in Bayesian
% optimisation of SVM hyperparameters, and how many times to repeat this
% process (best is taken)
SVM_options.hyperparam_retries = 1;
SVM_options.max_func_evals = 30;



%%%%%%%%%%%%%%%%%%%% GAUSSIAN PROCESS (GP) PARAMETERS %%%%%%%%%%%%%%%%%%%%

% Mean function for the GP, fit using regression
% Options are:  'none', 'constant', 'linear'
GP_options.basis = 'linear';

% Covariance function for the GP. Automatic Relevance Determination (ARD) 
% is currently always performed
% Options are:  'SquareExp', 'Matern32', 'Matern52'
GP_options.covType = 'Matern32';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

