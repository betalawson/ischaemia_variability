function emulator = buildEmulator(X, Y, GP_options)

% This function takes supplied training and test data and constructs an
% emulator that can then be used to make predictions at other locations in
% the parameter space. The inputs are the training data (inputs X and
% outputs Y), and options that define what to do with the Gaussian
% process regression. This version of the file simply calls MATLAB's built
% in Gaussian process fitting routine.

      
% Read out the desired covariance function and express it in
% MATLAB's desired format
switch GP_options.covType
    case 'Matern32'
        covFun = 'ardmatern32';
    case 'Matern52'
        covFun = 'ardmatern52';
    case 'SquareExp'
        covFun = 'ardsquaredexponential';
end        
        
% Each feature is emulated separately, but only if there are values
% for this output recorded for this classification. Otherwise, a
% 'dummy' GP is created
for i = 1:d
    fprintf('%g... ',i);
    if ~all(isnan(Y(:,i)))
        Y_GPs{i} = fitrgp(X,Y(:,i),'BasisFunction',GP_options.basis,'KernelFunction',covFun,'FitMethod','sr','PredictMethod','exact','Standardize',1);
    else
        Y_GPs{i} = 'dummy';
    end
end
fprintf('Done! \n');
      
%%% STORE ALL REQUIRED VARIABLES IN THE emulator STRUCT
emulator.covFun = covFun;
emulator.d = d;
emulator.Y_GPs = Y_GPs;
emulator.X = X;