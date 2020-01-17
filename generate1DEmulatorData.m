function generate1DEmulatorData(data1D)
% This function loads in the '1D' data, and attempts to build emulators to
% predict the base quantities of interest (APD, CV) as well as predict the
% level of susceptibility to conduction block

% Specify the number of training points to use (remaining will be used for
% testing)
N_train = 800;

% Specify which metrics are to be emulated
use_metrics = [1, 2, 3, 4, 6, 12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load in the data if it wasn't provided
if nargin < 1
    load('data1Dmat.mat','data1D');
end

% Define the inputs and outputs for the emulator
inputs = data1D.Rparams;
outputs = data1D.metrics(:,use_metrics);

% Define discrete outputs to predict
discrete = data1D.metrics(:,7);

% Define the class used for partitioning
classes = data1D.metrics(:,7);

% Dimensions of data
[N_data, N_outputs] = size(outputs);

% Set a random seed so that random selections of train/test data are
% consistent
rng(18);
order = randperm(N_data);
inputs = inputs(order,:);
outputs = outputs(order,:);
classes = classes(order,:);
discrete = discrete(order,:);

% Split into training and testing data
train_data.X = inputs(1:N_train,:);
train_data.Y = outputs(1:N_train,:);
train_data.C = classes(1:N_train,:);
train_data.D = discrete(1:N_train,:);
test_data.X = inputs(N_train+1:end,:);
test_data.Y = outputs(N_train+1:end,:);
test_data.C = classes(N_train+1:end,:);
test_data.D = discrete(N_train+1:end,:);

% Create the partitioned emulator (all outputs at once)
emulator = CreatePartitionedEmulator(train_data, test_data, 'predict');

% Predict the values at the test points - use inputs because test_data was
% potentially modified by removing variables for classification, etc.
% This data is added to the structures for the train/test data
[train_data.Ypred, train_data.Cpred] = PartitionedEmulatorPrediction(train_data.X, emulator);
[test_data.Ypred, test_data.Cpred] = PartitionedEmulatorPrediction(test_data.X, emulator);

% Gather all the data into a single structure
emulator_data.train_data = train_data;
emulator_data.test_data = test_data;
emulator_data.emulator = emulator;


% Save the emulator data
save('emulator1Ddata.mat', 'emulator_data');


end

