function createRunInfo(N_sets)

% Maximum number of operations to do in one parfor loop (just in case this helps with memory and such)
set_max = 100;

% Specify the parameters that control the healthy and ischaemic condition
% params:              [ [ATP]_i,  [K+]_o,  [K+]_i,  V_mod,  I_Na_factor,  I_CaL_factor,  I_NaCa_factor ] 
params_zero =  [   5.4  ,   5.4  ,   138.3,   0.0 ,     1.00    ,     1.00     ,      1.00      ];
params_one  =  [   2.6  ,   10.0 ,   125.0,   3.4 ,     0.75    ,     0.75     ,      1.00      ];    % Currently, I_NaCa is unaffected by acidosis

% Specify the ranges of the parameters (currently all parameters are
% designed around a base 0 to 1 range)
hypoxia_range = [0, 1.25];
hyperkalemia_range = [0, 1.25];
acidosis_range = [0, 1.25];

% Perform a Latin hypercube sampling over these ranges, for the requested
% number of parameter sets
X = lhsdesign(N_sets, 3);

% Convert these into values of the strength parameters using their ranges
hypoxia_strengths = hypoxia_range(1) + X(:,1) * (hypoxia_range(2) - hypoxia_range(1));
hyperkalemia_strengths = hyperkalemia_range(1) + X(:,2) * (hyperkalemia_range(2) - hyperkalemia_range(1));
acidosis_strengths = acidosis_range(1) + X(:,3) * (acidosis_range(2) - acidosis_range(1));

% Store these together for later convenience
ischaemia_strengths = [hypoxia_strengths, hyperkalemia_strengths, acidosis_strengths];

% Convert these strengths into the ischaemia parameters that take action in
% the model

% Hypoxia affects [ATP]_i
params(:,1) = params_zero(1) + hypoxia_strengths * ( params_one(1) - params_zero(1) );
% Hyperkalemia affects [K+]_o
params(:,2) = params_zero(2) + hyperkalemia_strengths * ( params_one(2) - params_zero(2) );
% Acidosis affects the remaining parameters ( [K+]_i, shifted reversal potential and current strengths for I_Na and I_CaL )
params(:,3:7) = params_zero(3:7) + acidosis_strengths * ( params_one(3:7) - params_zero(3:7) );    % Outer product, could also be .*

% Create the 1D fibre problem that will be used for steady state
% calculation
createFibreProblem('Fibre');
pause(2);
load('Fibre.mat','problem');

% Loop over each individual parameter set - these can be run independently
Vss = zeros(N_sets,1);
Sss = zeros(N_sets,8);      % Eight state variables for TT03   MIGHT BE ELEVEN
propagating = zeros(N_sets,1);

% Split up the problem into multiple smaller parfor loops for memory
% reasons
N_parsets = ceil(N_sets / set_max);

% Remind MATLAB that we loaded in this problem (funny parfor thing)
problem = problem;

% Enact loops
for m = 1:N_parsets
    
    % Start point
    start_pt = set_max*(m-1)+1;
    end_pt = min(N_sets, set_max*m);
    
    % Temporary variables
    prop_data = zeros(end_pt - start_pt + 1, 1);
    Vss_data = zeros(end_pt - start_pt + 1, 1);
    Sss_data = zeros(end_pt - start_pt + 1, 8);
    
    params_here = params(start_pt:end_pt,:);
    
    parfor k = 1:(end_pt-start_pt+1)
        % Pace the 1D Fibre at 1Hz until steady state, with this set of
        % ischaemia parameter values
        problem_here = problem;
        problem_here.extra_params = params_here(k,:);
        [prop_data(k), Vss_data(k), Sss_data(k,:)] = findSteadyState(problem_here);
    end
    
    % Store the data from the temporary variables (just keeping them out of
    % parfor)
    Vss(start_pt:end_pt) = Vss_data;
    Sss(start_pt:end_pt,:) = Sss_data;
    propagating(start_pt:end_pt) = prop_data;
    
end


% Now, manually append the steady state values for the extra state
% variables in Rafael's code (these are calculated directly from Vss)

% Append d_inf
Sss = [Sss, 1 ./ ( 1 + exp( -(Vss + 8)/7.5 ) )];
% Append r_inf (endo formulation)
Sss = [Sss, 1 ./ ( 1 + exp( -(Vss - 20)/6 ) )];
% Append xr2_inf
Sss = [Sss, 1 ./ ( 1 + exp( (Vss + 88)/24 ) )];

% Now write a text file that features the three parameter values, and then
% the voltage, and then the steady state
write_data = [ischaemia_strengths, params, propagating, Vss, Sss];

% Open file
file_obj = fopen(['param_SS_',num2str(N_sets),'configs.txt'],'w');

% Write parameter data to it. Writing is column-wise, so transpose the
% write data first
fprintf(file_obj, [ getFormatStr(size(ischaemia_strengths,2)), '| ', getFormatStr(size(params,2)), '| ', getFormatStr(1,'d'), '| ', getFormatStr(1), '| ', getFormatStr(size(Sss,2)), '\n'], write_data');

% Close file
fclose(file_obj);


function str = getFormatStr(N,type)

% Default type
if nargin < 2
    type = 'g';
end

% Initialise string
str = '';

% Add to the string the number of requested values
for k = 1:N
    str = [str,'%',type,' '];
end