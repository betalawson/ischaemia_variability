function [propagating, Vss, Sss, CV, trace] = findSteadyState(problem)
% This function simulates the monodomain equation on a mesh where the
% volume fraction (fraction of accessible material, as opposed to say
% collagenous occlusions due to fibrosis) and conductivity tensor are
% allowed to vary between different elements. Completely non-conductive
% elements may also be included in the mesh.
%
% A vertex-centred finite volume method is used to numerically integrate
% the monodomain equation, with bilinear interpolation allowing for
% integration of non-diagonal conductivity tensors, as well as non-local
% integration of the other terms.
%
% Input is the problem structure, as created by one of the createProblem
% type files.

% Monodomain parameters
chi = 1400;                               % Surface-to-volume ratio for tissue (cm^-1)
Cm = 1;                                   % Tissue capacitance per unit area (uF/cm²) - cell capacitance is defined in ionic model files

% Stimulus settings
stim_dur = 2;                             % Stimulus duration (ms)
stim_amp = 38;                            % Amplitude of stimulus per unit area (uA/cm²)
stim_delay = 1000;                        % Delay before initiating stimulus
stim_BCL = 1000;                          % Cycle length (in ms) with which to trigger stimuli


% Timestepping and solution methods
dt = 0.02;                                % Timestep (ms)
solve_exact = 0;                          % Require exact solves (direct methods) for the linear systems that result from the time and space discretisations
second_order = 0;                         % Uses second order timestepping. Threatens stability, but provides better accuracy for sufficiently low timestep

% Plotting
visualise = 0;                            % Flag for whether to visualise or not
save_anim = 0;                            % Flag for whether or not to save an animation (filename same as problem name, CAREFUL not to overwite!)
plot_interval = 5;                        % Time interval for plotting (ms)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the scale factor for diffusion (in monodomain model), denoted
% alpha to match paper
alpha = 1 / chi / Cm;

% Perform the encoding of the problem (stiffness and mass matrices) that
% only need to be calculated once. This code also outputs some extra
% calculated mesh information, including a list of which nodes are active
[K, M, mesh] = encodeProblem(problem.occ_map, problem.D_tensor, problem.Vfrac, problem.grid, alpha);
M = speye(size(M));

% Finish preparing the numerical method in terms of these matrices
[A_new, A_old, A_J] = prepareNumerics(K, M, dt, second_order);

% Read out the list of active nodes from mesh file for notational
% cleanliness
active = mesh.active;

% Read out the number of nodes to solve at
N = length(active);

% Read out stimulus sites, and convert to a vectorised form
stim_sites1 = problem.stim_sites1;
stim_sites1 = stim_sites1';
stim_sites1 = stim_sites1(:);
stim_sites2 = problem.stim_sites2;
stim_sites2 = stim_sites2';
stim_sites2 = stim_sites2(:);

% Read out cell models, and convert to a vectorised form
cell_models = problem.cell_models;
model_assignments = problem.model_assignments;
model_assignments = model_assignments';
model_assignments = model_assignments(:);

% Read out node X and Y locations, and also vectorise them
nodeX = problem.nodeX;
nodeY = problem.nodeY;
nodeX = nodeX'; nodeX = nodeX(:);
nodeY = nodeY'; nodeY = nodeY(:);

% Check to see if the problem structure has extra parameters provided, if
% so grab them out
if isfield(problem,'extra_params')
    extra_params = problem.extra_params;
else
    extra_params = [];
end

% Create video object if one is needed
if save_anim && visualise
    vid_obj = VideoWriter('outputVideo.avi');
    open(vid_obj);
end

% Initialise problem
[V, S] = initialiseProblem(cell_models, model_assignments, active);

% The old value for the state variables is also set to the current value
% for the first step. Also, the information that comes from the cell model
% (gating variable rate constants and steady states) is initialised as
% blank to show there is no old information for these
S_old = S;
Sinf = [];
invtau = [];
I_stim_old = zeros(N,1);
J_old = [];


% Designate three nodes, one used for checking propagation (opposite end of
% fibre) one for checking steady state (close to stimulus site) and one for
% recording a trace
Lx = max(nodeX(:)); Ly = max(nodeY(:));
SS_check_node = find( abs( nodeX - 0.05 * Lx ) <= problem.grid.dx/2, 1);
propagation_check_node = find( abs( nodeX - 0.9 * Lx ) <= problem.grid.dx/2, 1);
node_dists = sqrt( (nodeX - Lx * 0.5).^2 + (nodeY - Ly * 0.5).^2 );
active_node_dists = node_dists(active);
[~, I] = sort(active_node_dists, 'ascend');
active_nodes = find(active);
trace_node = active_nodes( I(1) );


% Set up the wavespeed measurement regions
start_sites = find( nodeX(:) >= 0.35 * Lx );
end_sites = find( nodeX(:) >= 0.95 * Lx );
measure_dist = min( nodeX(end_sites) ) - min( nodeX(start_sites) );
stim_times2 = []; % Turn off other stimulus to be safe


% Initialise the flag for waves successfully propagating to true (will
% become false if it is found a wave fails to propagate
stim_propagated = 1;
propagating = 1;
start_reached = 0;
end_reached = 0;

% Initialise trace structure
trace.t = [];
trace.V = [];

% Loop over time integrations
t = 0;
steady_state_found = 0;
SS_check_vals = [];
while ~steady_state_found && propagating
    
    % Increment time
    t = t + dt;
    
    
    %%% REACTION STEP HANDLING
    
    % Stimulate if this is a stimulus time
    I_stim = zeros(N,1);
    
    check_time = t - stim_delay - floor( (t - stim_delay) / stim_BCL ) * stim_BCL;
    
    if check_time <= stim_dur && check_time >= 0 && t >= stim_delay
        I_stim(stim_sites1) = -stim_amp;
        stim_propagated = 0;  % Set flag for this stimulus to zero, so it can be re-set to one if the wavefront does successfully propagate
        start_reached = 0; end_reached = 0;   % Set flags for conduction velocity measurement to zero, will be reset and CV recalculated for each pulse
    end
    
    % Also use the same check time variable to see if propagation is 
    % succeeding
    % (check if the wave made it to the other end of the tissue at a time 90% of the way to next stimulus)
    if check_time >= 0.9 * stim_BCL && ~stim_propagated
        propagating = 0;
        Vss = NaN;
        Sss = NaN(1,size(S,2));
        CV = NaN;
    end
        
    
       
    % Process reaction update - uses current voltage values and current
    % state variable values, S. Only processes active sites
    if second_order
        [I_ion, S_new, Sinf, invtau] = processReaction(V, S, S_old, Sinf, invtau, dt, I_stim, I_stim_old, cell_models, model_assignments, mesh, second_order, extra_params);
    else
        [I_ion, S_new] = processReaction(V, S, S_old, Sinf, invtau, dt, I_stim, I_stim_old, cell_models, model_assignments, mesh, second_order, extra_params);
    end

    
    % Calculate the total 'current density'
    J = (1/Cm) * ( I_ion(active) + I_stim(active) );
    
    
    
    %%% PROCESS UPDATE
    V_active = takeTimestep(V(active), J, J_old, A_new, A_old, A_J, solve_exact, second_order);
    V(active) = V_active;
    
    % Now update the stored current and old values for different variables
    S_old = S;  
    S = S_new;
    I_stim_old = I_stim;
    J_old = J;
    
    
    %%% VOLTAGE TRACE
    trace.t = [ trace.t, t ];
    trace.V = [ trace.V, V(trace_node) ];
    
    
    % Check if a steady state has been reached by storing the voltage and
    % state variables at the check locations, and then comparing them to
    % the past values. This happens just before the stimulus begins
    if abs( check_time - ( stim_BCL - 2*dt ) ) <= dt*0.6 && t > stim_delay
       
        % Append the values for the steady state check
        SS_check_vals = [SS_check_vals; [V(SS_check_node), S(SS_check_node,:)] ];
        % Compare these to the previous time's state
        if size(SS_check_vals,1) > 1
            dCheck = abs( SS_check_vals(end,:) - SS_check_vals(end-1,:) ) ./ SS_check_vals(end-1,:);
            if all(dCheck < 0.01)
                steady_state_found = 1;
                Vss = V(SS_check_node);
                Sss = S(SS_check_node,:);
            end
        end
        
    end
    
    
    % Check if the node for propagation checking has been excited
    if ~stim_propagated
        if V(propagation_check_node) > -35
            stim_propagated = 1;
        end
    end
    
    
    %%% CONDUCTION VELOCITY CALCULATION
    
    % Check if any of the start sites has been stimulated
    if ~start_reached && any(V(start_sites) > -35)
        start_reached = 1;
        start_time = t;
    end
    
    % Check if any of the end sites has been stimulated
    if ~end_reached && any(V(end_sites) > -35)
        end_reached = 1;
        end_time = t;
        CV = measure_dist / (end_time - start_time);
    end
    
    
    
    % Check plot frequency, plot if hit
    if visualise && ( t - floor(t / plot_interval) * plot_interval <= dt )
        
        % Visualise the current state
        visualiseState( V, problem.Nx, problem.Ny, problem.occ_map, t );
        
        % Write frame to video object if animation requested
        if save_anim
            frame = getframe(gcf);
            writeVideo(vid_obj, frame);
            % Otherwise, draw on screen so user can watch in real time
        else
            drawnow;
        end
        
    end    
    
end

% Close video object if one was created
if save_anim && visualise
    close(vid_obj);
end