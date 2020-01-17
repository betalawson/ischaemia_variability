function visualiseSim(simdata, filename)
% This function takes the simulation info in input structure 'simdata' and
% plays a movie of the activation times. If a filename is provided, it also
% creates a movie in that name

sim_maxtime = 1000;

% Assume no video is to be written before checking
writevideo = 0;

if nargin > 1  % Create movie
    
    % Create videowriting object
    vid_obj = VideoWriter([filename,'.avi']);
    vid_obj.FrameRate = 200;
    open(vid_obj);
    writevideo = 1;
    
end
    

% Find the last time at which an activation occurred
max_t = max(simdata.activation_times(:) + simdata.APDs(:));
max_t = min([max_t, sim_maxtime]);

% Specify the timestep
dt = 1;

% Initialise the plot matrix
plot_matrix = -simdata.occ;

% Loop over moments in time, plotting those that sites that were active at
% that moment in time
N_activations = zeros( size(simdata.activation_times,1), 1 );
for k = dt:dt:max_t
    
    % Determine which locations were active at this moment in time
    [active_nodes,~] = find( abs(simdata.activation_times - k) < 0.6*dt );
    N_activations(active_nodes) = N_activations(active_nodes) + 1;
    
    % Add these to the plot matrix
    if ~isempty(active_nodes)
        %plot_matrix = indfill( linger*ones(length(active_nodes),1), simdata.mapping(active_nodes,:), plot_matrix);
        plot_matrix = indfill( simdata.APDs( ind2sub( size(simdata.APDs), active_nodes) ), simdata.mapping(active_nodes,:), plot_matrix);
    end
    
    % Visualise
    imagesc(flipud(plot_matrix));
    caxis([-1 1]);
    title(['Activated sites at time t = ',num2str(k)],'Fontsize',24);
    drawnow;
    
    if writevideo
        writeVideo(vid_obj,getframe(gcf));
    end
    
    % Subtract any values greater than zero (they will persist above zero
    % for a number of timesteps equal to linger
    plot_matrix(plot_matrix > 0) = plot_matrix(plot_matrix > 0) - dt;
    plot_matrix(plot_matrix > -1 & plot_matrix < 0) = 0;     % Set any values that have finished repolarising back to zero
    
    
end

if writevideo
    close(vid_obj);
end