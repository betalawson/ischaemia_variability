function FIGURE_ReentryExample(data)    
% This function plots two example simulations from the full dataset in 2D.
% One re-entry, and one that fails to form a re-entry

% Specify whether to append text on the left
append_text = 0;

% Specify here the simulation numbers to use
sims = [597;     % No re-entry
        597;
        2613;
        2613];    % Re-entry

% Specify here which moments in time to use as snapshots for each
% simulation
plot_times = [75, 150, 225, 300;
              575, 600, 625, 650;
              75, 150, 225, 300;
              325, 350, 375, 625];
          
% Define the sizes for axis positions
margin = 0.05;           % Margin around all sides of plot
xsep = -0.175;             % Separation between axes in x direction
ysep = 0.015;             % Separation between axes in y direction
ygap = 0.05;              % Extra gap inserted to further separate different sims (currently inserted manually between rows 1 and 2)
leftTextSpace = 0.15;     % Space for text placed on left of image
rightTextSpace = 0;       % Space for text placed on right of image
titleSpace = 0.025;        % Space above each plot for title text
          
% Append text
append1 = {'Hypoxia      = 0.8283';
           'Hyperkalemia = 0.6742';
           'Acidosis     = 0.0599'};

append2 = {'Hypoxia      = 0.8098';
           'Hyperkalemia = 0.0253';
           'Acidosis     = 0.7024'};       
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the data if it hasn't been supplied
if nargin < 1
    load('fulldata.mat','data');
end

% Remove space for text if not using it
if ~append_text
    leftTextSpace = 0;
end

% Count the number of simulations and times for each sim
[N_rows, N_cols] = size(plot_times);

% Initialise Figure
figure('units','Normalized','OuterPosition',[0 0 0.7 1]);

% Determine the size of the individual axes
dx = ( 1 - leftTextSpace - rightTextSpace - 2*margin - (N_cols-1)*xsep ) / N_cols;
dy = ( 1 - titleSpace * N_rows - 2*margin - (N_rows-1)*ysep - ygap ) / N_rows;

% Create axes
ax = cell(N_rows,N_cols);
for i = 1:N_rows
    for j = 1:N_cols
        
        % Set up axis positions
        xpos = margin + leftTextSpace + (j-1)*(dx+xsep);
        ypos = 1 - margin - titleSpace - (i-1)*(dy+ysep+titleSpace) - dy - (i > 2)*ygap;
        % Create axis object
        ax{i,j} = axes('Position', [xpos, ypos, dx, dy]);
        
    end
end

% Plot all simulations
counter = 0;
for m = 1:length(sims)
    
    % Plot all snapshots for this simulation
    for k = 1:size(plot_times,2)
        
        % Increment counter
        counter = counter + 1;
    
        % Plot the snapshot
        visualiseSnapshot(ax{m,k},data.sims{sims(m)}, plot_times(m,k));
        % Title with the time
        title(ax{m,k},['t = ',num2str(plot_times(m,k))], 'FontSize', 16);
        
    end
end

% Append the text
if append_text
    text(ax{1,1}, -700, 450, append1, 'FontSize',18, 'FontName','Consolas');
    text(ax{3,1}, -700, 450, append2, 'FontSize',18, 'FontName','Consolas');
end

%%% Label the attempted re-entry

% Set hold on for the axis in question
hold(ax{2,3},'on');

% Set up the ellipse
x_cent = 140;
y_cent = 215;
x_rad = 40;
y_rad = 40;
t = 0:(2*pi/100):2*pi;
r = sqrt( 1 ./ ( 1 / (y_rad)^2 + ( 1 / (x_rad)^2 - 1 / (y_rad)^2 ) * cos(t) ) );
x = x_cent + r .* cos(t);
y = y_cent + r .* sin(t);

% Plot the ellipse on this axis
plot(ax{2,3}, x, y, 'b', 'LineWidth', 2);

end

function visualiseSnapshot(ax, simdata, time)
% This function takes the input simulation, and plots its current state
% (level set plot of voltage) for the requested time point


% Define colours
colours = [ 0.4, 0.4, 0.4;      % Dark grey - fibrosis
            0.5, 0.7, 1.0;      % Light blue - unactivated
            1.0, 1.0, 0.5;      % Yellow - activated first
            1.0, 0.3, 0.3];     % Red - subsequent activations
            

% Read out the simulation data of interest
AT = simdata.activation_times;
APD = simdata.APDs;

% Calculate which states are currently repolarising (activated but with no
% completed repolarisation yet)
%active = any( time <= (AT + APD) & time >= AT, 2);
[active_now_j, active_now_i] = find( time <= (AT + APD) & time >= AT );

% Active activations
%active_ATs = AT(active,:);
% Find which number activation this is (search along each row individually)
%for k = 1:size(active_ATs,1)
%    activation_nums = find( active_ATs(k,:)

% Visualise fibrosis in one colour, unactivated sites in another, and
% activated sites in the remaining colours (first or subsequent activation)
vis_mat = -double(simdata.occ);
vis_mat = indfill( active_now_i, simdata.mapping(active_now_j,:), vis_mat );
vis_mat(vis_mat > 2) = 2;  % Set all subsequent activations to be the "second" for plotting

% Plot and use the requested colouring
imagesc(ax, flipud( vis_mat ) );
colormap(colours);

% Ensure that even if there are no re-entries present, that the colour map
% isn't effected
caxis(ax,[-1 2]);

% Turn the axis off, also ensure it's square
axis(ax,'off','equal');


end

