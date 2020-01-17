function FIGURE_plotData(data)

% Define parameter names
param_names = {'Hypoxia','Hyperkalemia','Acidosis'};

% Define fibrosis densities
fib_densities = [0.33, 0.36, 0.39];

% Range of parameters
range = 1.25;

% Specify the plot setup stuff - main figure
margin = 0.15;
xgap = -0.25;
ygap = 0.175;
titleSpace = -0.1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over the full dataset, just reading out which rows of the dataset
% were actually simulated in 2D or not
simulated = false(length(data.flags),1);
for k = 1:length(data.flags)
   
    % Check if the simulation data is present for this row
    if ~isempty(data.sims{k})
        simulated(k) = true;
    end
    
end

% Grab out the different classes of simulation that will be visualised
fail1D = ( data.flags == -1 & ~simulated );
fail2D = ( data.flags == -1 & simulated );
propagates = ( data.flags == 0 );
local_reentry = ( data.flags == 1 );
reentry_temp = ( data.flags == 2 );
reentry = ( data.flags == 3 );

% Read out dimensions
N_params = size(data.Rparams,2);
N_rows = length(fib_densities);
N_cols = 1/2 * N_params * (N_params - 1);

% Initialise the figure
figure('Units','Normalized','OuterPosition',[0 0 1 1]);

% Set up the axes for the figure
dx = (1 - 2*margin - (N_cols-1)*xgap ) / N_cols;
dy = (1 - 2*margin - N_rows*titleSpace - (N_rows-1)*ygap ) / N_rows;
for m = 1:N_cols
    for k = 1:N_rows
        xpos = margin + (m-1)*(dx+xgap);
        ypos = 1 - margin - k*titleSpace - (k-1)*(dy+ygap) - dy;
        ax{k,m} = axes('Position',[xpos, ypos, dx, dy]);
        hold( ax{k,m}, 'on' );
    end
end

% Add a tiny bit of shuffle to the parameter values (so they don't overlap)
data.Rparams = data.Rparams + 0.01 * ( rand( size(data.Rparams) ) - 0.5 );

%%% Plotting

% Loop over densities
for i = 1:length(fib_densities)
    
    % Grab out the subset of the data for this density
    fail1D_here = fail1D & (data.densities == fib_densities(i) );
    fail2D_here = fail2D & (data.densities == fib_densities(i) );
    propagates_here = propagates & (data.densities == fib_densities(i) );
    local_reentry_here = local_reentry & (data.densities == fib_densities(i) );
    reentry_temp_here = reentry_temp & (data.densities == fib_densities(i) );
    reentry_here = reentry & (data.densities == fib_densities(i) );
    
    % Loop over parameter combinations
    counter = 0;
    for j = 1:N_params-1
        for k = j+1:N_params
            
            counter = counter + 1;
            
            % Plot the different classes
            plot( ax{i,counter}, data.Rparams(fail1D_here,j), data.Rparams(fail1D_here,k), 'x', 'MarkerEdgeColor', [0.7 0.7 0.7], 'LineWidth', 2 );
            plot( ax{i,counter}, data.Rparams(fail2D_here,j), data.Rparams(fail2D_here,k), 'x', 'MarkerEdgeColor', [0.35 0.35 0.35], 'LineWidth', 2 );
            plot( ax{i,counter}, data.Rparams(propagates_here,j), data.Rparams(propagates_here,k), '.', 'MarkerEdgeColor', [0.20 0.20 0.95], 'MarkerSize',10 );
            plot( ax{i,counter}, data.Rparams(local_reentry_here,j), data.Rparams(local_reentry_here,k), 'd', 'MarkerEdgeColor', [0 0.8 0], 'MarkerFaceColor', [0 0.8 0], 'MarkerSize',5.5 );
            plot( ax{i,counter}, data.Rparams(reentry_here,j), data.Rparams(reentry_here,k), 'd', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize',7 );
            plot( ax{i,counter}, data.Rparams(reentry_temp_here,j), data.Rparams(reentry_temp_here,k), 'd', 'MarkerEdgeColor', [1 0 1], 'MarkerFaceColor', [1 0 1], 'MarkerSize',7);
            
            % Ensure axes are square
            axis( ax{i,counter}, 'square' );
            
            % Set limits
            axis( ax{i,counter}, [0 range 0 range] );
            
            % Increase fontsize
            set( ax{i,counter}, 'FontSize', 20 );
            % Turn on axis box
            box( ax{i,counter}, 'on' );
            
            % Turn off Xticks if this is not a bottommost axis
            if i < length(fib_densities)
                ax{i,counter}.XTick = [];
            end
            % Turn off Yticks if this is not a leftmost axis
            if counter > 1
                ax{i,counter}.YTick = [];
            end
            
            % If this is the bottom row, add X and Y labels
            if i == length(fib_densities)
                xlabel( ax{i,counter}, param_names{j}, 'FontSize', 20);
                ylabel( ax{i,counter}, param_names{k}, 'FontSize', 20);
            end
            
            % Append density text to the rightmost column
            if j == N_params-1 && k == N_params   % Check for final parameter combination (rightmost column)
                text( ax{i,counter}, range*1.10, range/2, ['$\phi = ',num2str(fib_densities(i)),'$'], 'Interpreter','LaTeX','FontSize', 24);
            end
            
        end
    end
    
    
end