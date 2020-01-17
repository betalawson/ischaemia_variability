function FIGURE_1DSlices(emulator_data)
% This function plots slices of two quantities of interest across the
% parameter space

% Specify the plot setup stuff - main figure
margin = 0.1;
xgap = -0.35;
ygap = 0.125;
titleSpace = 0.05;

% Specify which metrics to plot (emulated output numbers)
plot_metrics = [3, 6];

% Specify ranges for the varied parameters
ranges = [0, 0;
          1.25, 1.25];
% Define stepsize
step = 0.005;

% Specify values of the third parameter
slice_values = [0.15, 1.15];

% List the names of these outputs of interest
output_names = {'Wavelength (mm)','Block Susceptibility'};

% List the names of the parameters
param_names = {'Hypoxia', 'Hyperkalemia'};
slice_param_name = 'Acidosis';

% Define colours for the different parameters
colours = [0.30 0.40 0.90;    % Blue
           0.60 0.25 0.60;    % Deep purple
           1.00 0.40 0.40];   % Red
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the emulator data if it wasn't provided
if nargin < 1
    load('emulator1Ddata.mat','emulator_data');
end
    
% Load in extra colormaps
load('extra_colormaps.mat','plasma','viridis');

% Set up colormaps for the figures
clrmap1 = [0.65, 0.65, 0.65; viridis(50:end,:)];
clrmap2 = plasma;

% Read out dimensions
N_cols = length(slice_values);
N_rows = length(plot_metrics);


% Initialise the figure
figure('Units','Normalized','OuterPosition',[0 0 1 1]);

% Set up the axes for the figure
dx = (1 - 2*margin - (N_cols-1)*xgap ) / N_cols;
dy = (1 - 2*margin - N_rows*titleSpace - (N_rows-1)*ygap ) / N_rows;
for i = 1:N_cols
    for j = 1:N_rows
        xpos = margin + (i-1)*(dx+xgap);
        ypos = 1 - margin - j*titleSpace - (j-1)*(dy+ygap) - dy;
        ax{i,j} = axes('Position',[xpos, ypos, dx, dy]);
    end
end

% Plot the slices
for j = 1:N_rows
    
    % Initialise the maximum value for this quantity
    maxval = -1e10;
    
    for i = 1:N_cols
        
        % Set up coords
        [X,Y] = meshgrid( ranges(1,1):step:ranges(2,1), ranges(1,2):step:ranges(2,2) );
        points = [X(:), Y(:), slice_values(i) * ones(size(X(:)))];
        % Get emulator predictions
        Ypred = PartitionedEmulatorPrediction( points, emulator_data.emulator );
        % Set NaN values to large negative dummy value for plotting
        Ypred(isnan(Ypred)) = -100;
        % Predictions of block sensitivity are [0,1]
        Ypred( Ypred(:,6) < 0, 6 ) = 0;
        Ypred( Ypred(:,6) > 1, 6 ) = 1;
        
        % Visualise
        pcolor(ax{i,j},X,Y, reshape(Ypred(:,plot_metrics(j)),size(X,1),size(X,2) ) );
        shading(ax{i,j},'flat');
        
        % Specify colormap
        if j == 1
            colormap(ax{i,j}, clrmap1);
        else
            colormap(ax{i,j}, clrmap2);
        end
        
        % Colormap and color axis
        clrbar = colorbar(ax{i,j});
        if i ~= N_cols
            set(clrbar,'Visible','Off');
        end
        
        % Update maximum value based on the predicted values
        maxval = max([Ypred(:,plot_metrics(j)); maxval]);
                
        % Labels, title, etc
        xlabel(ax{i,j}, param_names{1},'Fontsize',20);
        ylabel(ax{i,j}, param_names{2},'Fontsize',20);
        title(ax{i,j},  {output_names{j},[slice_param_name,' = ',num2str(slice_values(i))]}, 'FontSize', 18);
        axis(ax{i,j},'square');
        box(ax{i,j},'off');
        hold(ax{i,j},'on');
        set(ax{i,j}, 'FontSize', 20);
        
    end
    
    % Only set colour axis after both columns have been processed, so that
    % a consistent scale can be assigned to each using the "maxval" that
    % encompasses both sets of predictions
    for i = 1:N_cols
        caxis(ax{i,j}, [ -0.01*maxval , maxval]);
    end
    
end