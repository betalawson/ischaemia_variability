function FIGURE_MainEffectPlots(main_effects)
% This function plots the main effects of the three parameters, for each of
% the quantities of interest


% Specify the plot setup stuff - main figure
margin = 0.1;
xgap = 0.075;
ygap = 0.175;
titleSpace = 0.025;
rightSpace = 0.05;

% Specify the number of rows to use
N_rows = 2;

% Specify which metrics to plot in the main figure
plot_outputs = 1:size(main_effects.output_vals{1},2);

% List the names of these outputs of interest
output_names = {'APD (ms)', 'Conduction Velocity (mm/ms)', 'Wavelength (mm)', 'Rest Potential (mV)', 'Action Potential Amplitude (mV)', 'Block Susceptibility'};

% List the names of the parameters
param_names = {'Hypoxia', 'Hyperkalemia', 'Acidosis'};

% Define colours for the different parameters
colours = [0.30 0.40 0.90;    % Blue
           0.60 0.25 0.60;    % Deep purple
           1.00 0.40 0.40];   % Red


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate the number of columns based on the number of rows
N_cols = ceil(length(plot_outputs) / N_rows);

% Count the number of parameters
N_params = length(main_effects.output_vals);

% Initialise the figure
figure('Units','Normalized','OuterPosition',[0 0 1 1]);

% Set up the axes for the main figure
dx = (1 - 2*margin - (N_cols-1)*xgap - rightSpace) / N_cols;
dy = (1 - 2*margin - N_rows*titleSpace - (N_rows-1)*ygap ) / N_rows;
for i = 1:N_cols
    for j = 1:N_rows
        xpos = margin + (i-1)*(dx+xgap);
        ypos = 1 - margin - j*titleSpace - (j-1)*(dy+ygap) - dy;
        ax{i,j} = axes('Position',[xpos, ypos, dx, dy]);
    end
end

for m = 1:length(plot_outputs)
    
    % Find which row and column this is
    col = m - floor((m-1)/N_cols)*N_cols;
    row = ceil(m/N_cols);
    
    % Hold figure
    hold(ax{col,row},'on');
    
    % Plot the individual main effect trends for each parameter
    for k = 1:N_params
        plot(ax{col,row}, main_effects.param_vals(k,:), main_effects.output_vals{k}(:,m), 'LineWidth',2,'Color',colours(k,:));
    end
    
    % Add title
    title(ax{col,row}, output_names{m}, 'Fontsize', 20);
    
    % Adjust font axis
    set(ax{col,row}, 'FontSize', 20);
    
    % Add axis labels
    xlabel(ax{col,row}, 'Extent of Ischaemia Component','FontSize',20);
    ylabel(ax{col,row}, 'Value','FontSize',20);
    
    % Set axis limits
    xlim(ax{col,row}, [ min(main_effects.param_vals(k,:)), max(main_effects.param_vals(k,:))]);
        
end

% Append legend
% Add a colorbar
axpos = get(ax{N_cols, 1},'Position');
legend_obj = legend(param_names,'Position', [axpos(1)+axpos(3)*4.2/4,  axpos(2)+axpos(4)*0.5, 0.1, 0.1], 'FontSize', 18);