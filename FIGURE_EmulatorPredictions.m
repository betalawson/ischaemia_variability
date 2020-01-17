function FIGURE_EmulatorPredictions(emulator_data)
% This function plots the predictions of the emulator passed in through the
% input "emulator_data", which is a structure that includes train and test
% data, as well as the trained emulator


% Specify the plot setup stuff - main figure
margin = 0.05;
xgap = 0.05;
ygap = 0.05;
titleSpace = 0.05;

% Specify the number of rows to use
N_rows = 2;

% Specify which metrics to plot in the main figure
plot_outputs = 1:size(emulator_data.train_data.Y,2);

% List the names of these outputs of interest
output_names = {'APD', 'Conduction Velocity', 'Wavelength', 'Rest Potential', 'Action Potential Amplitude', 'Block Susceptibility'};

% Calculate the number of columns based on the number of rows
N_cols = ceil(length(plot_outputs) / N_rows);

% Load the emulator data if it wasn't provided
if nargin < 1
    load('emulator1Ddata.mat','emulator_data');
end
    
% Grab out the individual parts of the provided data
test_data = emulator_data.test_data;


%%% EMULATOR VALIDATION

% Initialise the figure
figure('Units','Normalized','OuterPosition',[0 0 1 1]);

% Set up the axes for the main figure
dx = (1 - 2*margin - (N_cols-1)*xgap ) / N_cols;
dy = (1 - 2*margin - N_rows*titleSpace - (N_rows-1)*ygap ) / N_rows;
for i = 1:N_cols
    for j = 1:N_rows
        xpos = margin + (i-1)*(dx+xgap);
        ypos = 1 - margin - j*titleSpace - (j-1)*(dy+ygap) - dy;
        ax{i,j} = axes('Position',[xpos, ypos, dx, dy]);
    end
end

% Plot the check of emulator performance, predictions versus actual for the
% test data
for m = 1:length(plot_outputs)
   
    col = m - floor((m-1)/N_cols)*N_cols;
    row = ceil(m/N_cols);
    
    % Plot the data, and a line of equality
    hold(ax{col,row},'on');
    minval = min([test_data.Y(:,plot_outputs(m)); test_data.Ypred(:,plot_outputs(m))]);
    maxval = max([test_data.Y(:,plot_outputs(m)); test_data.Ypred(:,plot_outputs(m))]);
    ax_min = minval - 0.1*(maxval-minval);
    ax_max = maxval + 0.1*(maxval-minval);
    plot(ax{col,row}, [ax_min ax_max], [ax_min ax_max], 'k', 'LineWidth', 2);
    plot(ax{col,row}, test_data.Y(:,plot_outputs(m)), test_data.Ypred(:,plot_outputs(m)), 'b.', 'MarkerSize',12 );
    xlim(ax{col,row}, [ax_min ax_max]);
    ylim(ax{col,row}, [ax_min ax_max]);
    title(ax{col,row}, output_names{m}, 'Fontsize', 20);
    
end
