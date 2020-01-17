function FIGURE_1Dvs2D(data)
% This function plots figures demonstrating the correlation (or,
% relationship) between corresponding metrics in one dimension and two
% dimensions.


% Specify class colours
colours = [0.70 0.70 0.70;    % Light grey
    0.10 0.10 0.90;    % Dark Blue
    0.90 0.10 0.10];   % Dark Red

% Specify densities used
densities = [0.33, 0.36, 0.39];


% Load in the data if it hasn't been supplied
if nargin < 1
    load('fulldata.mat','data');
end

% First figure - Wavelength by Density
figure('units','Normalized','OuterPosition',[0 0 1 1]);
for k = 1:3
    subplot(1,3,k); hold on;
    % Plot the requested variables against one another, colourised according to simulation class (divide by 1000 for unit consistency)
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 0, 3), data.metrics2D( data.densities == densities(k) & data.flags == 0, 6)/1000, '.', 'MarkerSize', 16, 'MarkerEdgeColor', colours(1,:));
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 1, 3), data.metrics2D( data.densities == densities(k) & data.flags == 1, 6)/1000, '.', 'MarkerSize', 16, 'MarkerEdgeColor', [0 0.8 0]);
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 2, 3), data.metrics2D( data.densities == densities(k) & data.flags == 2, 6)/1000, '.', 'MarkerSize', 30, 'MarkerEdgeColor', colours(2,:));
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 3, 3), data.metrics2D( data.densities == densities(k) & data.flags == 3, 6)/1000, '.', 'MarkerSize', 30, 'MarkerEdgeColor', colours(3,:));
    % Adjust labelling, etc
    title([num2str(densities(k)*100),'% Obstacles'], 'FontSize', 22);
    xlabel('Wavelength (1D Fibre)','FontSize',20);
    ylabel('Wavelength (2D Fibrotic Tissue)','FontSize',20);
    set(gca,'FontSize',20);
    axis square;
end

% Second figure - Block Susceptibility by Density
figure('units','Normalized','OuterPosition',[0 0 1 1]);
for k = 1:3
    subplot(1,3,k); hold on;
       
    % Plot the requested variables against one another, colourised according to simulation class
    plot( data.metrics1D( data.densities == densities(k) & data.flags == -1, end), data.metrics2D( data.densities == densities(k) & data.flags == -1, 2), 'x', 'MarkerSize', 20, 'MarkerEdgeColor', [0,0,0]);
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 0, end), data.metrics2D( data.densities == densities(k) & data.flags == 0, 2), '.', 'MarkerSize', 16, 'MarkerEdgeColor', colours(1,:));
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 1, end), data.metrics2D( data.densities == densities(k) & data.flags == 1, 2), '.', 'MarkerSize', 16, 'MarkerEdgeColor', [0 0.8 0]);
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 2, end), data.metrics2D( data.densities == densities(k) & data.flags == 2, 2), '.', 'MarkerSize', 30, 'MarkerEdgeColor', colours(2,:));
    plot( data.metrics1D( data.densities == densities(k) & data.flags == 3, end), data.metrics2D( data.densities == densities(k) & data.flags == 3, 2), '.', 'MarkerSize', 30, 'MarkerEdgeColor', colours(3,:));
    
    % Adjust labelling, etc
    title([num2str(densities(k)*100),'% Obstacles'], 'FontSize', 22);
    xlabel('Block Susceptibility (1D Fibre)','FontSize',20);
    ylabel('Block % (2D Fibrotic Tissue)','FontSize',20);
    set(gca,'FontSize',20);
    axis square;
    set(gca,'yscale','log')
    ylim([5e-6 1e-1]);
    
end


end

