function FIGURE_dangerous_full(data, emulator_data)
% This function plots slices of two quantities of interest across the
% parameter space

% Specify the plot setup stuff - main figure
margin = 0.1;
xgap = -0.35;
ygap = 0.075;
titleSpace = 0;

% Specify which outputs of the emulator correspond to the quantities of
% interest
WL_output = 3;
block_output = 6;

% Specify ranges for the varied parameters
ranges = [0   , 0   , 0  ;
          1.25, 1.25, 1.25];

% Define stepsize
step = 0.025;

% Specify parameter combos
param_combos = [1, 2;
                1, 3;
                2, 3];
slice_params = [3; 2; 1];

% Define the densities
densities = [0.33, 0.36, 0.39];

% List the names of these outputs of interest
output_names = {'Wavelength (mm)','Block Susceptibility'};

% List the names of the parameters
param_names = {'Hypoxia', 'Hyperkalemia', 'Acidosis'};

% Define colours for the different parameters
colours = [0.35 0.35 0.35;    % Grey - non-propagating
    0.40 0.40 0.95;    % Blue - guaranteed safe
    0.40 1.00 0.40;    % White - ???
    0.90 0.30 0.30];   % Red - dangerous!

line_colours = [0.20 0.20 0.95;    % Blue - guaranteed safe
                0.5 0.5 0.5;    % Grey - non-propagating
                0.90 0.10 0.10];   % Red - dangerous!
            
line_alphas = [0.7, 0.3, 0.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the emulator data if it wasn't provided
if nargin < 1
    load('emulator1Ddata.mat','emulator_data');
end

% Loop over densities
for k = 1:length(densities)
    
    % Determine the cutoffs for the regions for this density
    overblock = max( data.metrics1D( data.densities == densities(k) & data.flags ~= -1 & ~isnan(data.flags), end) );
    dangerous = [ min( data.metrics1D( data.flags >= 2 & data.densities == densities(k), 3) ), max( data.metrics1D( data.flags >= 2  & data.densities == densities(k), 3) );
        min( data.metrics1D( data.flags >= 2 & data.densities == densities(k), end) ), max( data.metrics1D( data.flags >= 2  & data.densities == densities(k), end) )];
    
    % Create a plot showing the 3D parameter space as a heat map by
    % integrating over the third dimension. Also plot histograms along the
    % diagonal
    full_figs{k} = figure('Units','Normalized','OuterPosition',[0 0 1 1]);
    
    % Set up a full 3-D grid of points at which to evaluate the emulator
    x = ranges(1,1):step:ranges(2,1);
    y = ranges(1,2):step:ranges(2,2);
    z = ranges(1,3):step:ranges(2,3);
    
    [Y, X, Z] = meshgrid( x, y, z );
    points = [X(:), Y(:), Z(:)];
    
    % Evaluate the emulator's predictions at this point
    pred = PartitionedEmulatorPrediction( points, emulator_data.emulator );
    
    % Classify these predictions according to the "vulnerable window"
    pred_WL = pred(:,3);
    pred_block = pred(:,6);
    danger = pred_WL >= dangerous(1,1) & pred_WL <= dangerous(1,2) & pred_block >= dangerous(2,1) & pred_block <= dangerous(2,2);
    block = pred_block > overblock;
    safe = pred_block < dangerous(2,1);
    
    % Create 3-D arrays of which parameter values lead to which classes of
    % behaviour
    dangerC = reshape( double(danger), size(X) );
    blockC = reshape( double(block), size(X) );
    safeC = reshape( double(safe), size(X) );
    
    % Set up axes for the parameter plots
    N_rows = 3;
    N_cols = 3;
    dx = (1 - 2*margin - (N_cols-1)*xgap ) / N_cols;
    dy = (1 - 2*margin - N_rows*titleSpace - (N_rows-1)*ygap ) / N_rows;
    for m = 1:N_cols
        for n = 1:N_rows
            xpos = margin + (m-1)*(dx+xgap);
            ypos = 1 - margin - n*titleSpace - (n-1)*(dy+ygap) - dy;
            for i = 1:2
                pax{k,n,m,i} = axes('Position',[xpos, ypos, dx, dy]);
                hold( pax{k,n,m,i}, 'on');
                box(pax{k,n,m,i}, 'on');
            end
        end
    end
    
    % Plot the marginal densities of the three classes
    classes = {safe, block, danger};
    binheightmax = -1;    % Maximum bin height for plotting purposes
    for m = 1:3
       
        % Calculate kernel density estimates for the three 
        for n = 1:length(classes)
            h = histogram( pax{k,m,m,1}, points(classes{n},m), [ranges(1,m)-0.00000000001*(ranges(2,m)-ranges(1,m)):step*5:ranges(2,m)+0.00000000001*(ranges(2,m)-ranges(1,m))] );   % Slight shift to avoid weird bin lumping on edges
            h.FaceColor = line_colours(n,:);
            h.FaceAlpha = line_alphas(n);
            axis(pax{k,m,m,1},'square');
            
            % Update maximum bin height
            binheightmax = max([binheightmax; h.Values']);
            
        end
        
        % Set limits for the axis to the ranges
        xlim(pax{k,m,m,1},[ranges(1,m) ranges(2,m)]);
        
        % Set fontsize
        set(pax{k,m,m,1}, 'FontSize', 20);
        
        % Turn off the second axis for these
        pax{k,m,m,2}.Visible = 'Off';
        pax{k,m,m,2}.XTick = [];
        pax{k,m,m,2}.YTick = [];
        
    end
    
    for m = 1:size(param_combos,1)
        
        % Read out params
        param1 = param_combos(m,1);
        param2 = param_combos(m,2);
        paramS = slice_params(m);
        
        % Calculate extent of block/safety
        propensities = -squeeze(sum( safeC-blockC, paramS ) );   % Negative used because colormap starts at "safe"
        % Scale to make propensity [-1, 1]
        propensities = propensities / length(x);
        
        % Create a grid of values in 2D for these two parameters
        [X,Y] = meshgrid( ranges(1,param1):step:ranges(2,param1), ranges(1,param2):step:ranges(2,param2) );
        points = [];
        points(:,param1) = X(:);
        points(:,param2) = Y(:);
        
        % Scatter plot the propensity
        scatter( pax{k,param2, param1, 1}, Y(:), X(:), 10, propensities(:), 'filled' );
        
        % Use a colormap that varies from deep blue to white to dark grey
        colormap( pax{k,param2,param1,1}, [ [ (0.2:0.01:1)', (0.2:0.01:1)', ones( length( 0.2:0.01:1), 1) ];
            [ (1:-0.01:0.2)', (1:-0.01:0.2)', (1:-0.01:0.2)' ] ] );
        
        % Set limits to the parameter ranges
        xlim(pax{k,param2,param1,1},[ranges(1,param1) ranges(2,param1)]);
        ylim(pax{k,param2,param1,1},[ranges(1,param2) ranges(2,param2)]);
        
        % Color axis for [-1,1]
        caxis( pax{k,param2,param1,1}, [-1 1] );
        
        % Clean up plot visuals
        axis(pax{k,param2,param1, 1}, 'square');
        box(pax{k,param2, param1, 1}, 'off');
        set(pax{k,param2,param1,1}, 'FontSize', 20);
        
        % Turn off the now-unused axis
        axis( pax{k,param1,param2,1}, 'off' );
        pax{k,param1,param2,1}.Visible = 'off';
        
    end
    
    for m = 1:size(param_combos,1)
        
        % Read out params
        param1 = param_combos(m,1);
        param2 = param_combos(m,2);
        paramS = slice_params(m);
        
        % Calculate extent of block/safety
        propensities = squeeze(sum( dangerC, paramS ) );
        % Scale to make propensity [0, 1]
        propensities = propensities / length(x);
        propensities(abs(propensities) < 1e-6 ) = NaN;
        
        % Create a grid of values in 2D for these two parameters
        [X,Y] = meshgrid( ranges(1,param1):step:ranges(2,param1), ranges(1,param2):step:ranges(2,param2) );
        points = [];
        points(:,param1) = X(:);
        points(:,param2) = Y(:);
        
        % Scatter plot the propensity
        scatter( pax{k,param2, param1, 2}, Y(:), X(:), 10, -propensities(:), 'filled' );
        colormap(pax{k,param2,param1,2},'autumn');
        
        % Overlay the data points for re-entries
        plot( pax{k,param2,param1,2}, data.Rparams(data.flags==1&data.densities==densities(k),param1), data.Rparams(data.flags==1&data.densities==densities(k),param2), 'o', 'MarkerSize', 7, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[1 1 1], 'LineWidth',2);
        plot( pax{k,param2,param1,2}, data.Rparams(data.flags>=2&data.densities==densities(k),param1), data.Rparams(data.flags>=2&data.densities==densities(k),param2), 'o', 'MarkerSize', 7, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[1 0.5 0.5], 'LineWidth',2);
        
        % Set limits to the parameter ranges
        xlim(pax{k,param2,param1,2},[ranges(1,param1) ranges(2,param1)]);
        ylim(pax{k,param2,param1,2},[ranges(1,param2) ranges(2,param2)]);
        
        % Color axis for [0,1]
        caxis( pax{k,param2,param1,2}, [-1 0] );
        
        % Clean up plot visuals, including turning axis off altogether
        % (data is overlaid over above surface)
        axis(pax{k,param2,param1, 2}, 'square');
        box(pax{k,param2, param1, 2},'off');
        pax{k,param2,param1,2}.Visible = 'off';
        pax{k,param2,param1,2}.XTick = [];
        pax{k,param2,param1,2}.YTick = [];
        
        % Turn off the now-unused axis
        axis( pax{k,param1,param2,2}, 'off' );
        pax{k,param1,param2,2}.Visible = 'off';
                
    end
        
end

% Loop over densities again to do some minor extra modifications to the
% plots
for k = 1:length(densities)
    
    % Fix heights of histograms to be consistent
    for i = 1:3
        ylim(pax{k,i,i,1}, [0 binheightmax]);
    end
    
    % On left and bottom, append name
    for i = 1:3
        xlabel( pax{k,3,i,1}, param_names{i}, 'Fontsize', 20);
        ylabel( pax{k,i,1,1}, param_names{i}, 'Fontsize', 20);
    end
    
    % Turn off all Y labels for histograms
    for m = 1:3
        pax{k,m,m,1}.YTick = [];
    end
    
    % Add a label for the obstacle proportion
    axis(pax{k,1,3,1}, [0 1 0 1]);
    text(pax{k,1,3,1}, 0.5, 0.5, ['$\phi = $',num2str(densities(k))], 'Interpreter','LaTeX','Fontsize', 35, 'HorizontalAlignment', 'center', 'VerticalAlignment','middle');
    
end