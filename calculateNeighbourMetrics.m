function [block_sensitivity, block_count] = calculateNeighbourMetrics( simdata )
% This function calculates a measure of sensitivity to conduction block by
% checking how many neighbours there are with mismatched first activation
% times

% Define the threshold used for assigning the activation at two
% neighbouring sites as a 'mismatch' (evidence of block)
mismatch_threshold = 5;    % (in ms)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate a list of all connections between sites (this will double count
% them)
[A,~] = Occ2Adjacency(simdata.occ);
[j,i,~] = find(logical(A));

% Create a mapping that converts the elements of A into the correct
% positions in the full data structure
empty_locs = find(~simdata.occ);

% Read out only the first activation time for all locations
activation_times = simdata.activation_times(:,1);

% Read out the APDs for these activations
APDs = simdata.APDs(:,1);

% Re-order the activation times and APDs so they correspond to the ordering
% induced by A
[~,index] = sortrows(simdata.mapping,[1 2]);
activation_times = activation_times(index);
APDs = APDs(index);

% Count the number of connections (graph edges)
N_edges = length(i);

% Initialise the RVI values to a very high value everywhere (low values
% represent high risk)
RVI = 1e10 * ones( size(activation_times) );

% Loop over all edges, checking if their activation times are mismatched
block_count = 0;
valid_count = 0;
for k = 1:N_edges
    
    % Read out sites represented by this edge
    here = empty_locs( i(k) );
    neigh = empty_locs( j(k) );
        
    % Read out the activation times for the site and its neighbour
    AT_here = activation_times( here );
    AT_neigh = activation_times( neigh );
    RT_here = AT_here + APDs( here );
    
    % For counting of mismatches, only check each edge once by checking its
    % indices i and j
    if i(k) < j(k)

        % Mismatch is a difference in activation time over the threshold
        if abs( AT_here - AT_neigh ) > mismatch_threshold
            block_count = block_count + 1;
        end
        % In calculating a ratio, only count those connections where at least
        % one site was activated
        if ~( AT_here == -1 && AT_neigh == -1 )
            valid_count = valid_count + 1;
        end
    
    end
    
    % Calculate the re-entry vulnerability index (Orini et al. 2019) by
    % comparing the repolarisation time at the current site to the
    % activation times of its neighbours (minimum is used, so simply check
    % if this edge represents an update to that minimal value)
    
    % Only calculate where both sites were activated
    if AT_here ~= -1 && AT_neigh ~= -1
        RVI( index(here) ) = min( [ RVI( index(here) );  RT_here - AT_neigh ] );
    end
    
end

% Block sensitivity is the proportion of edges that showed a mismatch
block_sensitivity = block_count / valid_count;

% Set any RVI values that weren't changed from the default to NaN
RVI( RVI == 1e10 ) = NaN;