function A = indfill(vals, inds, varargin)
% This function fills in the provided values (vector) into the locations
% indexed by 'inds'. The user can optionally provide a matrix as an 
% additional argument, if so its size is used for attempting to store the 
% provided data. The values of the input matrix are also kept, except where
% the fill according to indices overwrites elements. If the user provides 
% an additional single value, that value is the value used for the 
% remaining values in the matrix

% Read out the number of optional arguments provided
N_args = length(varargin);

% Process the additional arguments
matrix_provided = 0;
defaultval_provided = 0;
for k = 1:N_args
    
    if numel(varargin{k}) > 1    % Check if this is a matrix type input
        [Ny, Nx] = size(varargin{k});
        A = varargin{k};
        matrix_provided = 1;
    else                         % or a scalar type input
        defaultval = varargin{k};
        defaultval_provided = 1;
    end
    
end

% Define a default value if one wasn't provided
if ~defaultval_provided
    defaultval = 0;
end

% Initialise a matrix if one wasn't provided as an input
if ~matrix_provided
    Nx = max(inds(:,1));
    Ny = max(inds(:,2));
    A = defaultval * ones(Ny,Nx);
end

% Fill the matrix with the list of values
A(sub2ind([Ny Nx], inds(:,2), inds(:,1))) = vals;



end

