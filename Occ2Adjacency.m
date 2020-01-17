function [A, mapping] = Occ2Adjacency(occ_map, moore)
% This function takes a map of occupied and unoccupied sites, and
% translates it to a graph structure.
%
%   INPUTS
% ------------
% occ_map:    A matrix of zeros and ones, with ones indicating occupied
%             sites
% (moore):    Specifies whether to use a moore neighbourhood (use "1"), or
%             a von Neumann neighbourhood (use "0"). Optional (default: 0).
%
%   OUTPUTS
% ------------
% A:          Adjacency matrix, including 1's on the diagonal
% mapping:    A matrix (same dimensions as input occ_map) 

% Read out the size of the input matrix
[Ny,Nx] = size(occ_map);

% Pad the occupancy map with occupied material (allows for a neater
% implementation)
occ_ext = [ones(1,Nx+2); [ones(Ny,1), occ_map, ones(Ny,1)]; ones(1,Nx+2)];

% Generate a list of unoccupied sites in the extended matrix
unocc_ext = find(~occ_ext);

% Assume a von Neumann neighbourhood if the "moore" argument was not
% supplied
if nargin < 2
    moore = 0;
end

% Initialise neighbours list
if moore
    neighbours = zeros( (Ny+2)*(Nx+2) , 8);
else
    neighbours = zeros( (Ny+2)*(Nx+2) , 4);
end

% For each site, add a neighbour for each cardinal direction (including
% Moore neighbourhood if requested)
neighbours(unocc_ext,1) = (unocc_ext + 1) .* ~occ_ext(unocc_ext + 1);
neighbours(unocc_ext,2) = (unocc_ext - 1) .* ~occ_ext(unocc_ext - 1);
neighbours(unocc_ext,3) = (unocc_ext + Ny + 2) .* ~occ_ext(unocc_ext + Ny + 2);
neighbours(unocc_ext,4) = (unocc_ext - Ny - 2) .* ~occ_ext(unocc_ext - Ny - 2);
if moore
    neighbours(unocc_ext,5) = (unocc_ext + 1 + Ny + 2) .* ~occ_ext(unocc_ext + 1 + Ny + 2);
    neighbours(unocc_ext,6) = (unocc_ext - 1 + Ny + 2) .* ~occ_ext(unocc_ext - 1 + Ny + 2);
    neighbours(unocc_ext,7) = (unocc_ext + 1 - Ny - 2) .* ~occ_ext(unocc_ext + 1 - Ny - 2);
    neighbours(unocc_ext,8) = (unocc_ext - 1 - Ny - 2) .* ~occ_ext(unocc_ext - 1 - Ny - 2);
end

% Loop over all sites, and fill in their neighbours
iv = []; jv = [];
for k = 1:size(neighbours,1)
   
    % Read out unoccupied neighbours here
    neigh_here = neighbours(k, neighbours(k,:)~=0 );
    
    % Build up the list of elements to fill in the sparse adjacency matrix
    iv = [iv; neigh_here']; jv = [jv;k*ones(length(neigh_here),1)];
    
end

% Fill in the sparse matrix
A = sparse(jv,iv,ones(length(jv),1),(Ny+2)*(Nx+2),(Ny+2)*(Nx+2));

% Cut out all the occupied sites
A = A(unocc_ext, unocc_ext);

% Record the mapping used here
mapping = zeros((Nx+2)*(Ny+2),1);             % Initialise
mapping(unocc_ext) = 1:length(unocc_ext);     % Fill out with numbering
mapping = reshape(mapping,Ny+2,Nx+2);         % Convert to matrix
mapping = mapping(2:Ny+1, 2:Nx+1);            % Cut out the "padding" used (all of this is bsaed off of occ_ext, not occ_map)

end