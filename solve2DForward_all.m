function [pred] = solve2DForward_all(N,mua)
%SOLVE2DFORWARD_ALL solves the 2D forward model of the diffusion eqn using 
%the finite volume method. This solves for all 4 sources where the sources
%are located in the middle of each bounary (west, east, north, south).
%   Inputs:
%       N: number of grid points
%       mua: N^2-by-1 vector of the mua profile
%   Outputs:
%       pred: 64-by-1 vector with 
%           the first 16 measurements correspond to the source on the west boundary
%           the second 16 measurements correspond to the source on the east boundary
%           the third 16 measurements correspond to the source on the north boundary
%           the fouth 16 measurements correspond to the source on the south boundary

global L x y dx dy Reff qin itermax tol Q

% all source locations
source_locs = ['W' 'E' 'N' 'S'];

% the prediction of the detector readings given a mua profile
pred = [];
for i=1:length(source_locs)
    source_loc = source_locs(i);
    [~,u] = solve2DForward(N,mua,source_loc);   % solve forward to get fluence
    P = Q*u;                                    % prediction of 1 source's partial currents
    pred = [pred;P];                            % prediction of all sources
end

end

