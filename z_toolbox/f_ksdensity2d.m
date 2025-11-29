function [densityMesh, Xmesh, Ymesh] = f_ksdensity2d(X, Y, xbin, ybin)
%F_KSDENSITY2D Computes 2D kernel density estimate for scatter plot data
%
%   This function calculates a smoothed density distribution of 2D point data
%   using kernel density estimation (KDE), returning both the density values and
%   corresponding grid coordinates.
%
%   Inputs:
%       X, Y   - Vectors of equal length containing x/y coordinates of data points
%       xbin   - Number of bins along x-axis (or vector of bin edges)
%       ybin   - Number of bins along y-axis (or vector of bin edges)
%
%   Outputs:
%       density - Matrix of density values at grid points (size: ybin × xbin)
%       Xmesh   - X-coordinates grid for density values (size: ybin × xbin)
%       Ymesh   - Y-coordinates grid for density values (size: ybin × xbin)
%
%   Algorithm:
%       1. Creates evaluation grid covering data range
%       2. Computes density using Gaussian kernels with optimal bandwidth
%       3. Returns reshaped results for direct plotting
%
%   Example:
%       [density, Xg, Yg] = f_ksdensity2d(randn(1000,1), randn(1000,1), 50, 50);
%       contourf(Xg, Yg, density);
%%
    % Create linearly spaced evaluation grid covering data range
    xi = linspace(min(X), max(X), xbin); % x-axis grid coordinates
    yi = linspace(min(Y), max(Y), ybin); % y-axis grid coordinates

    % Generate full 2D grid using meshgrid
    [Xmesh, Ymesh] = meshgrid(xi, yi); % Creates ybin × xbin matrices

    % Compute 2D kernel density estimate:
    % - Uses optimal bandwidth selection (Silverman's rule by default)
    % - Evaluates at all grid points simultaneously
    % - Returns both density values and used bandwidth (stored but unused here)
    [density, bandwidth] = ksdensity([X, Y], [Xmesh(:), Ymesh(:)]);

    % Reshape density vector into 2D matrix matching grid dimensions
    densityMesh = reshape(density, size(Xmesh)); % Convert to ybin × xbin matrix

end