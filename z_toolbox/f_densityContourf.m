function f_densityContourf(X, Y, xbin, ybin)
%F_DENSITYCONTOURF Creates a filled contour plot of 2D point density
%
%   This function visualizes the density distribution of 2D data points using
%   filled contours with a custom color gradient, providing a smooth
%   representation of point concentration.
%
%   Inputs:
%       X, Y   - Vectors of equal length containing x/y coordinates of data points
%       xbin   - Number of bins or bin edges along x-axis for density estimation
%       ybin   - Number of bins or bin edges along y-axis for density estimation
%
%   Algorithm:
%       1. Computes 2D kernel density estimate using f_ksdensity2d
%       2. Generates a smooth color gradient from white to dark green
%       3. Creates filled contours without border lines
%       4. Applies the custom colormap to the density plot
%
%   Example:
%       X = randn(1000,1);
%       Y = 0.5*X + randn(1000,1);
%       f_densityContourf(X, Y, 50, 50);
%       colorbar; title('Point Density Distribution');
%% Main Function
    % Compute 2D kernel density estimate
    % Returns:
    %   densityMesh - Matrix of density values
    %   Xmesh, Ymesh - Grid coordinates for density values
    [densityMesh, Xmesh, Ymesh] = f_ksdensity2d(X, Y, xbin, ybin);

    % Define custom color gradient (commented section shows alternative purple-yellow)
    colorList=[0.2700         0    0.3300
                0.2700    0.2300    0.5100
                0.1900    0.4100    0.5600
                0.1200    0.5600    0.5500
                0.2100    0.7200    0.4700
                0.5600    0.8400    0.2700
                0.9900    0.9100    0.1300];

    % Create color interpolation function
    colorFunc=f_colorFuncFactory(colorList);

    % Generate 256-color gradient from the defined color stops
    colorList=colorFunc(linspace(0, 1, 256)');
    colorList(1, :) = [1 1 1];

    % Create filled contour plot with:
    %   - Computed density values
    %   - No contour lines ('LineColor','none')
    %   - Automatic level selection based on data range
    contourf(Xmesh, Ymesh, densityMesh, 'LineColor', 'none');

    % Apply the custom colormap
    colormap(colorList);
end