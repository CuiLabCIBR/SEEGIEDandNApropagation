function colorFunc=f_colorFuncFactory(colorList)
%F_COLORFUNCFACTORY Creates a color interpolation function from a list of RGB colors.
%
%   This function generates a callable function that maps input values in [0,1]
%   to interpolated RGB colors based on the provided color list.
%
%   Input:
%       colorList - Nx3 matrix where each row is an RGB color (values in [0,1]).
%                  First row corresponds to X=0, last row to X=1.
%
%   Output:
%       colorFunc - Function handle that takes values X in [0,1] and returns
%                   interpolated RGB colors using piecewise cubic Hermite interpolation.
%
%   Example:
%       colors = [1 0 0;  % Red at X=0
%                 0 1 0]; % Green at X=1
%       getColor = f_colorFuncFactory(colors);
%       pink = getColor(0.5); % Returns [0.5 0.5 0] (yellow)
%
%   See also: interp1, pchip
    % Normalized positions for each color in the list (0 to 1)
    x = (0:size(colorList,1)-1)./(size(colorList,1)-1);

    % Separate RGB components into individual vectors
    y1 = colorList(:,1);
    y2 = colorList(:,2);
    y3 = colorList(:,3);

    % Create interpolation function
    % Uses piecewise cubic Hermite interpolation (pchip) for smooth transitions
    % that preserve monotonicity and avoid overshooting
    colorFunc=@(X)[interp1(x, y1, X, 'pchip'), interp1(x, y2, X, 'pchip'), interp1(x, y3, X, 'pchip')];
end