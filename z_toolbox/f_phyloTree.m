function f_phyloTree(seriesStruct, AllElements, ElementsColor, RadRange, depth)
% F_PHYLOTREE Recursively draws a phylogenetic tree with circular arcs
% 
% This function visualizes a phylogenetic tree structure by drawing partial
% circles (arcs) for each element. The tree is built recursively based on 
% the provided series structure. Each arc represents a leader element, and 
% child structures are drawn recursively within the parent's arc segment.
%
% Inputs:
%   seriesStruct   : Structure array containing tree data with fields:
%                    - Leaders: The leader element for this series
%                    - ChildrenSeriesStruct: Child structures (recursive)
%   AllElements    : Cell array of all possible element names
%   ElementsColor  : RGB matrix defining colors for each element in AllElements
%   RadRange       : 2-element vector [start_angle, end_angle] in radians
%                    defining the angular range for this recursion level
%   depth          : Current recursion depth (controls circle radius)
%
% The function produces a circular phylogenetic tree visualization without
% returning any values.

    % Handle empty input case
    if isempty(seriesStruct)
        return;
    end
    depth = depth + 1;

    % Extract all Leaders from the series structure
    Leaders = cell(length(seriesStruct), 1);
    for nS = 1:length(seriesStruct)
        Leaders{nS} = seriesStruct(nS).Leaders;
    end

    % Calculate radius based on current depth
    r = 0.01 + (depth - 1) * 0.2;

    % Create angular ranges for each leader element
    Rad = linspace(RadRange(1), RadRange(2), length(Leaders)+1);

    % Process each element in AllElements
    N = 0;
    for nE = 1:length(AllElements)
        % Check if current element is a leader in any series
        Lia = ismember(Leaders, AllElements{nE});
        if sum(Lia) == 1
            N = N + 1;
            RadCurrent = Rad(N); RadTarget = Rad(N+1);

            % Get color for current element
            colorRGB = ElementsColor(nE, :);

            % Draw partial circle (arc) for current element
            partialCircle(r, [RadCurrent RadTarget], colorRGB);

            % Add text label if arc is wide enough
            if RadTarget - RadCurrent > pi / 25
                rot = rad2deg(mean([RadCurrent RadTarget]));
                currentPos = sum([RadCurrent, RadTarget * 2]) / 3;
                x = cos(currentPos) * r;
                y = sin(currentPos) * r;
                % Create text object with formatting
                h = text(x, y, AllElements{nE}, 'FontSize', 16, ...
                    'FontName', 'Arial', 'FontWeight', 'bold');
                h.Rotation = rot;
            end

            % Recursive call for child structures
            childStruct = seriesStruct(Lia).ChildrenSeriesStruct;
            if ~isempty(childStruct)
                f_phyloTree(childStruct, AllElements, ...
                    ElementsColor, [RadCurrent, RadTarget], depth);
            end
        end
    end
    % Hide axes for cleaner visualization
    set(gca, 'Visible', 'off');
end
%% subfunction: partialCircle
function partialCircle(r, th, color)
% PARTIALCIRCLE Draws a filled partial circle (arc) segment
%
% Inputs:
%   r     : Radius of the inner circle
%   th    : 2-element vector [start_angle, end_angle] in radians
%   color : RGB color vector for filling the arc
    x = 0; 
    y = 0; 
    rFactor = 0.19;
    % Calculate inner and outer circle coordinates
    [x1, y1] = circle(x, y, r, th);
    [x2, y2] = circle(x, y, r + rFactor, th);
    % Create polygon coordinates for filling
    x = [fliplr(x1), x2];
    y = [fliplr(y1), y2]; 
    fill(x, y, color); 
end
%% subfunction: circle
function [xunit, yunit] = circle(x, y, r, th, plotting)
% CIRCLE Calculate coordinates for a circle or arc
%
% Inputs:
%   x, y    : Center coordinates of the circle
%   r       : Radius of the circle
%   th      : 2-element vector [start_angle, end_angle] (default: full circle)
%   plotting: Boolean flag to directly plot the circle (default: false)
%
% Outputs:
%   xunit, yunit: Coordinates of the circle/arc points
    hold on % Maintain current plot
    % Set default values for optional parameters
    if nargin <= 4
        plotting = false; 
    end
    if nargin <= 3
        th = [0 2 * pi];
    end
    % Create angular values for the circle/arc
    th = linspace(th(1), th(2), 100); 
    % Calculate circle coordinates
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    % Return if not plotting
    if ~plotting
        return;
    end
    % Plot the circle if requested
    plot(xunit, yunit, 'k'); 
end
