function seriesStruct = f_recursiveSeriesStruct(Series)
% RECURSIVELY PROCESSES HIERARCHICAL DATA STRUCTURES
% This function converts tabular hierarchical data into a nested struct format
% representing parent-child relationships through recursive decomposition.
% Inputs:
%   Series - Cell array
%            Column 1: Parent identifiers
%            Columns 2..end: Child identifiers at subsequent levels
%
% Output:
%   seriesStruct - Struct array with fields:
%       .Leaders: Parent identifier at current level
%       .Childrens: Unique direct children of current parent
%       .ChildrenSeriesStruct: Recursive struct of child hierarchies
%
% Termination Condition: 
%   Stops recursion when input has no children columns (size(Series,2) == 1)
%   or when input data is empty.
    %% TERMINATION CONDITION
    if isempty(Series)
        seriesStruct = [];
        return;
    end
    %% INITIALIZE STRUCTURE
    Leaders = unique(Series(:, 1));
    seriesStruct = struct('Leaders', cell(size(Leaders)), ...
                         'Childrens', [], ...
                         'ChildrenSeriesStruct', []);
    %% PROCESS EACH PARENT NODE
    for nL = 1:length(Leaders)
        % Identify rows belonging to current parent
        Lia = ismember(Series(:, 1), Leaders{nL});

        % Extract child data (remove parent column)
        ChildrenSeries = Series(Lia, 2:end);
        
        %% RECURSIVE PROCESSING
        % Conditionally process child data if valid
        if ~isempty(ChildrenSeries)
            childStruct = f_recursiveSeriesStruct(ChildrenSeries);
            seriesStruct(nL).Leaders = Leaders{nL};
            seriesStruct(nL).Childrens = unique(ChildrenSeries(:, 1));
            seriesStruct(nL).ChildrenSeriesStruct = childStruct;
        else
            childStruct = {};
            seriesStruct(nL).Leaders = Leaders{nL};
            seriesStruct(nL).Childrens = [];
            seriesStruct(nL).ChildrenSeriesStruct = childStruct;
        end
    end
end
