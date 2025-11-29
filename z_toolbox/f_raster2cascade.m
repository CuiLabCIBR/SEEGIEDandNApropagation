function Cascade = f_raster2cascade(raster, intervalL)
% Converts a raster plot (spike timings) into "cascades" of closely spaced events.
% A "cascade" is defined as a group of consecutive spikes where the interval 
% between adjacent spikes is <= intervalL.
%
% Inputs:
%   raster    : Nx2 matrix, where column 1 = spike times, column 2 = unit IDs.
%   intervalL : Maximum allowed time interval between spikes in a cascade.
%
% Output:
%   Cascade   : Struct array with fields:
%               - time: Array of spike times in the cascade.
%               - unit: Array of corresponding unit IDs.
%%% @author: Longzhou Xu
%%% @version: | Last modified: 2025-04-17
%%
    % Sort raster by spike times (ascending) and store sorted indices
    [rasterTime, idx] = sort(raster(:, 1), 'ascend');

    % Compute time differences between consecutive spikes
    rasterTimeDiff = diff(rasterTime);

    % Identify spikes that are part of a cascade (interval <= intervalL)
    continuePointer = rasterTimeDiff<=intervalL;

    % Mark cascade boundaries: 
    % A = 1  where a cascade starts, A = -1 where it ends
    A = diff([0; continuePointer; 0]);

    % Extract start/end indices of cascades
    CascadeStart = find(A==1);
    CascadeEnd = find(A==-1);

    % Build output struct for each cascade
    Cascade = []; N = 0;
    for nC = 1:length(CascadeStart)
        CS = CascadeStart(nC);
        CE = CascadeEnd(nC);
        unit = raster(idx(CS:CE), 2); % Unit IDs
        if length(unit)>=3
            N = N + 1;
            % Extract spike times and unit IDs for this cascade
            % (Using sorted indices I to preserve original unit IDs)
            Cascade(N).time = raster(idx(CS:CE), 1); % Spike times
            Cascade(N).unit = unit; % Unit IDs
        end
    end
end