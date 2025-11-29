function waveEventPointer = f_waveExtremaInflection(signal, type)
% F_WAVEEXTREMAINFLECTION Detects extrema and inflection points in a signal
%
%   Inputs:
%       signal - Input signal vector
%       type   - Detection type: 'extrema', 'inflection', or 'all'
%
%   Output:
%       wavePointer - Binary vector marking detected points (1 = detected, 0 = not detected)
%
% @author Longzhou Xu
% @version | Last modified: 2025-04-26
%%
    if strcmp(type, 'extrema')
        % EXTRMA DETECTION (local maxima/minima)
        
        % Compute first derivative (difference between consecutive points)
        signalDiff = diff(signal);
        signalDiff = f_replaceZero(signalDiff);
        signalDiff = [signalDiff, signalDiff(end)];      
        
        % Create direction indicator: 1=increasing, -1=decreasing
        Pointer = zeros(size(signalDiff));
        Pointer(signalDiff > 0) = 1;        % Mark increasing segments
        Pointer(signalDiff < 0) = -1;       % Mark decreasing segments
        
        % Find where direction changes (extrema points)
        waveEventPointer = zeros(size(signal));
        waveExtremaIdx = find(diff(Pointer) < 0)+1;
        waveEventPointer(waveExtremaIdx) = 1;    % Mark extrema locations
        waveExtremaIdx = find(diff(Pointer) > 0)+1;
        waveEventPointer(waveExtremaIdx) = -1;    % Mark extrema locations

    elseif strcmp(type, 'inflection')
        % INFLECTION POINT DETECTION (curvature changes)
        
        % Compute second derivative (difference of differences)
        signalDiff2 = diff(diff(signal));
        signalDiff2 = [signalDiff2(1), signalDiff2, signalDiff2(end)];
        signalDiff2 = f_replaceZero(signalDiff2);
        
        % Create curvature direction indicator: 1=convex, -1=concave
        Pointer = zeros(size(signalDiff2));
        Pointer(signalDiff2 > 0) = 1;       % Mark convex segments
        Pointer(signalDiff2 < 0) = -1;      % Mark concave segments
        
        % Find where curvature direction changes (inflection points)
        waveInflectionIdx = find(abs(diff(Pointer))>0);
        waveEventPointer  =zeros(size(signal));
        waveEventPointer(waveInflectionIdx) = 1;

    elseif strcmp(type, 'all')
        % DETECT BOTH EXTREMA AND INFLECTION POINTS
        
        % First derivative processing (for extrema)
        signalDiff = diff(signal);
        signalDiff = f_replaceZero(signalDiff);
        signalDiff = [signalDiff, signalDiff(end)];
        Pointer = zeros(size(signalDiff));
        Pointer(signalDiff > 0) = 1;
        Pointer(signalDiff < 0) = -1;
        waveExtremaIdx = find(abs(diff(Pointer)) > 0)+1;
        
        % Second derivative processing (for inflection points)
        signalDiff2 = diff(diff(signal));
        signalDiff2 = f_replaceZero(signalDiff2);
        signalDiff2 = [signalDiff2(1), signalDiff2, signalDiff2(end)];
        Pointer = zeros(size(signalDiff2));
        Pointer(signalDiff2 > 0) = 1;
        Pointer(signalDiff2 < 0) = -1;
        waveInflectionIdx = find(abs(diff(Pointer))>0)+1;
        
        % Combine both detections
        waveEventPointer = zeros(size(signal));
        waveEventPointer(waveExtremaIdx) = 1;
        waveEventPointer(waveInflectionIdx) = 1;
    end
end

function series = f_replaceZero(series)
% F_REPLACEZERO Replaces zero values in a series with neighboring averages
%
%   Input:
%       series - Input data series with possible zero values
%
%   Output:
%       series - Processed series with zeros replaced
%%
    % Count remaining zeros (loop until all zeros are replaced)
    L = length(find(series==0));
    while L > 0
        A = find(series==0);
        for n = 1:length(A)
            idx = A(n);
            if idx==1
                series(idx) = series(idx+1);
            elseif idx == length(series)
                series(idx) = series(idx-1);
            else
                series(idx) = (series(idx-1)+series(idx+1))/2;
            end
        end
        L = length(find(series==0));
    end
end