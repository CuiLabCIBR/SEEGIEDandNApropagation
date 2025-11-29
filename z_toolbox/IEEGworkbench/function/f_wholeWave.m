function newEventSampleinfo = f_wholeWave(signal, eventSampleinfo, waveN)
% f_wholeWave - Expands detected event boundaries to include complete waveform cycles
%
% This function takes preliminary IED detections and extends their start/end points
% to include complete waveform oscillations by analyzing local extrema patterns.
%
% Inputs:
%   signal          : Raw SEEG signal (1D vector)
%   eventSampleinfo : Preliminary event boundaries [N x 2 matrix of [start, end] samples]
%   waveN           : Number of adjacent extrema cycles to include (default=3)
%
% Output:
%   newEventSampleinfo : Adjusted event boundaries [N x 2 matrix of [start, end] samples]
%
% Methodology:
%   1. Locates all waveform extrema (peaks/valleys)
%   2. For each preliminary event:
%      a) Identifies the nearest extrema before start/after end
%      b) Expands boundaries to include 'waveN' full cycles
%   3. Merges overlapping extended events
%%% @author: Longzhou Xu
%%% @version: | Last modified: 2025-04-17
%% Initialize parameters and detect extrema
    % Extract initial event boundaries
    EventBegin =  eventSampleinfo(:, 1);
    EventEnd = eventSampleinfo(:, 2);

    % Detect all extrema (peaks and valleys) in the signal
    waveExtremaPointer = f_waveExtremaInflection(signal, 'extrema');
    waveExtremaTime = find(waveExtremaPointer>0);

    % Calculate typical wavelength (4x median inter-extrema distance)
    waveLength = ceil(4*max(diff(waveExtremaTime)));

    % Initialize binary pointer to mark expanded events
    EventCount = length(EventBegin);
    pointer = zeros(size(signal));

%% Process each preliminary event  
    for nEvent = 1:EventCount
        oldBegin = EventBegin(nEvent);
        oldEnd = EventEnd(nEvent);

        % Mark original event region
        pointer(oldBegin:oldEnd)=1;

        % Define search window around the event (±waveLength)
        sampleList = (oldBegin-waveLength):(oldEnd+waveLength);
        sampleList(sampleList<=0) = []; % Remove negative indices
        sampleList(sampleList>length(signal)) = [];  % Remove exceed signal length
        
        % Get extrema within search window
        waveExtremaPointer_soi = waveExtremaPointer(sampleList);
        waveExtremaIndex_soi = find(waveExtremaPointer_soi)+sampleList(1)-1;

        %% Find optimal expansion points
        % Calculate distance from event start to preceding extrema
        CBegin = waveExtremaIndex_soi - oldBegin;
        CBegin(CBegin>0)=min(CBegin); % Only consider extrema BEFORE start
        [~, CUTBegin] = max(CBegin);  % Select closest following extremum

        % Calculate distance from event end to following extrema
        CEnd = waveExtremaIndex_soi - oldEnd;
        CEnd(CEnd<0)=max(CEnd);     % Only consider extrema AFTER end
        [~, CUTEnd] = min(CEnd);    % Select closest following extremum
        
        %% Adjust expansion based on waveform presence
        % Check if original event contains any extrema
        waveExtremaPointer_old = waveExtremaPointer(oldBegin:oldEnd);
        if sum(waveExtremaPointer_old)==0
            CUTBegin = CUTBegin - (waveN+1);
            CUTBegin(CUTBegin<1) = 1;
            CUTEnd = CUTEnd + (waveN+1);
            CUTEnd(CUTEnd>length(waveExtremaIndex_soi)) = length(waveExtremaIndex_soi);
        else
            CUTBegin = CUTBegin - waveN;
            CUTBegin(CUTBegin<1) = 1;
            CUTEnd = CUTEnd + waveN; 
            CUTEnd(CUTEnd>length(waveExtremaIndex_soi)) = length(waveExtremaIndex_soi);
        end

        %% Apply boundary expansion
        % New start includes preceding cycles
        newBegin = waveExtremaIndex_soi(CUTBegin);
        pointer(newBegin:oldBegin) = 1;
        % New end includes following cycles
        newEnd = waveExtremaIndex_soi(CUTEnd);
        pointer(oldEnd:newEnd) = 1;
    end

    %% Generate final event boundaries
    % Merge overlapping regions and get continuous events
    [NewEventBegin, NewEventEnd] = f_eventDetection(pointer, 0, 0);

    % Format output matrix
    newEventSampleinfo(:, 1) = NewEventBegin;
    newEventSampleinfo(:, 2) = NewEventEnd;
end