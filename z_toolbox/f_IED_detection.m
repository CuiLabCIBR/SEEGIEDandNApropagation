function IEDs = f_IED_detection(Data, delta1, delta2, slidingTwin, IEDtimeScale)
% F_IED_DETECTION Detect Interictal Epileptiform Discharges (IEDs) in iEEG data
%
%   Inputs:
%       Data            : Fieldtrip style EEG data
%       delta1          : Threshold for peak prominence
%       delta2          : Threshold for peak-to-peak
%       slidingTwin     : Time window size for dynamic thresholding (seconds)
%       IEDtimeScale    : Time window for IED event detection around negative peaks (seconds)
%
%   Output Parameters
%       IEDs          : Structured detection results
%           .time           : Time vector reference
%           .fsample        : Sampling frequency
%           .iEEGsinfo      : Trial metadata
%           .label          : Channel labels
%           .cfg            : Configuration parameters
%           .IEDraster      : Cell array of detected IED events
%           .IEDsinfo       : Event time intervals per trial
%           .IEDchanPartici : Channel participation matrix
%
%   Edit by xlzhope, 20250812
%   Check by deepseek, 20250812
%% ---- Initialization & Status ----
    disp("===============================================================");
    disp('****IED Detection Start!!!!');
    disp('.-.. / --- / -. / --. / --.. / .... / --- / ..- // -..- / ..-');
    disp(['****', char(datetime('now'))]);
    disp("---------------------------------------------------------------");

%% ---- Primary Detection Pipeline ----
    % Initialize output structure with input metadata
    IEDs.time  = Data.time;
    fsample = Data.fsample;
    IEDs.fsample = Data.fsample;
    IEDs.iEEGsinfo = Data.sampleinfo;
    IEDs.label = Data.label;

    % Configure detection parameters
    IEDs.cfg.peakProminenceThreshold = delta1;
    IEDs.cfg.Pos2NegPeakChangeThreshold = delta2;
    IEDs.cfg.slidingTwin = slidingTwin;
    IEDs.cfg.IEDtimeScale = IEDtimeScale;

    % Execute core detection algorithm
    IEDraster  = f_find_IED(Data, delta1, delta2, slidingTwin, IEDtimeScale);
    IEDs.IEDraster = IEDraster;

%% ---- Temporal Event Processing ----
    trialCount = length(IEDraster);
    chanCount = length(IEDs.label);
    for nTrial = 1:trialCount
        sinfo = IEDs.iEEGsinfo(nTrial, :);
        signalLength = sinfo(2) - sinfo(1) + 1;

        % Initialize temporal activity pointer
        pointer = zeros(chanCount, signalLength);
        currIEDraster = IEDraster{nTrial};
        
        % Expand detection windows for temporal context
        for nIED = 1:size(currIEDraster, 1)
            IEDchanIdx = currIEDraster(nIED, 2);
            IEDtimeIdx = currIEDraster(nIED, 1);
            startIdx = max(1, IEDtimeIdx - ceil(IEDtimeScale*fsample));
            endIdx = min(signalLength, IEDtimeIdx + ceil(IEDtimeScale*fsample));
            pointer(IEDchanIdx, startIdx:endIdx) = 1;
        end

        % Detect event boundaries
        [eventStart, eventEnd, ~] = f_eventDetection(sum(pointer, 1)>0, 0, 0);
        IEDs.IEDsinfo{nTrial} = [eventStart', eventEnd'];

        % Create channel participation matrix
        channelParticipation = false(chanCount, 1);
        for nEvent = 1:length(eventStart)
            activeChannels = sum(pointer(:, eventStart(nEvent):eventEnd(nEvent)), 2) > 1;
            channelParticipation(:, nEvent) = activeChannels;
        end
        IEDs.IEDchanPartici{nTrial} = channelParticipation;
    end
%% ---- Completion Status ----
    disp("---------------------------------------------------------------");
    disp('****IED Detection Completed!!!!');
    disp('.-.. / --- / -. / --. / --.. / .... / --- / ..- // -..- / ..-');
    disp(['****', char(datetime('now'))]);
    disp("===============================================================");
end
