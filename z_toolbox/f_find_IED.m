function IEDraster = f_find_IED(Data, delta1, delta2, slidingTwin, IEDtimeScale)
%F_FIND_IED Detect Interictal Epileptiform Discharges (IEDs) in multi-channel EEG data.
%
%   This function analyzes SEEG signals to identify IED events using a 
%   Peak Prominence Detection Pipeline:
%   1. Apply Z-score normalization to the input signal;
%   2. Calculate the dynamic boundary using sliding time window: 
%       Upper boundary = rolling_mean + delta1 * rolling_std;
%       Lower boundary = rolling_mean + delta1 * rolling_std;
%   3. Find positive and negative peaks whose value exceed the boundary;
%   4. Filter peaks whose prominence is greater than delat2 * rolling_std;
%   5. Identify negative peaks where:
%       i. At least one positive peak exists within the IEDtimeScale
%       window;
%       ii. The amplitude difference between positive and negative peaks
%       exceeds delta2 × rolling_std;
%   6. Merge adjacent negative peaks within the IEDtimeScale window;
%
%   Inputs:
%       Data            : Fieldtrip style EEG data                    
%       delta1          : Threshold for peak prominence 
%       delta2          : Threshold for positive-to-negative peak amplitude difference 
%       slidingTwin     : Time window size for dynamic thresholding (seconds)
%       IEDtimeScale    : Time window for IED event detection around negative peaks (seconds)
%
%   Output:
%       IEDraster   : Cell array (one per trial). Each cell is N×2 matrix:
%                     - Column 1: IED event time (samples)
%                     - Column 2: Channel index
%
%   Note:
%   It should be noted that proir to this processing, the SEEG data inputs
%   undergo preprocessing through a 10-60Hz bandpass filter. This step 
%   enhance the detection accuracy of IED by isolating the frequency range
%   most characteristic of pathological neural activity.
%   @Ref: R. Janca et al., Detection of interictal epileptiform discharges 
%   using signal envelope distribution modelling: application to epileptic 
%   and non-epileptic intracranial recordings. Brain Topogr 28, 172-183 (2015).
%
%   Edit by xlzhope, 20250812
%   Check by deepseek, 20250812
    %% ==== Parameter Initialization ====
    fsample = Data.fsample;                 % Get sampling frequency
    slidingTwin = ceil(slidingTwin*fsample);    % Convert window size to samples
    IEDtimeScale = ceil(IEDtimeScale*fsample);  % Convert window size to samples
    trialCount = length(Data.trial);        % Number of trials
    chanNum = size(Data.trial{1}, 1);       % Number of channels
    IEDraster = cell(trialCount, 1);        % Preallocate output cell array

    %% ==== Main Detection Loop ====
    for nTrial = 1:trialCount         % Process each trial
        IEDraster{nTrial} = []; 
        for nChan = 1:chanNum         % Process each channel

            % --- Apply Z-score normalization to the input signal ---
            signal = Data.trial{nTrial}(nChan, :);
            signal = zscore(signal);

            % --- Calculate the dynamic boundary ---
            signalMeanCurve = movmean(signal, slidingTwin);
            signalStdCurve = movstd(signal, slidingTwin);
            TDCurve1 = signalMeanCurve + delta1 * signalStdCurve;
            TDCurve2 = signalMeanCurve - delta1 * signalStdCurve;

            % --- Threshold Crossing Detection ---
            sigPointer = zeros(size(signal));
            sigPointer(signal>=TDCurve1) = 1;
            sigPointer(signal<=TDCurve2) = 1;
            
            % --- Peak Detection ---
            [~, locsP, ~, promP] = findpeaks(signal);
            [~, locsN, ~, promN] = findpeaks(-1*signal);
            TDCurve3 = delta1*signalStdCurve;

            % Filter positive peaks by prominence
            promCurveP = zeros(size(signal));
            promCurveP(locsP) = promP;
            promPointerP = (promCurveP >= TDCurve3) & (sigPointer == 1); % Validated peaks

            % Filter negative peaks by prominence
            promCurveN = zeros(size(signal));
            promCurveN(locsN) = promN;
            promPointerN = (promCurveN >= TDCurve3) & (sigPointer == 1); % Validated peaks
            
            % --- IED Candidate Identification ---
            signalL = length(signal);
            pos2negChange = zeros(1, signalL);      % Amplitude difference storage
            locsNcandi = find(promPointerN==1);     % Valid negative peak locations
            
            % --- Calculate Amplitude Difference ---
            for nN = 1:length(locsNcandi)
                % Negative Peak Value
                negPksLocs = locsNcandi(nN);
                negPksValue = signal(negPksLocs);
                % Find nearby Postive Peak
                winB = negPksLocs - IEDtimeScale;
                winB(winB<=0) = 1;
                winE = negPksLocs + IEDtimeScale;
                winE(winE>=signalL) = signalL;
                promPointerP_win = promPointerP(winB:winE);
                signal_win = signal(winB:winE);
                if sum(promPointerP_win) >= 1
                    posPksValue = signal_win(promPointerP_win==1);
                    pos2negChange(negPksLocs) = max(posPksValue-negPksValue);
                end
            end

            % --- Apply Amplitude Difference Threshold ---
            TDCurve4 = delta2*signalStdCurve;        % Amplitude difference threshold
            pos2negChange2 = pos2negChange;
            pos2negChange2(pos2negChange2 < TDCurve4) = 0;  % Thresholding
            pos2negChange2 = mergeNearEvent(pos2negChange2, 2*IEDtimeScale);

            % ==== Output IED raster ====
            IEDtime = find(pos2negChange2 >= TDCurve4);   
            IEDchan = nChan*ones(size(IEDtime));
            IEDraster{nTrial} = [IEDraster{nTrial}; [IEDtime', IEDchan']];

            % % Visualization IED detection pipeline
            % t = Data.time{nTrial};
            % figure; winB = 50000; winE = 110000;
            % subplot(5, 1, 1);
            % plot(t(winB:winE), signal(winB:winE), 'Color', [0, 0, 1]); hold on;
            % plot(t(winB:winE), TDCurve1(winB:winE), 'Color', [1, 0, 0]); 
            % plot(t(winB:winE), TDCurve2(winB:winE), "Color", [1, 0, 0]); hold off;
            % xlim([t(winB), t(winE)]); ylabel('Signal');
            % ylim([min(signal(winB:winE)), max(signal(winB:winE))]);
            % set(gca, 'box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 2);
            % 
            % subplot(5,1,2); 
            % plot(t(winB:winE), promCurveP(winB:winE), 'Color', [0, 0, 1]); hold on;
            % plot(t(winB:winE), TDCurve3(winB:winE), 'Color', [1, 0, 0]); hold off;
            % xlim([t(winB), t(winE)]); ylabel('Pos. Prominence');
            % ylim([min(promCurveP(winB:winE)), max(promCurveP(winB:winE))]);
            % set(gca, 'box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 2);
            % 
            % subplot(5,1,3); 
            % plot(t(winB:winE), promCurveN(winB:winE), 'Color', [0, 0, 1]); hold on;
            % plot(t(winB:winE), TDCurve3(winB:winE), 'Color', [1, 0, 0]); hold off;
            % xlim([t(winB), t(winE)]); ylabel('Neg. Prominence');
            % ylim([min(promCurveN(winB:winE)), max(promCurveN(winB:winE))]);
            % set(gca, 'box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 2);
            % 
            % subplot(5,1,4);
            % plot(t(winB:winE), pos2negChange(winB:winE), 'Color', [0, 0, 1]); hold on;
            % plot(t(winB:winE), TDCurve4(winB:winE), 'Color', [1, 0, 0]); hold off;
            % xlim([t(winB), t(winE)]); ylabel('Amp. Diff.');
            % ylim([min(pos2negChange(winB:winE)), max(pos2negChange(winB:winE))]);
            % set(gca, 'box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 2);
            % 
            % subplot(5,1,5);
            % plot(t(winB:winE), pos2negChange2(winB:winE), 'Color', [0, 0, 1]); hold on;
            % plot(t(winB:winE), TDCurve4(winB:winE), 'Color', [1, 0, 0]); hold off;
            % xlim([t(winB), t(winE)]); xlabel('Time (s)'); ylabel('Amp. Diff.');
            % ylim([min(pos2negChange2(winB:winE)), max(pos2negChange2(winB:winE))]);
            % set(gca, 'box', 'off', 'FontName', 'Arial', 'FontSize', 16, 'LineWidth', 2);
        end
    end
end
%% ==== Helper Function: Merge Nearby Events ====
function Event = mergeNearEvent(Event, winSizeTD)
% MERGENEAREVENT Combine proximate IED events to avoid duplicates.
%   Inputs:
%       Event     : Binary event vector (1 = IED occurrence)
%       winSizeTD : Merge window size (samples)
%   Output:
%       Event     : Merged event vector
    eventIdx = find(Event>0);
    if ~isempty(eventIdx)
        eventIdxDiff = diff([eventIdx(1)-winSizeTD-1, eventIdx]);
        Pointer = zeros(size(eventIdxDiff));
        Pointer(eventIdxDiff<winSizeTD) = 1;
        [B, E, ~] =  f_eventDetection(Pointer, 0, 0);
        for iB = 1:length(B)
            winIdx = [B(iB)-1, B(iB):E(iB)];
            EventWin = Event(eventIdx(winIdx));
            EventWin(EventWin<max(EventWin)) = 0;
            Event(eventIdx(winIdx)) = EventWin;
        end
    end
end