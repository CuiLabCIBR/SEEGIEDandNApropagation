% Performs parallel processing of cross-correlation analyses during non-IED
% periods in EEG data across multiple subjects, identifying IED-free segments,
% applying bandpass filtering, and computing channel-wise maximum correlation
% coefficients and delays with sliding window analysis.
clc; clear; close all; addpath z_toolbox;
addpath 'H:\project_manger\xlz_toolbox\fieldtrip-20231220';
ft_defaults; 
subG = {'sub-0001', 'sub-0002', 'sub-0003', 'sub-0004', 'sub-0005', ...
                'sub-0007', 'sub-0008', 'sub-0009', 'sub-0010', 'sub-0011', ...
                'sub-0012', 'sub-0013', 'sub-0014', 'sub-0015', 'sub-0016', ...
                'sub-0017', 'sub-0019', 'sub-0020', 'sub-0021', 'sub-0022', ...
                'sub-0023', 'sub-0028', 'sub-0029', 'sub-0031', 'sub-0033', ...
                'sub-0034', 'sub-0035', 'sub-0036', 'sub-0039', 'sub-0046', ...
                'sub-0047', 'sub-0050', 'sub-0052', 'sub-0054', 'sub-0062', ...
                'sub-0063', 'sub-0066', 'sub-0067', 'sub-0074', 'sub-0077', ...
                'sub-0087', 'sub-0089', 'sub-0092', 'sub-0093', 'sub-0094', ...
                'sub-0113', 'sub-0115'};
parfor nS = 1:length(subG)
    subID = subG{nS}; % Current subject ID
    T = 100; 
    NAxcorr(subID, T);
end
%% ===== Non-IED Cross-Correlation Function =====
function NAxcorr(subID, T)
    lags = 0:T;
    %% ===== Data Loading =====
    % Find preprocessed EEG files for current subject
    iEEGfileDir = dir(fullfile('step2_SEEG_TCAR_BPF', subID, [subID, '_task-rest_run-*_TCAR_BPF.mat']));
    % Load IED detection results for current subject
    load(fullfile('step3_IEDs', [subID, '_task-awake_ref-TCA_IEDs-diffTD.mat']));
    %% ==== Preprocessing for each run ====
    % Initialize Storage Variables
    NAxcDelay_allRun = []; NAxcR_allRun = []; NA2 = 0;
    for nRun = 1:length(IEDs)
        disp(['****** ', subID]);
        %% Apply bandpass filter (10-60 Hz) to SEEG data
        load(fullfile(iEEGfileDir(nRun).folder, iEEGfileDir(nRun).name));
        cfg = [];
        cfg.detrend = 'yes';
        cfg.demean = 'yes';
        cfg.baselinewindow = 'all';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [10, 60];
        cfg.bpfiltord = 3;
        data = ft_preprocessing(cfg, data);
        signals = data.trial{1};
        %% Identify Non-IED Periods
        % Retrieve IED timing information
        IEDsinfo = IEDs{nRun}{6}.IEDsinfo{1};
        % Create normal activity pointer to identify non-IED periods (1 = non-IED, 0 = IED)
        NApointer = ones(1, size(signals, 2));
        for nIED = 1:size(IEDsinfo, 1)
            NApointer(IEDsinfo(nIED, 1) : IEDsinfo(nIED, 2)) = 0;
        end
        % Detect continuous non-IED segments
        [Event_begin, Event_end, ~] = f_eventDetection(NApointer, 0, 0);
        EventL = Event_end-Event_begin;
        NAidx = find(EventL > 500); % Select segments longer than 500 samples (0.5s in 1000Hz)
        Event_begin = Event_begin(NAidx);
        Event_end = Event_end(NAidx);
        %% ===== Calculate Cross-correlation of Each Non-IED Segment =====
        NAxcDelay = zeros(size(signals, 1), size(signals, 1), length(Event_begin), 'single');
        NAxcR = zeros(size(signals, 1), size(signals, 1), length(Event_begin), 'single');
        for nNA = 1:length(Event_begin)
            % Define analysis window centered at the middle of the non-IED segment
            M = ceil((Event_begin(nNA) + Event_end(nNA)) ./ 2);
            winB = M - 100; winE = M + 100;
            % Preallocate matrix for shifted signals during cross-correlation
            NAsignals = zeros(size(signals, 1), (winE-winB+1), length(lags));
            for t = 1:length(lags)
                NAsignals(:, :, t) = signals(:, winB+lags(t):winE+lags(t));
            end
            % Calculate cross-correlation for each channel pair
            for nChan1 = 1:size(signals, 1)
                NAsignals1 = signals(nChan1, winB:winE)';
                for nChan2 = 1:size(signals, 1)
                    NAsignals2 = squeeze(NAsignals(nChan2, :, :));
                    R = corr(NAsignals1, NAsignals2);
                    [maxR, b] = max(R);
                    NAxcDelay(nChan1, nChan2, nNA) = lags(b);
                    NAxcR(nChan1, nChan2, nNA) = maxR;
                end
            end
        end
        % Accumulate results across runs
        NA2 = NA2 + length(Event_begin);
        NA1 = NA2 - length(Event_begin) + 1;
        NAxcDelay_allRun(:, :, NA1:NA2) = NAxcDelay;
        NAxcR_allRun(:, :, NA1:NA2) = NAxcR;
    end
    % ===== Save Results for Current Subject =====
    saveFolder = 'step6_NAxcorr'; mkdir(saveFolder);
    saveName = [subID, '_NAnon6d0IED_xcorr', num2str(T), 'ms_Rmax_delay.mat'];
    save(fullfile(saveFolder, saveName), 'NAxcDelay_allRun', 'NAxcR_allRun');
end