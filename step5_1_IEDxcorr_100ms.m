% Calculate the maximum cross-correlation and corresponding time delays 
% between electrode channels during IED events by analyzing signal windows 
% around detected spikes and computing correlations at various time lags。
clc; clear; close all;
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
%% ==== Main Processing Loop =====
parfor nS = 1:length(subG) % Parallel loop through subjects
    subID = subG{nS}; % Current subject ID
    T = 100; % Time window size (100ms)
    IEDxcorr(subID, T);
end
%% ===== IED Cross-Correlation Processing Function =====
function IEDxcorr(subID, T)
    % Calculate cross-correlations for IED events in a single subject
    % Inputs:
    %   subID - Subject identifier string
    %   T     - Time window size in milliseconds
    lags = 0 : T;
    %% ===== Data Loading =====
    % Find preprocessed EEG files for current subject
    iEEGfileDir = dir(fullfile('step2_SEEG_TCAR_BPF',  subID, [subID, '_task-rest_run-*_TCAR_BPF.mat']));
    % Load IED detection results for current subject
    load(fullfile('step3_IEDs', [subID, '_task-awake_ref-TCA_IEDs-diffTD.mat']));
    %% ===== Initialize Storage Variables =====
    chanCount = length(IEDs{1}{1}.label);
    IEDxcDelay = cell(chanCount, chanCount);
    IEDxcR = cell(chanCount, chanCount);
    %% ===== Process Each Recording Session =====
    for nRun = 1:length(IEDs)
        disp(['********', subID]);
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
        %% Calculate the maximum cross-correlation and corresponding time delays
        IEDsinfo = IEDs{nRun}{6}.IEDsinfo{1};
        IEDraster = IEDs{nRun}{6}.IEDraster{1};
        signals = data.trial{1};
        signalL = size(signals, 2);
        for nIED = 1:size(IEDsinfo, 1)
                B = IEDsinfo(nIED, 1); E = IEDsinfo(nIED, 2);
                % 1. Find spikes occurring within current IED event
                idx = find(IEDraster(:,1)>B & IEDraster(:,1)<E);
                IEDtime = IEDraster(idx, 1);
                IEDchan = IEDraster(idx, 2);
                for n = 1:length(IEDtime)
                        % 2. Skip spikes near signal boundaries
                        if (IEDtime(n) < 201) || (IEDtime(n) > (signalL-200))
                            continue;
                        end
                        % 3. Extract Signal Window
                        winB = IEDtime(n)-100; winE = IEDtime(n)+100;
                        % Extract signal from primary channel
                        IEDchan1 = IEDchan(n);
                        IEDsignals1 = signals(IEDchan1, winB:winE);
                        % Get unique channels involved in current IED event
                        IEDchan2 = unique(IEDchan);
                        IEDchan2(IEDchan2==IEDchan1) = [];
                        % 4. Precompute Shifted Signals
                        IEDsignals2 = [];
                        for t = 1:length(lags)
                            IEDsignals2(:, :, t) = signals(IEDchan2, winB+lags(t):winE+lags(t));
                        end
                        % 5. Calculate Cross-Correlations
                        for m = 1:length(IEDchan2)
                                IEDchan3 = IEDchan2(m);
                                % Prepare signals for correlation
                                x = IEDsignals1';
                                y = squeeze(IEDsignals2(m, :, :));
                                % Calculate correlation between signals at different lags
                                r = corr(x, y);
                                % Find maximum correlation and corresponding lag
                                [rMax, b]= max(r);
                                % Store results
                                IEDxcDelay{IEDchan1, IEDchan3} = [IEDxcDelay{IEDchan1, IEDchan3}, lags(b)];
                                IEDxcR{IEDchan1, IEDchan3} = [IEDxcR{IEDchan1, IEDchan3}, rMax];
                        end
                end
        end
    end
    % ===== Save Results for Current Subject =====
    saveFolder = 'step5_IEDxcorr';
    mkdir(saveFolder);
    saveName = [subID, '_6d0IED_xcorr', num2str(T), 'ms_Rmax_delay.mat'];
    save(fullfile(saveFolder, saveName), 'IEDxcDelay', 'IEDxcR');
end