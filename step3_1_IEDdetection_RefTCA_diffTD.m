% Perform automated detection of IEDs from SEEG data across multiple 
% subjects using bandpass filtering (10-60 Hz) and a parallelized 
% peak-detection algorithm with varying thresholds.
clc; clear; close all;  addpath z_toolbox;
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
seegPath = 'step2_SEEG_TCAR_BPF';
% ==== Main Processing Loop =====
for nS = 1:length(subG)
    subID = subG{nS}; 
    % Locate preprocessed EEG files for the subject
    seegDir = dir(fullfile(seegPath, subID, [subID, '_task-rest_run-*_TCAR_BPF.mat']));
    IEDs = cell(length(seegDir), 1);
    for nRun = 1:length(seegDir)
        load(fullfile(seegDir(nRun).folder, seegDir(nRun).name))
        % --- Bandpass Filtering ---
        cfg = [];
        cfg.detrend = 'yes';
        cfg.demean = 'yes';
        cfg.baselinewindow = 'all';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [10, 60]; 
        cfg.bpfiltord = 3; 
        data = ft_preprocessing(cfg, data);
        % --- IED detection ----
        pksDelta = 3; 
        slidingTwin = 30;
        IEDtimeScale = 0.1;
        TD = 5:0.2:9;
        IEDtd = cell(length(TD), 1);
        parfor nTD = 1:length(TD)
            disp([num2str(nTD), '----', subID]);
            p2pDelta = TD(nTD);
            IEDtd{nTD} = f_IED_detection(data, pksDelta, p2pDelta, slidingTwin, IEDtimeScale);
        end
        IEDs{nRun} = IEDtd;
    end
    % Save results
    savefolder = 'step3_IEDs';
    mkdir(savefolder);
    savename = [subID, '_task-awake_ref-TCA_IEDs-diffTD.mat'];
    save(fullfile(savefolder, savename), 'IEDs', 'TD');
end