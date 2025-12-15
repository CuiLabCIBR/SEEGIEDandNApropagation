% Perform tissue-specific common average rereferencing on SEEG data 
% by separately processing white matter, gray matter, and subcortical 
% electrode contacts before reintegrating them for each subject.
clc; clear; close all; addpath z_toolbox;
addpath H:\project_manger\xlz_toolbox\fieldtrip-20231220;% add fieldtrip
ft_defaults;
addpath z_toolbox\IEEGworkbench;
iEEGworkbench_initial;
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
iEEGPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
for nS = 1:length(subG)  % Iterate through subjects
    subID = subG{nS};
    % Locate SEEG file  and load contact info
    iEEGfile = dir(fullfile(iEEGPath, subID, [subID, '_task-rest_run-*_ieeg_miniPrep_BPfiltering.mat']));
    load(fullfile(iEEGPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));
    parfor nRun = 1:length(iEEGfile)
        SEEG_ICAR_BPFiltering(subID, iEEGfile, nRun, contAnat)
    end
end
%% subfunction for parfor
function SEEG_ICAR_BPFiltering(subID, iEEGfile, nRun, contAnat)
        load(fullfile(iEEGfile(nRun).folder, iEEGfile(nRun).name)); % Load EEG data (FieldTrip format)
        %% ===== White Matter Rereferencing =====
        WMcontact = {};
        a = strcmp('White', contAnat.CWS);
        if any(a)
            WMcontact = contAnat.label(a == 1);
            cfg = [];
            cfg.channel = WMcontact;
            wmdata = ft_preprocessing(cfg, data);
            refMethod = 'CommonAverage';
            wmdata = f_reref_SEEG(wmdata, refMethod);
        end    
        %% ===== Gray Matter Rereferencing =====
        GMcontact = {};
        a = strcmp('Cortex', contAnat.CWS);
        if any(a)
            GMcontact = contAnat.label(a == 1);
            cfg = [];
            cfg.channel = GMcontact;
            gmdata = ft_preprocessing(cfg, data);
            refMethod = 'CommonAverage';
            gmdata = f_reref_SEEG(gmdata, refMethod);
        end
        %% ===== Subcortical Rereferencing =====
        SCcontact = {};
        a = strcmp('Subcortex', contAnat.CWS);
        if any(a)
            refMethod = 'CommonAverage';
            scdata = f_reref_SEEG(data, refMethod);% CAR on full dataset
            SCcontact = contAnat.label(a == 1);
            cfg = [];
            cfg.channel = SCcontact;
            scdata = ft_preprocessing(cfg, scdata);% Extract subcortical subset
        end
        %% ===== Data Reintegration =====
        % Recombine tissue-specific datasets while preserving spatial grouping
        if ~isempty(SCcontact) && ~isempty(GMcontact) && ~isempty(WMcontact)
            data.label = [gmdata.label; wmdata.label; scdata.label];
            data.trial{1} = [gmdata.trial{1}; wmdata.trial{1}; scdata.trial{1}];
        elseif isempty(SCcontact) && ~isempty(GMcontact) && ~isempty(WMcontact)
            data.label = [gmdata.label; wmdata.label];
            data.trial{1} = [gmdata.trial{1}; wmdata.trial{1}];
        elseif ~isempty(SCcontact) && ~isempty(GMcontact) && isempty(WMcontact)
            data.label = [gmdata.label; scdata.label];
            data.trial{1} = [gmdata.trial{1}; scdata.trial{1}];
        elseif ~isempty(SCcontact) && isempty(GMcontact) && ~isempty(WMcontact)
            data.label = [wmdata.label; scdata.label];
            data.trial{1} = [wmdata.trial{1}; scdata.trial{1}];
        else
            error('There were only one-tissue contacts!');
        end
        %% ===== Bandpass Filtering =====
        cfg = [];
        cfg.detrend = 'yes';
        cfg.demean = 'yes';
        cfg.baselinewindow = 'all';
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [3, 160]; 
        cfg.bpfiltord = 3;
        data = ft_preprocessing(cfg, data);
        %% ===== Save Processed Data =====
        saveName = f_strsplit(iEEGfile(nRun).name, '_');
        saveName = [saveName{1}, '_', saveName{2}, '_', saveName{3},  '_TCAR_BPF.mat'];
        saveFolder = fullfile('step2_SEEG_TCAR_BPF', subID);
        mkdir(saveFolder);
        save(fullfile(saveFolder, saveName), 'data');
end