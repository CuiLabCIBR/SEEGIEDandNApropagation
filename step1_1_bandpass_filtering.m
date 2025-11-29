% Performing bandpass filtering (1-170 Hz) on SEEG data, excluding 
% pre-identified bad channels using Fieldtrip toolbox.
clc; clear; close all; addpath z_toolbox;
addpath(fullfile('z_toolbox', 'IEEGworkbench'));
iEEGworkbench_initial;
addpath H:\project_manger\xlz_toolbox\fieldtrip-20231220;% add fieldtrip
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
% ===== Data Path Configuration =====
iEEGPath = 'E:\seegDataset_xuanwu\step7_ageOver16_awakeSEEG_300s_miniPrep';
% ===== Main Processing Loop =====
for ns = 1:length(subG) % Iterate through subjects
    subID = subG{ns};   % Current subject ID
    disp(subID);
    % Find subject's SEEG files
    iEEGfile = dir(fullfile(iEEGPath, subID, [subID, '_task-rest_run-*_ieeg_miniPrep.mat']));
    for nRun = 1:length(iEEGfile) % Process each recording session
        % Load EEG data and bad channel markers
        load(fullfile(iEEGfile(nRun).folder, iEEGfile(nRun).name));
        load(fullfile(iEEGfile(nRun).folder, [subID, '_badchannels.mat']));
        % ===== Preprocessing =====
        % Select channels
        cfg = [];
        cfg.channel = [{'all'}, badchannels];
        data = ft_preprocessing(cfg, data);
        % Bandpass Filtering
        data = f_filter_bandpass(data, [1, 170]);
        % ===== Save Processed Data =====
        saveName = f_strsplit(iEEGfile(nRun).name, '.');
        saveName = [saveName{1}, '_BPfiltering.mat'];
        saveFolder = fullfile('step_1st_BPfiltering', subID);
        mkdir(saveFolder);
        save(fullfile(saveFolder, saveName), 'data');
    end
end