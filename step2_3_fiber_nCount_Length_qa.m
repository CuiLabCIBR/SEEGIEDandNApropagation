% Extract fiber connectivity matrices (count, length, and 
% quality metrics) between electrode contacts from precomputed DSI 
% track data.
clc; clear; close all;
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
fiberPath = 'step2_fiberTrack';
seegPath = 'step2_SEEG_TCAR_BPF';
elecPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
% Initialize cell arrays for storing fiber metrics 
fiberCount = cell(length(subG), 1);
fiberLength = cell(length(subG), 1);
fiberqa = cell(length(subG), 1);
% ==== Main Processing Loop =====
for nS = 1:length(subG)
    subID = subG{nS};
    % Load preprocessed SEEG data and anatomical info of electrode contacts
    iEEGfileDir = dir(fullfile(seegPath, subID, [subID, '_task-rest_run-*_TCAR_BPF.mat']));
    load(fullfile(iEEGfileDir(1).folder, iEEGfileDir(1).name));
    load(fullfile(elecPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));
    labels = data.label;
    %% ----- Fiber Count Matrix -----
    fiberName = [subID, '_ses-001_run-001_space-ACPC_desc-preproc_dwi.tt.gz.', ...
        subID, '_space-ACPC_roi.ncount.pass.connectivity.mat'];
    if isempty(dir(fullfile(fiberPath, 'DSI', subID, fiberName)))
        fiberCount{nS} = [];
        continue;
    end
    load(fullfile(fiberPath, 'DSI', subID, fiberName), 'connectivity');
    % Build fiber count matrix for electrode pairs
    fiberCount{nS} = zeros(length(labels), length(labels));
    for nChan1 = 1:length(labels)
        for nChan2 = 1:length(labels)
            Lia1 = ismember(contAnat.label, labels{nChan1});
            Lia2 = ismember(contAnat.label, labels{nChan2});
            fiberCount{nS}(nChan1, nChan2) = connectivity(Lia1, Lia2);
        end
    end
    %% ----- Fiber Length Matrix -----
    fiberName = [subID, '_ses-001_run-001_space-ACPC_desc-preproc_dwi.tt.gz.', ...
        subID, '_space-ACPC_roi.mean_length.pass.connectivity.mat'];
    load(fullfile(fiberPath, 'DSI', subID, fiberName), 'connectivity');
    fiberLength{nS} = zeros(length(labels), length(labels));
    for nChan1 = 1:length(labels)
        for nChan2 = 1:length(labels)
            Lia1 = ismember(contAnat.label, labels{nChan1});
            Lia2 = ismember(contAnat.label, labels{nChan2});
            fiberLength{nS}(nChan1, nChan2) = connectivity(Lia1, Lia2);
        end
    end
    %% ----- Fiber Quality Matrix -----
    fiberName = [subID, '_ses-001_run-001_space-ACPC_desc-preproc_dwi.tt.gz.', ...
        subID, '_space-ACPC_roi.qa.pass.connectivity.mat'];
    load(fullfile(fiberPath, 'DSI', subID, fiberName), 'connectivity');
    fiberqa{nS} = zeros(length(labels), length(labels));
    for nChan1 = 1:length(labels)
        for nChan2 = 1:length(labels)
            Lia1 = ismember(contAnat.label, labels{nChan1});
            Lia2 = ismember(contAnat.label, labels{nChan2});
            fiberqa{nS}(nChan1, nChan2) = connectivity(Lia1, Lia2);
        end
    end
end
% ===== Save Results =====
saveName = 'FiberMatrix_ncount_length_qa_allsubject.mat';
save(fullfile(fiberPath, saveName), 'subG',  'fiberCount',  'fiberLength',  'fiberqa');