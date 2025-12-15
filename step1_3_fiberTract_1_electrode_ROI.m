% This MATLAB code processes neuroimaging and electrode data for a cohort 
% of subjects, handling MRI and DWI data transfer, electrode-to-ROI mapping,
% and ROI image construction for subsequent analysis.
clc; clear; close all; addpath z_toolbox;
addpath 'H:\project_manger\xlz_toolbox\fieldtrip-20231220';% add fieldtrip
ft_defaults;
subG = {'sub-0001', 'sub-0002', 'sub-0003', 'sub-0004', 'sub-0005', ...
                'sub-0008', 'sub-0009', 'sub-0010', 'sub-0011', ...
                'sub-0012', 'sub-0013', 'sub-0014', 'sub-0015', 'sub-0016', ...
                'sub-0017', 'sub-0019', 'sub-0020', 'sub-0021', 'sub-0022', ...
                'sub-0023', 'sub-0028', 'sub-0031', 'sub-0033', ...
                'sub-0035', 'sub-0036', ...
                'sub-0050',  ...
                'sub-0067', 'sub-0074', ...
                'sub-0087', 'sub-0089', 'sub-0092', 'sub-0093', 'sub-0094'};
% ===== Data Path Configuration =====
bsPath = 'E:\brainstrom_db\seeg_awake_xuanwu';
qsiPrepPath = 'E:\seegDataset_xuanwu\step3_DSI\qsiprep_250530_SyNSDC';
eleLocPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
op = 'step2_fiberTrack\DSI';
% ==== Main Processing Loop =====
for nS = 1:length(subG)
    subID = subG{nS}; disp(subID);
    % ---- MRI Data Loading ----
    % Load MRI data from Brainstorm database
    bsMRIdir1 = dir(fullfile(bsPath, 'anat', subID, 'subjectimage_MRI.mat'));
    bsMRIdir2 = dir(fullfile(bsPath, 'anat', subID, 'subjectimage_MRI_T1.mat'));
    if ~isempty(bsMRIdir1)
        smri = load(fullfile(bsMRIdir1.folder, bsMRIdir1.name));
    elseif ~isempty(bsMRIdir2)
        smri = load(fullfile(bsMRIdir2.folder, bsMRIdir2.name));
    end
    % ---- Image Format Conversion ----
    AnatFilename = [subID, '-space-brainstorm_T1.nii'];
    T1SavePath = fullfile(op, subID);
    mkdir(T1SavePath);
    bsT1wFile = fullfile(T1SavePath, AnatFilename);
    out_mri_nii(smri, bsT1wFile, 'uint8');
    % ---- QSIPrep Data Transfer ----
    fileTypes = {'anat', [subID '_ses-001_space-ACPC_desc-preproc_T1w.nii.gz']; % T1 file
                 'dwi',  [subID '_ses-001_run-001_space-ACPC_desc-preproc_dwi.nii.gz']; % qsiprep DWI image
                 'dwi',  [subID '_ses-001_run-001_space-ACPC_desc-preproc_dwi.bvec']; % b vector file
                 'dwi',  [subID '_ses-001_run-001_space-ACPC_desc-preproc_dwi.bval']; % b value file
                 'dwi',  [subID '_ses-001_run-001_space-ACPC_desc-preproc_dwi.b_table.txt']}; % b-table file
    for ft = 1:size(fileTypes, 1)
        src = fullfile(qsiPrepPath, subID, 'ses-001', fileTypes{ft,1}, fileTypes{ft,2});
        dest = fullfile(T1SavePath, fileTypes{ft,2});
        if exist(src, 'file')
            copyfile(src, dest);  % Transfer preprocessed files
        else
            warning('File not found: %s', src);
        end
    end
    % ---- Electrode-to-ROI Mapping ----
    % Load electrode contact anatomy data
    load(fullfile(eleLocPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']), 'contAnat');
    worldXYZ = cell2mat(contAnat.World);
    % Read MRI volumes
    bsT1 = ft_read_mri(bsT1wFile);
    qpT1 = ft_read_mri(fullfile(T1SavePath, [subID '_ses-001_space-ACPC_desc-preproc_T1w.nii.gz']));
    % ROI voxel
    [roiVoxelXYZ, roiCoreWorldXYZ] = f_xyz2roi_withCoreg(bsT1, qpT1, worldXYZ, 2);
    % ---- ROI Image Construction ----
    roiImageAnat = zeros(size(qpT1.anatomy));
    chanLabel = contAnat.label;
    for n = 1:length(chanLabel)
        for v = 1:size(roiVoxelXYZ{n}, 1)
            x = roiVoxelXYZ{n}(v, 1);
            y = roiVoxelXYZ{n}(v, 2);
            z = roiVoxelXYZ{n}(v, 3);
            roiImageAnat(x, y, z) = n;
        end
    end
    qpT1.anatomy = int16(roiImageAnat);
    chanROIimageFile = fullfile(op, subID, [subID, '_space-ACPC_roi.nii']);
    ft_write_mri(chanROIimageFile, qpT1, 'dataformat', 'nifti2');
    % ---- ROI Label Export ----
    ROIlabelFile = fullfile(op, subID, [subID, '_space-ACPC_roi.txt']);
    fileID = fopen(ROIlabelFile, 'w');
    if fileID == -1, error('Wrong txt file'); end
    for n = 1:length(chanLabel)
        fprintf(fileID, '%s\n', [num2str(n), ' ', chanLabel{n}]);
    end
    fclose(fileID);
end