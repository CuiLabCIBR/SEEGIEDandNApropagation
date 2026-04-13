clc; clear; close all;
addpath('H:\project_manger\xlz_toolbox\cifti-matlab');
addpath('H:\project_manger\xlz_toolbox\fieldtrip-20231220\external\gifti');
addpath H:\project_manger\xlz_toolbox\fieldtrip-20231220;
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

% ===== Path Configuration =====
ContAnatPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
electrodePath = 'E:\IEDpropagation\step2_fiberTrack';
gmElectrodeFile_L = fullfile(electrodePath, 'Yeo7_Electrode_wbview_all_L.txt');
gmElectrodeFile_R = fullfile(electrodePath, 'Yeo7_Electrode_wbview_all_R.txt');

%% import surface template gii file
Lsurf_pial_file = fullfile(electrodePath, 'S1200.L.pial_MSMAll.32k_fs_LR.surf.gii');
Rsurf_pial_file = fullfile(electrodePath, 'S1200.R.pial_MSMAll.32k_fs_LR.surf.gii');
Lsurf_pial=gifti(Lsurf_pial_file); 
Rsurf_pial=gifti(Rsurf_pial_file);
surfMNI_pial = [Lsurf_pial.vertices; Rsurf_pial.vertices];
Lsurf_veryinflated_file = fullfile(electrodePath, 'S1200.L.very_inflated_MSMAll.32k_fs_LR.surf.gii');
Rsurf_veryinflated_file = fullfile(electrodePath, 'S1200.R.very_inflated_MSMAll.32k_fs_LR.surf.gii');
Lsurf_veryinflated=gifti(Lsurf_veryinflated_file); 
Rsurf_veryinflated=gifti(Rsurf_veryinflated_file);
surfMNI_veryinflated = [Lsurf_veryinflated.vertices; Rsurf_veryinflated.vertices];

%% import yeo7atlas dlabel file
yeo7dlabel = cifti_read(fullfile(electrodePath, 'Yeo2011_7Networks_N1000.dlabel.nii'));
cdata = yeo7dlabel.cdata;

% ==== Main Processing Loop =====
for nSub = 1:length(subG)
    subID = subG{nSub};
    fprintf('Processing %s (%d/%d)...\n', subID, nSub, numel(subG));

    % ---- Load Electrode Contact Data ----
    load(fullfile(ContAnatPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));

    for nChan = 1:size(contAnat, 1)
        label = contAnat.label{nChan};
        LR = contAnat.LR{nChan};
        Yeo7 = contAnat.Yeo7{nChan};
        MNI = contAnat.MNI{nChan};
        if ~strcmp(Yeo7, 'none')
            % map the channel to surface
            distance = sqrt((MNI(1)-surfMNI_pial(:, 1)).^2+(MNI(2)-surfMNI_pial(:, 2)).^2+(MNI(3)-surfMNI_pial(:, 3)).^2);
            distance(cdata==0) = max(distance);
            [a, b] = min(distance);
            if a > 6
                disp([subID, ' chan ', num2str(nChan), ' distance is ', num2str(a), ' !!!!']);
            end
            coorn = surfMNI_veryinflated(b, :);
            if  strcmp(LR, 'L')
                % edit the left brain gm workbench node file
                RGB = [rand(1)*255 rand(1)*255 rand(1)*255];
                atlas = [subID, Yeo7, 'chan', num2str(nChan)];
                fileID = fopen(gmElectrodeFile_L, 'a+');
                fprintf(fileID,'%s\n', atlas);
                fprintf(fileID,'%1.0f', RGB(1));
                fprintf(fileID,'% 1.0f', RGB(2:3));
                fprintf(fileID,'% 5.2f', coorn(1:2));
                fprintf(fileID,'% 5.2f\n', coorn(3));
                fclose(fileID);
            end
            % edit the right brain gm workbench node file
            if strcmp(LR, 'R')
                RGB = [rand(1)*255 rand(1)*255 rand(1)*255];
                atlas = [subID, Yeo7, 'chan', num2str(nChan)];
                fileID = fopen(gmElectrodeFile_R, 'a+');
                fprintf(fileID,'%s\n', atlas);
                fprintf(fileID,'%1.0f', RGB(1));
                fprintf(fileID,'% 1.0f', RGB(2:3));
                fprintf(fileID,'% 5.2f', coorn(1:2));
                fprintf(fileID,'% 5.2f\n', coorn(3));
                fclose(fileID);
            end
        end
    end
end
%% create foci file
wb_command = 'H:\project_manger\xlz_toolbox\workbench\win64\wb_command.exe';
fociL_file = fullfile(electrodePath, 'Yeo7_Electrode_wbview_all.L.32k_fs_LR.foci');
fociR_file = fullfile(electrodePath, 'Yeo7_Electrode_wbview_all.R.32k_fs_LR.foci');
delete(fociR_file); delete(fociL_file);
cmd = [wb_command ' -foci-create -class foci ' gmElectrodeFile_L ' ' Lsurf_veryinflated_file ' ' fociL_file];
system(cmd)
cmd = [wb_command ' -foci-create -class foci ' gmElectrodeFile_R ' ' Rsurf_veryinflated_file ' ' fociR_file];
system(cmd)