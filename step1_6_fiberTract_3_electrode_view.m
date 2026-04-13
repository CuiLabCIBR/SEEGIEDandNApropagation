% Converting SEEG electrode MNI coordinates to voxel space and generating 
% a CSV file formatted for DSI Studio visualization, including electrode
% geometry and orientation information.
clc; clear; close all;
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

colorS = {'4291267352', '4286858482', '4282708561', '4294057291', '4281893063', ...
    '4293708349', '4282780454', '4291525487', '4280272425', '4293651044', ...
    '4288245501', '4280448417', '4292467262', '4282526581', '4294700112', ...
    '4281262725', '4285520709', '4294127706', '4288240579', '4286347386', ...
    '4293060342', '4284946925', '4290646036', '4283846629', '4294122871', ...
    '4289359740', '4282793124', '4294704902', '4281234087', '4287019510', ...
    '4290842039', '4281692102', '4294046486', '4289294450', '4286194847', ...
    '4293637241', '4293637241', '4280684293', '4291626192', '4285560049', ...
    '4288803834', '4292504641', '4281065113', '4291454669', '4288012374', ...
    '4292388433', '4286157484'};

% ===== Path Configuration =====
ContAnatPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
fiberTrackPath = 'step2_fiberTrack';
% ===== Reference brain template ====
templateFile = 'ICBM152_adult_ICBM152_adult.T1W.nii.gz';
DSIstudioTemplate = ft_read_mri(fullfile(fiberTrackPath, templateFile));  % Load standard brain

% ===== Data Initialization =====
NTotal = 0;  % Global electrode counter
elecInfo_DSIstudio = {};  % Master cell array for DSI Studio output

% ==== Main Processing Loop =====
for nSub = 1:length(subG)
    subID = subG{nSub};
    fprintf('Processing %s (%d/%d)...\n', subID, nSub, numel(subG));

    % ---- Load Electrode Contact Data ----
    contAnatFile = fullfile(ContAnatPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']);
    load(contAnatFile);

    % ---- Parse Electrode Information ----
    elecInfo = cell(size(contAnat, 1), 3);
    for nChan = 1:size(contAnat, 1)
        chanName = contAnat.label{nChan};

        % Extract electrode name
        chanTitle = join(regexp(chanName, '\D', 'match'), '');
        chanTitle = chanTitle{1};

        % Extract contact index
        chanIndex = join(regexp(chanName, '\d', 'match'), '');
        chanIndex = str2num(chanIndex{1});

        % Convert MNI coordinates to voxel space
        MNI = contAnat.MNI{nChan};
        homogenousCoords = [MNI, 1]'; % Add homogeneous dimension
        vox_coor = (DSIstudioTemplate.transform \ homogenousCoords)';
        vox_coor = vox_coor(1:3);

        % Store parsed data
        elecInfo{nChan, 1} = chanTitle;
        elecInfo{nChan, 2} = chanIndex;
        elecInfo{nChan, 3} = vox_coor;
    end

    % ---- Calculate Electrode Geometry info for DSI studio visualization ----
    elecU = unique(elecInfo(:,1));
    for neu = 1:length(elecU)
        N = 0;  % Contact counter per electrode

        % Find start/end contacts of each electrode
        for nChan = 1:size(elecInfo, 1)
            if strcmp(elecInfo{nChan, 1}, elecU{neu})
                N = N + 1;
                if elecInfo{nChan, 2} == 1
                    elec_start = elecInfo{nChan, 3};
                end
                if elecInfo{nChan, 2} == 4
                    elec_end = elecInfo{nChan, 3};
                end
            end
        end

        % Calculate orientation vector
        elec_orientation = elec_end - elec_start;
        elec_orientation = elec_orientation./sqrt(sum(elec_orientation.^2));

        % Set DSI Studio CSV parameters
        elec_location = elec_start;
        elec_Contacts = N;
        elec_name = [subID, '-', elecU{neu}];
        elec_type = ['SEEG Electrode:', num2str(elec_Contacts), ' Contacts'];
        elec_length = 10;

        % Append to DSI Studio output matrix
        NTotal = NTotal + 1;
        elecInfo_DSIstudio{NTotal, 1} = elec_name;
        elecInfo_DSIstudio{NTotal, 2} = elec_type;
        elecInfo_DSIstudio{NTotal, 3} = [num2str(elec_location(1)), ' ', num2str(elec_location(2)), ' ', num2str(elec_location(3))];
        elecInfo_DSIstudio{NTotal, 4} = [num2str(elec_orientation(1)), ' ', num2str(elec_orientation(2)), ' ', num2str(elec_orientation(3))];
        elecInfo_DSIstudio{NTotal, 5} = num2str(elec_length(1));
        elecInfo_DSIstudio{NTotal, 6} = colorS{nSub};
    end
end

% ===== Output Results =====
% Save as CSV for DSI Studio import
writecell(elecInfo_DSIstudio, fullfile(fiberTrackPath, 'AllElectrode_dsistudio_view_diffSub.dv.csv'));
% Backup MATLAB format
save(fullfile(fiberTrackPath, 'AllElectrode_dsistudio_view.mat'), 'elecInfo_DSIstudio');