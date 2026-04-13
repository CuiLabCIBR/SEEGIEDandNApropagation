% Process subject-level IED cross-correlation results to compute
% mean inter-channel time delays and mean Fisher z-transformed
% cross-correlation values, after filtering invalid observations.
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
% Preallocate cell arrays for all subjects
IEDxcDelay_as = cell(length(subG), 1);
IEDxcFC_as = cell(length(subG), 1);

%% Process each subject
for nS = 1:length(subG)
        subID = subG{nS}; disp(subID);
        % Load cross-correlation results for the current subject
        load(fullfile('step5_IEDxcorr',  [subID, '_6d0IED_xcorr100ms_Rmax_delay.mat']));
        % Determine the number of channels
        chanCount = size(IEDxcDelay, 1);
        % Initialize output matrices for the current subject
        xcDelay_mean = zeros(chanCount, chanCount);
        xcFC_mean = zeros(chanCount, chanCount);

        % Process each channel pair
        for nChan1 = 1:chanCount
                for nChan2 = 1:chanCount

                        % Extract delay values and correlation coefficients
                        A = abs(IEDxcDelay{nChan1, nChan2});
                        B = IEDxcR{nChan1, nChan2};

                        % Remove invalid observations:
                        % 1) delays >= 99
                        % 2) correlations <= 0.3
                        B(A >= 99) = []; A(A >= 99) = [];
                        A(B <= 0.3) = []; B(B <= 0.3) = [];

                        % Apply Fisher z-transform to correlation coefficients
                        z = 0.5 * log((1 + B) ./ (1 - B));

                        % Compute mean delay for the current channel pair
                        xcDelay_mean(nChan1, nChan2) = mean(A, 'omitnan');

                        % Compute mean Fisher z-transformed correlation
                        xcFC_mean(nChan1, nChan2) = mean(z, 'omitnan');
                end
        end
        % Set diagonal elements (self-connections) to zero
        xcDelay_mean = xcDelay_mean - diag(diag(xcDelay_mean));
        xcFC_mean = xcFC_mean - diag(diag(xcFC_mean));
        % Store subject-level results
        IEDxcDelay_as{nS} = xcDelay_mean;
        IEDxcFC_as{nS} = xcFC_mean;
end
%% Save aggregated results
saveFolder = 'step5_IEDxcorr_TD-6d0_IT-100ms_fitting';
mkdir(saveFolder);
saveName = 'op1_IEDxcorr_delay_FC_allSubject.mat';
save(fullfile(saveFolder, saveName),  'IEDxcDelay_as',  'IEDxcFC_as');