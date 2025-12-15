% Processes cross-correlation data during non-IED EEG periods, calculating
% mean delays and Fisher z-transformed correlation coefficients between 
% channel pairs while excluding invalid data points and self-connections.
clc;clear;close all;
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
% Initialize cell arrays to store results for all subjects
Delay_allSubject = cell(length(subG), 1);
xcR_allSubject = cell(length(subG), 1);
%% Process each subject
for nS = 1:length(subG)
    subID = subG{nS}; disp(subID);
    % Load cross-correlation data for current subject and time window
    load(fullfile('step6_NAxcorr', [subID, '_NAnon6d0IED_xcorr100ms_Rmax_delay.mat']));
    % Get number of channels from loaded data
    chanCount = size(NAxcDelay_allRun, 1);
    % Initialize matrices for current subject
    Delay = zeros(chanCount, chanCount);
    DelayNorm = zeros(chanCount, chanCount); 
    xcR = zeros(chanCount, chanCount);
    % Process each channel pair
    for nChan1 = 1:chanCount
        for nChan2 = 1:chanCount
            % Extract delay and correlation values for current channel pair
            A = squeeze(NAxcDelay_allRun(nChan1, nChan2, :));
            B = squeeze(NAxcR_allRun(nChan1, nChan2, :));
            % For 100ms window: remove delays ≥99ms
            B(A>=99) = []; A(A>=99) = [];
            % Remove non-positive correlations and corresponding delays
            A(B<=0.3) = []; B(B<=0.3) = [];
            % Apply Fisher z-transform to correlation coefficients
            z = 0.5 * log((1 + B) ./ (1 - B));
            % Calculate mean delay (unweighted)
            Delay(nChan1, nChan2) = mean(A, 'omitnan');
            % Calculate mean Fisher z-transformed correlation
            xcR(nChan1, nChan2) = mean(z, 'omitnan');
        end
    end
    % Set diagonal elements to zero (self-connections)
    Delay = Delay - diag(diag(Delay));
    xcR = xcR - diag(diag(xcR));
    % Store results for current subject
    Delay_allSubject{nS} = Delay;
    xcR_allSubject{nS} = xcR;
end
% Save results
saveFolder = 'step6_NAxcorr_polyfit';
mkdir(saveFolder);
saveName = 'NAnon6d0IED_xcorr100ms_Rmax_delay_allSubject.mat';
save(fullfile(saveFolder, saveName), 'Delay_allSubject',  'xcR_allSubject');