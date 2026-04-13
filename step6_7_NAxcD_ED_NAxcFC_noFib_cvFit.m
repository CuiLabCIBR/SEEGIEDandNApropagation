% Repeated cross-validation analysis of non-IED cross-correlation delay prediction.
%
% For each subject, the script:
% 1) loads subject-level non-IED cross-correlation measures,
% 2) extracts valid channel pairs,
% 3) evaluates single-predictor and combined models using repeated 2-fold CV,
% 4) compares the contribution of Euclidean distance (ED) and non-IED
%    cross-correlation strength (xcFC) by using shuffled-predictor controls.
clc; clear; close all;  addpath z_toolbox;
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

%% Load subject-level traveling-wave results
 load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op2_IEDtw_TD-6d0_IT-100ms.mat'));

%% Load subject-level cross-correlation results
load(fullfile('step6_NAxcorr_TD-6d0_IT-100ms_fitting',  'op1_NAxcorr_delay_FC.mat'));

%% Load subject-level structural connectivity metrics
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'));

%% Process each subject
for nS = 1:length(subG)
        subID = subG{nS}; disp(subID);

        % Load the Euclidean distance matrix between channels
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']), 'DistanceWorld');
        
        %% Extract subject-specific matrices
        twR = IEDtw(nS).indirectTransferRateMatrix;
        xcD = NAxcDelay_as{nS};
        xcFC = NAxcFC_as{nS};
        fibN = fiberCount{nS};

        % Skip subject if the fiber-count matrix is unavailable
        if isempty(fibN)
            continue; 
        end

        %% Build a mask to retain valid channel pairs
        % Exclude:
        % 1) diagonal elements
        % 2) weak XC strength values
        % 3) invalid XC delay values
        % 4) channel pairs with zero fiber count
        % 5) invalid fiber-QA values
        Mask = ones(size(xcD));
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask < 0.1) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(xcFC) | isinf(xcFC) | xcFC <= zThd) = 0;
        Mask(isnan(xcD) | isinf(xcD)) = 0;
        Mask(fibN == 1) = 0;
        
        %% Extract valid observations
        ED = DistanceWorld(Mask > 0);
        xcD = xcD(Mask > 0);
        xcFC = log(xcFC(Mask > 0));

        %% Repeated 2-fold cross-validation
        % Each repetition uses a new random split generated inside fit_lsqcurvefit_cv
        stats  = {};
        parfor i = 1:1000
            disp(subID);
            stats{i} = cvfit(xcD, ED, xcFC);
        end
        
        %% Store subject-level repeated CV results
        NAxcDelay_fitting_nofib{nS} = stats;
end
%% Save all subject-level CV results
saveFolder = 'step6_NAxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op6_delay_ED_FC_noFib_fitting.mat';
save(fullfile(saveFolder, saveName), 'NAxcDelay_fitting_nofib');

%% Local function: repeated CV for one subject
function stats = cvfit(xcD, ED, xcFC)
    kfold = 2;
    %% Model 1: xcDelay ~ ED
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x./ p(3)));
    x = ED; 
    y = xcD;
    beta0 = [1; 1; 1]; 
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_ED_R2 = R2;

     %% Model 2: xcDelay ~ IEDtwR
    modelfun = @(p, x) p(1) + p(2) .* x;
    x = xcFC; 
    y = xcD;
    beta0 = [1; -1];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_xcFC_R2 = R2;

    %% Model 3: Full model 
    % xcDelay ~ ED + xcFC
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2) ;
    x = [ED, xcFC];
    y = xcD;
    beta0 = [1; 1; 1; -1];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_R2 = R2;
end