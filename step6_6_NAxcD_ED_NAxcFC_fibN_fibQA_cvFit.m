% Repeated cross-validation analysis of non-IED cross-correlation delay prediction.
%
% For each subject, the script:
% 1) loads subject-level non-IED cross-correlation and structural connectivity measures,
% 2) extracts valid channel pairs,
% 3) evaluates several single-predictor and full models using repeated 2-fold CV,
% 4) stores subject-level CV performance for downstream comparison.
%
% Predictors:
% - ED   : Euclidean distance
% - xcFC : non-IED cross-correlation strength (Fisher z-transformed)
% - fibN : fiber count
% - fibQA: fiber QA
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

%% Load subject-level structural connectivity results (fiber count)
load(fullfile('step2_fiberTrack',  'FiberMatrix_ncount_length_qa_allsubject.mat'));

%% Process each subject
for nS = 1:length(subG)
        subID = subG{nS}; disp(subID);

        % Load the Euclidean distance matrix between channels
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']), 'DistanceWorld');

        %% Extract transfer matrix
        twR = IEDtw(nS).indirectTransferRateMatrix;
        xcD = NAxcDelay_as{nS};
        xcFC = NAxcFC_as{nS};
        fibN = fiberCount{nS};
        fibQA = fiberqa{nS};

        % Skip subject if the fiber-count matrix is unavailable
        if isempty(fibN)
            continue; 
        end

        %% Build a validity mask for channel pairs
        % The following channel pairs are excluded:
        % 1) diagonal elements (self-self pairs)
        % 2) pairs with low transfer rate (< 0.1)
        % 3) pairs with weak XC strength (threshold corresponding to r = 0.3 in Fisher-z space)
        % 4) pairs containing NaN or Inf values in XC delay / XC strength
        Mask = twR;
        Mask = Mask - diag(diag(Mask));
        Mask(isnan(Mask)) = 0;
        Mask(Mask < 0.1) = 0;
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(xcFC) | isinf(xcFC) | xcFC <= zThd) = 0;
        Mask(isnan(xcD) | isinf(xcD)) = 0;
        Mask(fibN == 0) = 0;

        %% Extract valid observations
        ED = DistanceWorld(Mask > 0);
        xcD = xcD(Mask > 0);
        xcFC = log(xcFC(Mask > 0));
        fibN = fibN(Mask > 0);
        fibN = log(fibN);
        fibQA = fibQA(Mask > 0);

        %% Repeated 2-fold cross-validation
        % Each repetition uses a new random split generated inside fit_lsqcurvefit_cv
        stats  = {};
        parfor i = 1:1000
            disp(subID);
            stats{i} = cvfit(xcD, ED, xcFC, fibN, fibQA);
        end

        %% Store subject-level repeated CV results
        NAxcDelay_fitting{nS} = stats;
end
%% Save all analysis results
saveFolder = 'step6_NAxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op5_delay_ED_FC_fibN_fibQA_fitting.mat';
save(fullfile(saveFolder, saveName), 'NAxcDelay_fitting');

%% Local function: repeated CV for one subject
function stats = cvfit(xcD, ED, xcFC, fibN, fibQA)
    
    kfold = 2;
    %% Model 1: xcDelay ~ ED
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x./ p(3)));
    x = ED; y = xcD;
    beta0 = [0; 1; 1]; 
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_ED_R2 = R2;

     %% Model 2: xcDelay ~ IEDtwR
    modelfun = @(p, x) p(1) + p(2) .* x;
    x = xcFC; y = xcD;
    beta0 = [0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_xcFC_R2 = R2;

    %% Model 3: xcDelay ~ fibN
    modelfun = @(p, x) p(1) + p(2) .* x;
    x = fibN; y = xcD;
    beta0 = [0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_fibN_R2 = R2;

    %% Model 4: xcDelay ~ fibQA
    modelfun = @(p, x) p(1) + p(2) .* x;
    x = fibQA; y = xcD;
    beta0 = [0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_fibQA_R2 = R2;

    %% Model 5: xcDelay ~ ED + xcFC + fibN + fibQA
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2) + p(5) .* x(:, 3) + p(6) .* x(:, 4);
    x = [ED, xcFC, fibN, fibQA];
    y = xcD;
    beta0 = [0; 1; 1; 0; 0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_R2 = R2;
    
    %% Model 6: xcDelay ~ xcFC + fibN + fibQA
    modelfun = @(p, x) p(1) + p(2) .* x(:, 1) + p(3) .* x(:, 2) + p(4) .* x(:, 3);
    x = [xcFC, fibN, fibQA];
    y = xcD;
    beta0 = [0; 0; 0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_noED_R2 = R2;

    %% Model 7: xcDelay ~ ED + fibN + fibQA
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2) + p(5) .* x(:, 3);
    x = [ED, fibN, fibQA];
    y = xcD;
    beta0 = [0; 1; 1; 0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_noFC_R2 = R2;

    %% Model 8: xcDelay ~ ED + xcFC + fibQA
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2) + p(5) .* x(:, 3);
    x = [ED, xcFC, fibQA];
    y = xcD;
    beta0 = [0; 1; 1; 0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_noFibN_R2 = R2;

    %% Model 9: xcDelay ~ ED + xcFC + fibN
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2) + p(5) .* x(:, 3);
    x = [ED, xcFC, fibN];
    y = xcD;
    beta0 = [0; 1; 1; 0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_noFibQA_R2 = R2;
end