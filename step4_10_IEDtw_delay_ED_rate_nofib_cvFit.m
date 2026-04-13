clc; clear; close all; addpath z_toolbox;
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
% Load subject-wise structural connectivity metrics
% fiberCount{nS} : fiber count matrix
% fiberLength{nS}: fiber length matrix
% fiberqa{nS}: fiber QA matrix
load(fullfile('step2_fiberTrack', 'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
% Load subject-wise IED traveling wave measures
% IEDtw(nS).indirectTransferDelayMatrix : delay matrix
% IEDtw(nS).indirectTransferRateMatrix  : rate matrix
load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op2_IEDtw_TD-6d0_IT-100ms.mat'));

%% Store repeated CV results for all subjects
for nS = 1:length(subG)
    subID = subG{nS};
    % Load Euclidean distance matrix for current subject
    load(fullfile('step2_channelDistance', [subID,'_EuclideabDistanceMatrix.mat']));
    % Extract structural connectivity matrices for current subject
    fibN = fiberCount{nS};
    % Skip subject if fiber count matrix is empty
    if isempty(fibN)
        continue;
    end

    % Extract IED traveling wave matrices
    IEDtwD = IEDtw(nS).indirectTransferDelayMatrix;
    IEDtwR = IEDtw(nS).indirectTransferRateMatrix;

    %% Build analysis mask:
    % 1. remove diagonal
    % 2. replace NaN with 0
    % 3. keep only edges with traveling-wave rate >= 0.1
    % 4. remove edges with zero fiber count
    Mask = IEDtwR;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    Mask(Mask < 0.1) = 0;
    Mask(fibN > 0 ) = 0;

    % Vectorize edge-wise measures using the final mask
    ED = DistanceWorld(Mask>0);
    IEDtwD = IEDtwD(Mask>0);
    IEDtwR = IEDtwR(Mask>0);
    IEDtwR = log(IEDtwR);
    
    %% Repeated 2-fold cross-validation
    % Each repetition uses a new random partition generated inside fit_lsqcurvefit_cv
    stats  = {};
    parfor i = 1:1000
        disp(subID);
        stats{i} = cvfit(IEDtwD, ED, IEDtwR);
    end
    
    % Save subject-level repeated CV results
    IEDtwDelay_fitting_nofib{nS} = stats;
end
%%
saveFolder = 'step4_IEDtw_TD-6d0_IT-100ms_fitting';
saveName = 'op10_IEDtw_delay_ED_rate_nofib_fitting.mat';
save(fullfile(saveFolder, saveName), 'IEDtwDelay_fitting_nofib');

%%
function stats = cvfit(IEDtwD, ED, IEDtwR)
    
    kfold = 2;
    %% Model 1: IEDtwD ~ ED
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x./ p(3)));
    x = ED;  y = IEDtwD;
    beta0 = [0; 1; 1]; 
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_ED_R2 = R2;

     %% Model 2: IEDtwD ~ IEDtwR
    modelfun = @(p, x) p(1) + p(2) .* x;
    x = IEDtwR; y = IEDtwD;
    beta0 = [0; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_rate_R2 = R2;
 
    %% Model 3: Full model 
    % IEDtwD ~ ED + IEDtwR
    modelfun = @(p, x) p(1) + p(2) .* (1 - exp(-x(:, 1) ./ p(3))) + p(4) .* x(:, 2);
    x = [ED, IEDtwR]; y = IEDtwD;
    beta0 = [0; 1; 1; 0];
    [R2, ~, ~, ~] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold);
    disp(['R2 = ', num2str(R2)]);
    stats.delay_R2 = R2;
end