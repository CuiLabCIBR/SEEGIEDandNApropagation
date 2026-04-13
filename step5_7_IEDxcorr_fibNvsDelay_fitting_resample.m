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
%% Load subject-level structural connectivity results (fiber count)
load(fullfile('step2_fiberTrack',  'FiberMatrix_ncount_length_qa_allsubject.mat'), 'fiberCount');
%% Load subject-level cross-correlation results
load(fullfile('step5_IEDxcorr_TD-6d0_IT-100ms_fitting',  'op1_IEDxcorr_delay_FC_allSubject.mat'));
%% Process each subject
for nS = 1:length(subG)
        subID = subG{nS}; disp(subID);
        % Extract subject-specific data
        xcD = IEDxcDelay_as{nS};
        xcFC = IEDxcFC_as{nS};
        fibN = fiberCount{nS};
        % Skip subject if fiber-count matrix is unavailable
        if isempty(fibN)
            continue; 
        end

        %% Construct a mask for valid channel pairs
        % Exclude:
        % 1) diagonal elements
        % 2) weak XC strength values
        % 3) invalid XC delay values
        % 4) channel pairs with zero fiber count
        Mask = xcFC;
        Mask = Mask - diag(diag(Mask));
        zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
        Mask(isnan(xcFC) | isinf(xcFC) | xcFC <= zThd) = 0;
        Mask(isnan(xcD) | isinf(xcD)) = 0;
        Mask(fibN == 0) = 0;

        %% Extract valid observations
        xcD = xcD(Mask>0);
        fibN = log(fibN(Mask>0));
        xcD = xcD(:);
        fibN = fibN(:);

        %% --- Uniform subsampling across xcD bins, then correlation fibN vs xcD ---
        % Parameters you can tune
        nBins   = 20;      % number of xcD bins (equal-width)
        nRepeat = 1000;     % how many random re-samplings
        minPerBinTarget = 20;  % optional: ensure at least this many per bin; otherwise reduce bins
        % Helper: choose bin number so that min count per non-empty bin is not too small
        % (optional but recommended)
        for tryBins = nBins:-1:5
            edges = linspace(min(xcD), max(xcD), tryBins+1);
            binID = discretize(xcD, edges);
            counts = accumarray(binID(~isnan(binID)), 1, [tryBins 1], @sum, 0);
            validBins = find(counts > 0);
            if isempty(validBins); continue; end
            nPerBin = min(counts(validBins));
            if nPerBin >= minPerBinTarget
                nBins = tryBins;
                break;
            end
        end

        % Recompute with final nBins
        edges = linspace(min(xcD), max(xcD), nBins+1);
        binID = discretize(xcD, edges);
        counts = accumarray(binID(~isnan(binID)), 1, [nBins 1], @sum, 0);
        validBins = find(counts > 0);
        nPerBin = min(counts(validBins));
        
        %%
        rhoRep = nan(nRepeat, 1);
        pRep   = nan(nRepeat, 1);
        for r = 1:nRepeat
            sel = false(size(xcD));
            for b = validBins(:)'
                idxB = find(binID == b);
                % sample exactly nPerBin from each bin
                pick = idxB(randperm(numel(idxB), nPerBin));
                sel(pick) = true;
            end
            xcD_u  = xcD(sel);
            fibN_u = fibN(sel);
            % Spearman is usually safer for skewed/heteroscedastic relationships
            [rhoRep(r), pRep(r)] = corr(fibN_u, xcD_u, 'Type', 'Spearman', 'Rows', 'complete');
        end

        % Summarize for this subject
        rho_mean = mean(rhoRep, 'omitnan');
        rho_std = std(rhoRep, 'omitnan');
        
        % Store outputs (create arrays outside loop)
        rho_uniform_mean(nS) = rho_mean;
        rho_uniform_std(nS) = rho_std;
end
saveFolder = 'step5_IEDxcorr_TD-6d0_IT-100ms_fitting';
saveName = 'op6_fibN-vs-xcD_fitting_uniformResample.mat';
save(fullfile(saveFolder, saveName), 'rho_uniform_mean', 'rho_uniform_std');