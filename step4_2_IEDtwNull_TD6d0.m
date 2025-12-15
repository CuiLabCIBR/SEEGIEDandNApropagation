% Constructs a null model for IED traveling wave propagation by performing 
% 1000 permutation tests with randomly shuffled channel sequences to 
% generate statistical baselines for direct and indirect transfer 
% matrices of delays and rates
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
% Construct IED traveling wave Null model
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDs_meth-none_TD-6d0.mat'), 'IED6d0')
for nS = 1:length(subG)
        subID = subG{nS};
        IEDs = IED6d0{nS};
        Size = length(IEDs{1}.label);
        % Initialize null model matrices for 1000 permutations:
        IEDtwD_DNull = zeros(Size, Size, 1000, 'single'); % - Direct transfer delay matrix
        IEDtwR_DNull = zeros(Size, Size, 1000, 'single'); % - Direct transfer rate matrix  
        IEDtwD_INull = zeros(Size, Size, 1000, 'single'); % - Indirect transfer delay matrix
        IEDtwR_INull = zeros(Size, Size, 1000, 'single'); % - Indirect transfer rate matrix
        % Parallel loop for 1000 permutation tests
        parfor nPM = 1:1000
            disp([subID, ' ', num2str(nPM)]);
            [IEDtfDelayD, IEDtfRateD, IEDtfDelayI, IEDtfRateI] = subparfor(IEDs, Size);
            IEDtwD_DNull(:, :, nPM) = IEDtfDelayD;
            IEDtwR_DNull(:, :, nPM) = IEDtfRateD;
            IEDtwD_INull(:, :, nPM) = IEDtfDelayI;
            IEDtwR_INull(:, :, nPM) = IEDtfRateI;
        end
        % Create output structure with results
        IEDtwNull.label = IEDs{1}.label;
        IEDtwNull.directTransferDelayMatrix = IEDtwD_DNull;
        IEDtwNull.directTransferRateMatrix = IEDtwR_DNull;
        IEDtwNull.indirectTransferDelayMatrix = IEDtwD_INull;
        IEDtwNull.indirectTransferRateMatrix = IEDtwR_INull;
        % Save results for current subject
        saveFolder = fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtwNull');
        mkdir(saveFolder);
        saveName = [subID, '_IEDtwNull_meth-None_TD-6d0.mat'];
        save(fullfile(saveFolder, saveName), "IEDtwNull");
end
%% SUBPARFOR Calculate IED transfer matrices for null model analysis
function [IEDtwD_D, IEDtwR_D, IEDtwD_I, IEDtwR_I] = subparfor(IEDs, Size)
        intervalL = 0.1;
        fsample = 1000;
        intervalL = intervalL*fsample;
        % Extract and process IED cascades with random permutation
        N = 0; spikeTime = {}; chanList = {};
        for nRun = 1:length(IEDs)
            raster = IEDs{nRun}.IEDraster{1};
            if ~isempty(raster)
                Cascade = f_raster2cascade(raster, intervalL);
                for nC = 1:length(Cascade)
                    N = N + 1;
                    spikeTime{N} = Cascade(nC).time;
                    unit = Cascade(nC).unit;
                    idx = randperm(length(unit));
                    unitPM = unit(idx);
                    chanList{N} = unitPM;
                end
            end
        end
        % Calculate transfer matrices using permuted data
        [delayDirect, delayIndirect] = f_transferMatrix(spikeTime, chanList, Size, intervalL);
        IEDtwD_D = zeros(Size, Size); IEDtwR_D = zeros(Size, Size);
        IEDtwD_I = zeros(Size, Size); IEDtwR_I = zeros(Size, Size);
        for nChan1 = 1:Size
            for nchan2 = 1:Size
                delayTemp1 = delayDirect{nChan1, nchan2};
                delayTemp2 = delayIndirect{nChan1, nchan2};
                if ~isempty(delayTemp1)
                    IEDtwD_D(nChan1, nchan2) = mean(delayTemp1);
                    IEDtwR_D(nChan1, nchan2) = length(delayTemp1);
                end
                if ~isempty(delayTemp2)
                    IEDtwD_I(nChan1, nchan2) = mean(delayTemp2);
                    IEDtwR_I(nChan1, nchan2) = length(delayTemp2);
                end
            end
        end
        % Normalize transfer rates by number of runs and trials
        IEDtwR_D = IEDtwR_D./(nRun*5);
        IEDtwR_I = IEDtwR_I./(nRun*5);
end