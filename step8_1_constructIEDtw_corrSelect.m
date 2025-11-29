% Analyzes IED propagation patterns across brain regions by calculating 
% indirect transfer delays and rates using spatiotemporal correlation 
% filtering.
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
% ===== IED Traveling Wave Parameters =====
intervalL = 0.1; fsample = 1000;
intervalL = intervalL*fsample;
% ===== IED Detection Threshold Parameters =====
TDstr = {'5d0', '5d2', '5d4', '5d6', '5d8', ...
        '6d0', '6d2', '6d4', '6d6', '6d8', ...
        '7d0', '7d2', '7d4', '7d6', '7d8', ...
        '8d0', '8d2', '8d4', '8d6', '8d8', '9d0'};
TD = 5:0.2:9;
% ==== Main Processing Loop =====
for nTD = 1:length(TD)       % Loop through threshold values
    for nS = 1:length(subG)  % Loop through subjects
        subID = subG{nS}; disp([TDstr{nTD}, ' ', subID]);
        % Load IED data and Euclidean distance matrix for current subject
        load(fullfile('step3_IEDs', [subID, '_task-awake_ref-TCA_IEDs-diffTD.mat']));
        load(fullfile('step2_channelDistance', [subID, '_EuclideabDistanceMatrix.mat']));
        %% Extract IED cascade information from all recording runs
        Size = length(IEDs{1}{1}.label);
        totalRuns = length(IEDs);
        CascadeAll = [];
        % Process each recording run to extract IED cascades
        for nRun = 1:totalRuns
            raster = IEDs{nRun}{nTD}.IEDraster{1};
            if ~isempty(raster)
                % Convert raster format to cascade format
                Cascade = f_raster2cascade(raster, intervalL);
                CascadeAll = [CascadeAll, Cascade];
            end
        end
        %% Identify IED traveling wave using Spearman correlation
        % Process each cascade to find valid triplets
        r = []; spikeTime = {}; chanList = {};
        for nCa = 1:length(CascadeAll)
            IEDtwChan = CascadeAll(nCa).unit;
            IEDtwTime = CascadeAll(nCa).time;
            if length(IEDtwChan) >= 4
                Tlist = abs(IEDtwTime(2:end) - IEDtwTime(1));
                Dlist = DistanceWorld(IEDtwChan(2:end), IEDtwChan(1));
                r(nCa) = corr(Tlist, Dlist, 'type', 'Spearman');
            else
                r(nCa) = 0;
            end
            if r(nCa)>0.3
                spikeTime{end+1} = IEDtwTime;
                chanList{end+1} = IEDtwChan;
            end
        end
        %% ===== Calculate Transfer Matrices =====
        IEDtwD = cell(Size, Size);
        for ntw = 1:length(spikeTime)
            for n = 2:length(chanList{ntw})
                IEDtwD{chanList{ntw}(1), chanList{ntw}(n)} = ...
                    [IEDtwD{chanList{ntw}(1), chanList{ntw}(n)}, ...
                    spikeTime{ntw}(n) - spikeTime{ntw}(1)];
            end
        end
        IEDtransferDelays = zeros(Size, Size);
        IEDtransferRate = zeros(Size, Size);
        for n1 = 1:size(IEDtwD, 1)
            for n2 = 1:size(IEDtwD, 2)
                IEDtransferDelays(n1, n2) = mean(IEDtwD{n1, n2}, 'omitnan');
                IEDtransferRate(n1, n2) = length(IEDtwD{n1, n2});
            end
        end
        IEDtransferRate = IEDtransferRate./(totalRuns*5);
        %% Store Results
        IEDtw(nS).label = IEDs{1}{1}.label;
        IEDtw(nS).travelingWaveSpikeTime = spikeTime;
        IEDtw(nS).travelingWaveChanList = chanList;
        IEDtw(nS).indirectTransferDelayMatrix = IEDtransferDelays;
        IEDtw(nS).indirectTransferRateMatrix = IEDtransferRate;
    end
    % Save Results for Current Threshold
    saveFolder = 'step8_IEDtravelingWave_corrSelect'; mkdir(saveFolder);
    saveName = ['IEDtw_meth-corrSelect_TD-', TDstr{nTD}, '.mat'];
    save(fullfile(saveFolder, saveName), "IEDtw");
end