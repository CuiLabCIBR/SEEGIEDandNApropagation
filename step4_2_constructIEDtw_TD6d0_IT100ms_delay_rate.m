% Analyzes the propagation patterns of IEDs across electrode channels in 
% SEEG data, calculating both direct and indirect transfer delays and rates 
% under varying detection thresholds.
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
intervalL = 0.1;
fsample = 1000;
intervalL = intervalL*fsample;
load(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op1_IEDs_TD-6d0.mat'))
% Loop through subjects
for nS = 1:length(subG) 
        subID = subG{nS};
        IEDs = IED6d0{nS};

        %% Extract IED cascade information
        Size = length(IEDs{1}.label);
        N = 0; 
        spikeTime = {}; 
        chanList = {};
        totalRuns = length(IEDs);
        for nRun = 1:totalRuns
            raster = IEDs{nRun}.IEDraster{1};
            if ~isempty(raster)
                % Convert raster data to cascade format
                Cascade = f_raster2cascade(raster, intervalL);
                % Extract time and channel information from each cascade
                for nC = 1:length(Cascade)
                    N = N + 1;
                    spikeTime{N} = Cascade(nC).time;
                    chanList{N} = Cascade(nC).unit;
                end
            end
        end

        %% Calculate mean delays and count rates for each channel pair
        % Calculate Transfer Matrices
        [delaysDirect, delaysIndirect] = f_transferMatrix(spikeTime, chanList, Size, intervalL);
        IEDtransferDelaysD = zeros(Size, Size);
        IEDtransferRateD = zeros(Size, Size);
        IEDtransferDelaysI = zeros(Size, Size);
        IEDtransferRateI = zeros(Size, Size);
        % Populate transfer matrices with calculated values
        for nChan1 = 1:Size
                for nChan2 = 1:Size
                        delayTemp1 = delaysDirect{nChan1, nChan2};
                        delayTemp2 = delaysIndirect{nChan1, nChan2};
                        % Process direct transfers
                        if ~isempty(delayTemp1)
                            IEDtransferDelaysD(nChan1, nChan2) = mean(delayTemp1);
                            IEDtransferRateD(nChan1, nChan2) = length(delayTemp1);
                        end
                        % Process indirect transfers
                        if ~isempty(delayTemp2)
                            IEDtransferDelaysI(nChan1, nChan2) = mean(delayTemp2);
                            IEDtransferRateI(nChan1, nChan2) = length(delayTemp2);
                        end
                end
        end
        % Normalize transfer rates by total recording time (runs * 5 minutes)
        IEDtransferRateD = IEDtransferRateD./(totalRuns*5);
        IEDtransferRateI = IEDtransferRateI./(totalRuns*5);

        %% Store results in structure
        IEDtw(nS).label = IEDs{1}.label;
        IEDtw(nS).travelingWaveSpikeTime = spikeTime;
        IEDtw(nS).travelingWaveChanList = chanList;
        IEDtw(nS).directTransferDelayMatrix = IEDtransferDelaysD;
        IEDtw(nS).directTransferRateMatrix = IEDtransferRateD;
        IEDtw(nS).indirectTransferDelayMatrix = IEDtransferDelaysI;
        IEDtw(nS).indirectTransferRateMatrix = IEDtransferRateI;
end
save(fullfile('step4_IEDtw_TD-6d0_IT-100ms_fitting', 'op2_IEDtw_TD-6d0_IT-100ms.mat'), "IEDtw");