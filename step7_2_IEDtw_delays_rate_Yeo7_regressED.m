% Calculates the comparison of IED propagation delays within the Yeo 7 network 
% and between networks, after removing those under 20mm and applying distance
% regression.
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
DistPath = 'step2_channelDistance';
AtlasPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
Yeo7Labels = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn',  'Limbic', 'Cont', 'DM', ...
                        'Subcortex', 'Cortex', 'White', 'none'};
% ==== Main Analysis Loop ====
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtw_meth-None_TD-6d0.mat'));
for nS = 1:length(subG)
    subID = subG{nS};
    % Load Euclidean distance matrix between channels
    load(fullfile(DistPath, [subID,'_EuclideabDistanceMatrix.mat']), 'DistanceWorld');
    % Load anatomical contact information
    load(fullfile(AtlasPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));
    %% Extract IED conduction delays and rate
    chanLabels = IEDtw(nS).label;
    delays = IEDtw(nS).indirectTransferDelayMatrix;
    rate = IEDtw(nS).indirectTransferRateMatrix;
    Mask = rate;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    Mask(Mask<0.1) = 0;
    Mask(DistanceWorld<20) = 0;
    delays(Mask==0) = nan;
    DistanceWorld(Mask==0) = nan;
    rate(Mask==0) = nan;
    %% Regress Euclidean distance
    x = DistanceWorld(Mask>0);
    y = delays(Mask>0);
    polyModel = f_polyfit(x, y, 1);
    delays(Mask>0) = polyModel.residual + polyModel.p(end);
    y = rate(Mask>0);
    polyModel = f_polyfit(x, y, 1);
    rate(Mask>0) = polyModel.residual + polyModel.p(end);
    %% the Yeo7 network label of each channel
    Yeo7 = {};
    for nChan = 1:length(chanLabels)
            Lia = ismember(contAnat.label, chanLabels{nChan});
            Yeo7Temp = contAnat.Yeo7(Lia);
            if strcmp(Yeo7Temp, 'none')
                CWS = contAnat.CWS(Lia);
                Yeo7{nChan, 1} = CWS{1};
            else
                Yeo7{nChan, 1} = Yeo7Temp{1};
            end
    end
    %% Calculate delay and distance in Yeo7 level
    delaysYeo7 = cell(11, 11); 
    rateYeo7 = cell(11, 11);
    for nChan1 = 1:size(delays, 1)
        for nChan2 = 1:size(delays, 2)
            Lia1 = find(ismember(Yeo7Labels, Yeo7{nChan1}));
            Lia2 = find(ismember(Yeo7Labels, Yeo7{nChan2}));
            delaysYeo7{Lia1, Lia2} = [delaysYeo7{Lia1, Lia2}, delays(nChan1, nChan2)];
            rateYeo7{Lia1, Lia2} = [rateYeo7{Lia1, Lia2}, rate(nChan1, nChan2)];
        end
    end
    % Compute Network-Level Averages
    for nY1 = 1:7
        for nY2 = 1:7
            delaysYeo7Mean(nY1, nY2, nS) = mean(delaysYeo7{nY1, nY2}, 'omitnan');
            rateYeo7Mean(nY1, nY2, nS) = mean(rateYeo7{nY1, nY2}, 'omitnan');
        end
    end
end
saveFolder = 'step7_IEDtw_Yeo7'; mkdir(saveFolder);
saveName = 'IEDtw_delays_rate_Yeo7_rED.mat';
save(fullfile(saveFolder, saveName), 'delaysYeo7Mean', 'rateYeo7Mean');