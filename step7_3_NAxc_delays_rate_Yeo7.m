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
load(fullfile('step6_NAxcorr_polyfit', 'NAnon6d0IED_xcorr100ms_Rmax_delay_allSubject.mat'), ...
    'Delay_allSubject', 'xcR_allSubject');
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDtw_meth-None_TD-6d0.mat'));
DistPath = 'step2_channelDistance';
AtlasPath = 'step1_AwakeSEEG_BPfiltering_Age16Up';
Yeo7Labels = {'Vis', 'SomMot', 'DorsAttn', 'SalVentAttn',  'Limbic', 'Cont', 'DM', ...
                        'Subcortex', 'Cortex', 'White', 'none'};
for nS = 1:length(subG)
    subID = subG{nS};
    load(fullfile(DistPath, [subID,'_EuclideabDistanceMatrix.mat']), 'DistanceWorld');
    load(fullfile(AtlasPath, subID, [subID, '_ROI-3mm_contactAnatomy.mat']));
    chanLabels = IEDtw(nS).label;
    NAxcD = Delay_allSubject{nS};
    NAxcR = xcR_allSubject{nS};
    Mask = NAxcR;
    Mask = Mask - diag(diag(Mask));
    Mask(isnan(Mask)) = 0;
    zThd = 0.5 * log((1 + 0.3) ./ (1 - 0.3));
    Mask(isnan(NAxcR) | isinf(NAxcR) | NAxcR <= zThd) = 0;
    Mask(isnan(NAxcD) | isinf(NAxcD)) = 0;
    NAxcD(Mask==0) = nan;
    NAxcR(Mask==0) = nan;
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
    %% Map delay and FC matrix to Yeo7 network
    NAxcD_Yeo7 = cell(11, 11); NAxcR_Yeo7 = cell(11, 11);
    for nChan1 = 1:size(NAxcD, 1)
        for nChan2 = 1:size(NAxcD, 2)
            nIdx1 = find(ismember(Yeo7Labels, Yeo7{nChan1}));
            nIdx2 = find(ismember(Yeo7Labels, Yeo7{nChan2}));
            NAxcD_Yeo7{nIdx1, nIdx2} = [NAxcD_Yeo7{nIdx1, nIdx2}, NAxcD(nChan1, nChan2)];
            NAxcR_Yeo7{nIdx1, nIdx2} = [NAxcR_Yeo7{nIdx1, nIdx2}, NAxcR(nChan1, nChan2)];
        end
    end
    for nY1 = 1:7
        for nY2 = 1:7
            NAxcD_Yeo7Mean(nY1, nY2, nS) = mean(NAxcD_Yeo7{nY1, nY2}, 'omitnan');
            NAxcR_Yeo7Mean(nY1, nY2, nS) = mean(NAxcR_Yeo7{nY1, nY2}, 'omitnan');
        end
    end
end
saveFolder = 'step7_IEDtw_Yeo7'; mkdir(saveFolder);
saveName = 'NAxc_delays_rate_Yeo7.mat';
save(fullfile(saveFolder, saveName),  'NAxcD_Yeo7Mean', 'NAxcR_Yeo7Mean');