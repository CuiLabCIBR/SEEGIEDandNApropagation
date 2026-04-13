% Analyzes the propagation patterns of IEDs across electrode channels in 
% SEEG data, calculating both direct and indirect transfer delays and rates 
% under varying detection thresholds.
clc; clear; close all; addpath z_toolbox
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
intervalLlist = [0.05, 0.06, 0.07, 0.08, 0.09,  0.1, 0.11, 0.12, 0.13, 0.14, 0.15];
TDstr = {'5d0', '5d2', '5d4', '5d6', '5d8',   '6d0', '6d2', '6d4', '6d6', '6d8', ...
            '7d0', '7d2', '7d4', '7d6', '7d8',   '8d0', '8d2', '8d4', '8d6', '8d8', '9d0'};
fsample = 1000;
for nS = 1:length(subG) % Loop through subjects
    subID = subG{nS}; 
    load(fullfile('step3_IEDs', [subID, '_task-awake_ref-TCA_IEDs-diffTD.mat']));
    chanSize = length(IEDs{1}{1}.label);
    for nIL = 1:length(intervalLlist)
        disp("======================================")
        intervalL = intervalLlist(nIL)*fsample;
        for nTD = 1:length(TDstr)      % Loop through threshold values
                %% Extract IED cascade information
                N = 0; 
                spikeTime = {}; 
                chanList = {};
                totalRuns = length(IEDs);
                for nRun = 1:totalRuns
                    raster = IEDs{nRun}{nTD}.IEDraster{1};
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
                [~, delaysIndirect] = f_transferMatrix(spikeTime, chanList, chanSize, intervalL);
                IEDtwDelay = zeros(chanSize, chanSize);
                IEDtwRate = zeros(chanSize, chanSize);
                % Populate transfer matrices with calculated values
                for nChan1 = 1:chanSize
                        for nChan2 = 1:chanSize
                                delayTemp = delaysIndirect{nChan1, nChan2};
                                % Process indirect transfers
                                if ~isempty(delayTemp)
                                    IEDtwDelay(nChan1, nChan2) = mean(delayTemp);
                                    IEDtwRate(nChan1, nChan2) = length(delayTemp);
                                end
                        end
                end
                % Normalize transfer rates by total recording time (runs * 5 minutes)
                IEDtwRate = IEDtwRate./(totalRuns*5);

                %% display view check
                disp([IEDtwDelay(1, 2), IEDtwRate(1, 2)]); 

                %% save
                IEDtwDelay1{nS}{nIL, nTD} = IEDtwDelay;
                IEDtwRate1{nS}{nIL, nTD} = IEDtwRate;
        end
    end
    chanLabel{nS} = IEDs{1}{1}.label;
end
% save
IEDtwDelay = IEDtwDelay1;
IEDtwRate = IEDtwRate1;
savefolder = 'step3_IEDtravelingWave';
mkdir(savefolder);
save(fullfile(savefolder, 'IEDtw_Delay_Rate_diff_TD_interval.mat'), ...
    "subG", "intervalLlist", "TD", "chanLabel", "IEDtwDelay", "IEDtwRate");
