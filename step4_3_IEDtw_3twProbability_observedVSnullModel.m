% Analyzes IED propagation patterns by calculating observed probabilities 
% of 3-channel sequences and comparing them against null models with 
% randomized channel orders
clc; clear; close all; addpath z_toolbox;
intervalL = 0.1; fsample = 1000;
intervalL = intervalL*fsample;
load(fullfile('step4_IEDtw_polyfit_TD6d0', 'IEDs_meth-none_TD-6d0.mat'), 'IED6d0')
%% Construct IED cascade from raster data
% This section processes IED raster data for each subject and converts it 
% into a cascade representation showing propagation patterns
for nS = 1:length(IED6d0)
    N = 0; IEDs = IED6d0{nS};
    for nRun = 1:length(IEDs)
        % Extract IED raster data
        raster = IEDs{nRun}.IEDraster{1};
        if ~isempty(raster)
            N = N + 1;
            % Convert raster data to cascade representation
            IEDcascade{nS}{N} = f_raster2cascade(raster, intervalL);
        end
    end
end
%% Calculate observed TW3 probability
% This section calculates the probability of each 3-channel traveling wave 
% sequence based on the observed IED propagation patterns
for nS = 1:length(IEDcascade)
    IEDcascade_eachsub = IEDcascade{nS}; % Get cascade data for current subject
    N = 0; TWchanList = {};
    % Extract all channel sequences from cascade data
    for nRun = 1:length(IEDcascade_eachsub)
        IEDcascade_eachrun = IEDcascade_eachsub{nRun};
        for nC = 1:length(IEDcascade_eachrun)
            N = N + 1;
            TWchanList{N, 1} = IEDcascade_eachrun(nC).unit;
        end
    end
    % Calculate probability of each 3-element IED traveling wave
    N = 0; TW3Str = {};
    for nTW = 1:length(TWchanList)
        for nC = 1:length(TWchanList{nTW})-2
            N = N + 1;
            TW3Str{N, 1} = [num2str(TWchanList{nTW}(nC)), ...
                ' ', num2str(TWchanList{nTW}(nC+1)), ...
                ' ', num2str(TWchanList{nTW}(nC+2))];
        end
    end
    % Calculate frequency and probability of each unique 3-channel sequence
    C = tabulate(TW3Str);
    TW3probability{nS} = cell2mat(C(:, 2))./length(TW3Str);
end
%% Calculate null model TW3 probability
for nS = 1:length(IEDcascade)
% This section creates a null model by randomizing channel sequences to
% compare against observed propagation patterns (1000 permutations)
    IEDcascade_eachsub = IEDcascade{nS};
    TW3probability_null_eachsub = cell(1000, 1);
    parfor nPM = 1:1000
        disp(['Null: ', num2str(nS), '-', num2str(nPM)]);
        TW3probability_null_eachsub{nPM} = subparfor(IEDcascade_eachsub);
    end
    TW3probability_null{nS} = TW3probability_null_eachsub;
end
% Save results for statistical comparison
savefolder = 'step4_IEDtw_polyfit_TD6d0';
savename = 'IEDtw_3elements_probability_observedVSnull.mat';
save(fullfile(savefolder, savename), 'TW3probability', 'TW3probability_null');
%% Subfunction for null model calculation
function TW3probability_null = subparfor(IEDcascade_eachsub)
    % This subfunction calculates null model probabilities by randomizing
    % channel sequences while preserving sequence lengths
    N = 0; TWchanList = {};
    for nRun = 1:length(IEDcascade_eachsub)
        IEDcascade_eachrun = IEDcascade_eachsub{nRun};
        for nC = 1:length(IEDcascade_eachrun)
            N = N + 1;
            unit = IEDcascade_eachrun(nC).unit;
            idx = randperm(length(unit));
            unitPM = unit(idx);
            TWchanList{N, 1} = unitPM;
        end
    end
    % Calculate probability of each 3-element sequence in null model
    N = 0; TW3Str = {};
    for nTW = 1:length(TWchanList)
        for nC = 1:length(TWchanList{nTW})-2
            N = N + 1;
            TW3Str{N, 1} = [num2str(TWchanList{nTW}(nC)), ...
                ' ', num2str(TWchanList{nTW}(nC+1)), ...
                ' ', num2str(TWchanList{nTW}(nC+2))];
        end
    end
    % Calculate frequency and probability of each unique randomized sequence
    C = tabulate(TW3Str);
    % Return probability distribution for null model
    TW3probability_null = cell2mat(C(:, 2))./length(TW3Str);
end