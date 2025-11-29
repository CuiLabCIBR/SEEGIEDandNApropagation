function travelingWave = f_IED_travelingWave(IEDs, intervalL)

%%% @author: Longzhou Xu
%%% @version: | Last modified: 2025-04-28
%%
    % Extract channel labels and sampling frequency
    label = IEDs.channelLabel;
    fsample = IEDs.fsample;
    intervalL = intervalL*fsample;

    % Initialize IED event list and output structure
    IEDlist = IEDs.IEDlist;
    trialCount = length(IEDlist);
    travelingWave = [];

    % Process each trial
    for nT = 1:trialCount
        IEDlist1 = IEDlist{nT};  % IED events for current trial
       
        % Construct IED raster: [negPeakTime, channelIndex]
        raster = []; channelList = {};
        for nIED = 1:size(IEDlist1, 2)
            for nChan = 1:length(label) % Match IED to channel label
                if strcmp(IEDlist1(nIED).channelLabel, label{nChan})
                    % Store negative peak time and channel index
                    raster(nIED, 1) = IEDlist1(nIED).negPeakTime;
                    raster(nIED, 2) = nChan;
                    channelList{nIED} =  IEDlist1(nIED).channelLabel;
                end
            end
        end

        % Detect cascades of temporally proximate IEDs
        Cascade = f_raster2cascade(raster, intervalL);

        % Store results for current trial
        travelingWave(nT).raster = raster;
        travelingWave(nT).Cascade = Cascade;
    end
end