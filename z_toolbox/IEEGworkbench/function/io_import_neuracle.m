function [data_fieldtrip, event] = io_import_neuracle(foldname)
    addpath_eeglab;
    eegdata = pop_importNeuracle( {'data.bdf', 'evt.bdf'}, foldname);
    data_fieldtrip.fsample = eegdata.srate;
    data_fieldtrip.sampleinfo = [1, eegdata.pnts];
    data_fieldtrip.trial{1} = eegdata.data;
    data_fieldtrip.time{1} = eegdata.times;
    chanlocs = struct2table(eegdata.chanlocs);
    data_fieldtrip.label = chanlocs.labels;
    event = eegdata.event;
    data_fieldtrip.event = event;
end

