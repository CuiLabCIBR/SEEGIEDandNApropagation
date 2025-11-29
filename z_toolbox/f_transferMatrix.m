function  [delays_direct, delays_indirect] = f_transferMatrix(spikeTime, chanList, Size, interval)
% Computes the transfer matrix between channels based on spike cascade data.
% The transfer matrix quantifies directional relationships between channels 
% (e.g., neurons or EEG electrodes) in terms of event counts and time delays.
%
% Inputs:
%   Cascade : Struct array with fields:
%             - time : Array of spike timestamps for each cascade.
%             - unit : Array of channel indices for each spike in the cascade.
%   Size    : Scalar, size of the output matrix (e.g., number of channels).
%%% @author Longzhou Xu
%%% @version | Last modified: 2025-05-08
%% Main Function
        % Direct IED transfer 
        delays_direct = cell(Size, Size);  % Cumulative time delays (i -> j)
        for nTW = 1:length(spikeTime)      % Loop through each cascade
            spikeTime1 = spikeTime{nTW};   % Spike timestamps
            chanList1 = chanList{nTW};    % Channel indices
            % Iterate over consecutive spike pairs (S -> T)
            for n = 1:length(chanList1)-1
                S = chanList1(n);        % Source channel
                T = chanList1(n+1);      % Target channel
                % transfer delay sum
                delays_direct{S, T} = [delays_direct{S, T}, spikeTime1(n+1)-spikeTime1(n)];  % Add delay
            end
        end
    
        % Indirect IED transfer
        delays_indirect = cell(Size, Size);
        for nTW = 1:length(spikeTime)      % Loop through each cascade
            spikeTime1 = spikeTime{nTW};   % Spike timestamps
            chanList1 = chanList{nTW};    % Channel indices
            % Iterate over consecutive spike pairs (S -> T)
            for n1 = 1:length(chanList1)-1
                for n2 = n1+1:length(chanList1)
                    S = chanList1(n1);        % Source channel
                    T = chanList1(n2);      % Target channel
                    I = spikeTime1(n2)-spikeTime1(n1);
                    % transfer delay sum
                    if I <= interval
                        delays_indirect{S, T} = [delays_indirect{S, T}, I];  % Add delay
                    end
                end
            end
        end
end