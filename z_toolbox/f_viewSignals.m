function [AX, timeWinBE, timeWin, signalsWin] = f_viewSignals(AX, data, timeWinBE, VoltageScale)
% Main function
    chanLabels = data.label;
    chanCount = length(chanLabels);
    
    winL = timeWinBE(2) - timeWinBE(1);
    time = data.time{1};
    if winL>length(time)
        timeWinBE = [timeWinBE(1), length(time)];
    end
    timeWin = time(timeWinBE(1):timeWinBE(2));

    signals = data.trial{1};
    signalsWin = signals(:, timeWinBE(1):timeWinBE(2));
    signalsWin = zscore(signalsWin')'.*0.3.*VoltageScale;
    y = repmat((1:chanCount)', 1, length(timeWin));
    signalsWin = signalsWin + y;
    colorData = lines(chanCount);
    for nChan = 1:chanCount
        line(AX, timeWin, signalsWin(nChan, :), 'Color', colorData(nChan, :), 'LineStyle', '-', 'Tag', chanLabels{nChan}, 'LineWidth', 2);
%         line(AX, timeWin, signalsWin(nChan, :), 'Color', [0.51, 0.51, 0.51], 'LineStyle', '-', 'Tag', chanLabels{nChan}, 'LineWidth', 2);
    end
    ylim([min(signalsWin(:))-0.5*std(signalsWin(1, :)), max(signalsWin(:))+0.5*std(signalsWin(end, :))]);
    xlim([timeWin(1), timeWin(end)]);
    xlabel('Time (s)');

    % channel label
    N = ceil(chanCount/20);
    AX.YTick = 1:N:chanCount;
    AX.YTickLabel = chanLabels(1:N:chanCount);

end