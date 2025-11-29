function IEDsinfoWin_time = f_viewIED(AX, IEDsinfo, timeWinBE, timeWin)



%% Main Function
    pointer1 = IEDsinfo(:, 1) >= timeWinBE(1) & IEDsinfo(:, 1) <=timeWinBE(2);
    pointer2 = IEDsinfo(:, 2) >= timeWinBE(1) & IEDsinfo(:, 2) <=timeWinBE(2);
    pointer3 = pointer1+pointer2;
    pointer = pointer3==2;
    IEDsinfoWin = IEDsinfo(pointer, :);
    IEDsinfoWin = IEDsinfoWin - timeWinBE(1) + 1;
    IEDsinfoWin_time = timeWin(IEDsinfoWin(:));
    IEDsinfoWin_time = reshape(IEDsinfoWin_time, size(IEDsinfoWin));
    for n = 1:size(IEDsinfoWin_time, 1)
        xline(AX, IEDsinfoWin_time(n, 1), 'Color', [0.51, 0.51, 0.51], 'LineStyle', '--', 'Tag', 'IEDbegin', 'LineWidth', 2);
        xline(AX, IEDsinfoWin_time(n, 2), 'Color', [0.51, 0.51, 0.51], 'LineStyle', '-', 'Tag', 'IEDend', 'LineWidth', 2);
    end
end