function [IEDrasterWin_time, IEDrasterWin_value] = f_viewIEDnegPeak(AX, IEDraster, timeWinBE, timeWin, signalsWin)
    pointer = IEDraster(:, 1) >= timeWinBE(1) & IEDraster(:, 1) <=timeWinBE(2);
    IEDrasterWin = IEDraster(pointer, :);
    IEDrasterWin_Tidx = IEDrasterWin(:, 1);
    IEDrasterWin_Tidx = IEDrasterWin_Tidx - timeWinBE(1) + 1;
    IEDrasterWin_time = timeWin(IEDrasterWin_Tidx);
    IEDrasterWin_Cidx = IEDrasterWin(:, 2);
    IEDrasterWin_value = zeros(1, length(IEDrasterWin_Cidx));
    for n = 1:length(IEDrasterWin_Cidx)
        IEDrasterWin_value(n) = signalsWin(IEDrasterWin_Cidx(n), IEDrasterWin_Tidx(n));
    end
    scatter(AX, IEDrasterWin_time, IEDrasterWin_value, 20, [0, 0, 0], 'filled', 'Marker', 'o');
end