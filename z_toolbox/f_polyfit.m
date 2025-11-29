function [polyModel, Curve, Ftest] = f_polyfit(x, y, order)
%
%%
    SST = sum((y - mean(y)).^2);% Total sum of squares

    %
    [polyP, polyS, polyMu] = polyfit(x, y, order);
    
    % Adjusted Rsquare
    y_pred = polyval(polyP, x, polyS, polyMu);
    residual = y - y_pred;
    SSE = sum((residual).^2); % Sum of squares error
    SSR = SST - SSE;            % The sum of squares due to regression
    n = length(y);
    k = length(polyP) - 1;
    adjR2 = 1 - (SSE/(n-k-1)) / (SST/(n-1));
    polyModel.p = polyP;
    polyModel.S = polyS;
    polyModel.mu = polyMu;
    polyModel.residual = residual;
    polyModel.SSE = SSE;
    polyModel.SSR = SSR;
    polyModel.SST = SST;
    polyModel.adjR2 = adjR2;

    % The curve of fit
    xfit = linspace(min(x), max(x), 1000);
    yfit = polyval(polyP, xfit, polyS, polyMu);
    Curve.x = xfit;
    Curve.y = yfit;

    % F test
    df_model = k;              
    df_error = n - k - 1;       
    Fstat = (SSR/df_model) / (SSE/df_error);
    p_value = 1 - fcdf(Fstat, df_model, df_error);
    Ftest.Fstat = Fstat;
    Ftest.df_model = df_model;
    Ftest.df_error = df_error;
    Ftest.p_value = p_value;
end