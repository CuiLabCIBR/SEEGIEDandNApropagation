function [R2, R2adj, RMSE] = fit_Rsquare(SSres, SStot, df_res, df_tot)
%F_RSQUARE Compute R^2, adjusted R^2, and RMSE from sums of squares and DOF.
%
%   [R2, R2adj, RMSE] = f_Rsquare(SSres, SStot, df_res, df_tot)
%
%   This helper computes three common regression fit metrics given:
%     - SSres : Residual sum of squares (a.k.a. SSE) = sum((y - yHat).^2)
%     - SStot : Total sum of squares (a.k.a. SST)       = sum((y - mean(y)).^2)
%     - df_res: Residual degrees of freedom          = n - k
%     - df_tot: Total degrees of freedom             = n - 1
%
%   OUTPUTS
%   -------
%   R2    : Coefficient of determination. Proportion of variance in y explained
%           by the model relative to an intercept-only baseline.
%           R2 = 1 - SSres / SStot
%
%   R2adj : Adjusted R^2. Complexity-penalized R^2 that accounts for degrees of freedom:
%           R2adj = 1 - (SSres/df_res) / (SStot/df_tot)
%           This reduces the tendency of R2 to increase as you add parameters.
%
%   RMSE  : Root mean squared error. Typical magnitude of residuals in y-units:
%           RMSE = sqrt(SSres / df_res)
%
%   NOTES / EDGE CASES
%   ------------------
%   - If SStot == 0 (i.e., y is constant), R2 is undefined; we return NaN.
%   - If df_res <= 0 or df_tot <= 0, adjusted R2 and/or RMSE are undefined.
%   - SSres and SStot should be nonnegative; negative values indicate upstream
%     numerical/logic issues.

     % Coefficient of determination, proportion of variance explained
    if SStot > 0
        R2 = 1 - SSres / SStot;       % R-squared
    else
        R2 = NaN;
    end
    % Adjusted R-squared: penalizes model complexity
    if df_tot > 0 && df_res > 0 && SStot > 0
        R2adj = 1 - (SSres/df_res) / (SStot/df_tot);
    else
        R2adj = NaN;
    end

    % Mean squared error (MSE) and root MSE (RMSE)
    if df_res > 0
        RMSE = sqrt(SSres / df_res);
    else
        RMSE = NaN;
    end
end