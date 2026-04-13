function [model_Fval, model_pval] = fit_modelFtest(SSres, SStot, df_model, df_res)
%F_MODEL_FTEST Overall (global) F-test for regression model significance.
%
%   [model_Fval, model_pval] = f_model_Ftest(SSres, SStot, df_model, df_res)
%
%   This function performs the classical overall F-test that evaluates whether
%   a regression model explains a significant amount of variance in y compared
%   to an intercept-only (mean-only) model.
%
%   In standard linear regression:
%       SStot = sum((y - mean(y)).^2)          % total variability
%       SSres = sum((y - yHat).^2)             % residual (unexplained) variability
%       SSreg = SStot - SSres                  % explained variability by the model
%
%   The F statistic compares mean squares:
%       MS_model = (SStot - SSres) / df_model
%       MS_error = SSres / df_res
%       F        = MS_model / MS_error
%
%   Under H0 (the model does not improve over intercept-only; equivalently all
%   slope-like coefficients are zero), F follows an F distribution with
%   (df_model, df_res) degrees of freedom.
%
%   INPUTS
%   ------
%   SSres    : Residual sum of squares (SSE), must be >= 0
%   SStot    : Total sum of squares (SST), must be >= 0
%   df_model : Model degrees of freedom (typically k - 1, excluding intercept)
%   df_res   : Residual degrees of freedom (typically n - k)
%
%   OUTPUTS
%   -------
%   model_Fval : F statistic value (scalar)
%   model_pval : Right-tail p-value P(F >= model_Fval | H0)
%
%   NOTES / EDGE CASES
%   ------------------
%   - If df_model <= 0 or df_res <= 0, the test is undefined -> returns NaN.
%   - If SSres == 0, then MS_error = 0 and F tends to Inf (perfect fit).
%     Here we compute directly; fcdf(Inf,...) -> 1, so p -> 0.
%   - If SStot < SSres due to numerical issues, (SStot - SSres) becomes negative.
%     In that case the statistic is not meaningful; consider guarding upstream.

    % Validate degrees of freedom and nonnegativity of SSE
    if (df_model > 0) && (df_res > 0) && (SSres >= 0)
        % Mean square explained by the model relative to intercept-only baseline
        MS_model = (SStot - SSres) / df_model;

        % Mean square error (residual mean square)
        MS_error = SSres / df_res;

        % F statistic: ratio of explained mean square to residual mean square
        model_Fval = MS_model / MS_error;

        % Right-tail p-value from the F distribution
        model_pval = 1 - fcdf(model_Fval, df_model, df_res);
    else
        model_Fval = NaN;
        model_pval = NaN;
    end
end