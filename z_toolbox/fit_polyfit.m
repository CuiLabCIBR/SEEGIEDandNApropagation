function [beta, residual, stats, Curve] = fit_polyfit(x, y, order, alpha)
%F_POLYFIT Polynomial regression with confidence intervals and diagnostics.
%
%   [beta, residual, stat, Curve] = f_polyfit(x, y, order)
%
%   This function fits a polynomial of a specified order to (x, y) data using
%   ordinary least squares (via MATLAB POLYFIT). It returns coefficient estimates,
%   residuals, common goodness-of-fit metrics, coefficient-level inferential stats,
%   and fitted-curve confidence/prediction intervals on a dense x-grid.
%
%   INPUTS
%   ------
%   x     : Independent variable (vector). Will be reshaped to a column vector.
%   y     : Dependent variable (vector). Will be reshaped to a column vector.
%   order : Polynomial order (non-negative integer). Number of parameters k = order+1.
%
%   OUTPUTS
%   -------
%   beta     : Polynomial coefficients (descending powers; same as POLYFIT output).
%   residual : Residuals y - yHat at the original x locations.
%   stat     : Structure with diagnostics and inference results.
%   Curve    : Structure with fitted curve and interval bands on a dense grid.
    %% ====== POLYNOMIAL FITTING (ORDINARY LEAST SQUARES) ======
    % POLYFIT solves the least-squares problem using a QR factorization of the Vandermonde design matrix. 
    % S includes QR-related information; 
    % mu is used for centering/scaling x to improve numerical stability.
    [beta, S, mu] = polyfit(x, y, order);

   % Calculate predicted values and residuals
    yHat = polyval(beta, x, S, mu);
    residual = y - yHat;

   %% ====== SUMS OF SQUARES & DOF ======
    % Total sum of squares (SStot): variance around the mean
    % Residual sum of squares (SSres): unexplained variance after fitting
    SStot = sum((y - mean(y)).^2);       % Total sum of squares
    SSres = sum(residual.^2);              % Residual sum of squares
    
    % Parameter count and degrees of freedom
    n = numel(y);
    k = order + 1;                                 % Number of parameters (order + intercept)
    df_tot = n - 1;                                 % Degrees of freedom for total variation
    df_model = k - 1;                           % Degrees of freedom for model
    df_res = n - k;                                % Degrees of freedom for residuals

    %% ====== COVARIANCE OF COEFFICIENTS (CovB) ======
    % residual variance estimate 
    if df_res > 0
        sigma2 = SSres / df_res; % kept for completeness; not used further
    else
        sigma2 = NaN; 
    end
    % For linear regression, Cov(beta) = sigma^2 * inv(X'X).
    % POLYFIT provides the upper-triangular R from QR such that X = Q*R.
    % Then X'X = R'*R and inv(X'X) = inv(R)*inv(R)'.
    % Use triangular solves instead of inv() for numerical stability.
    invR = inv(S.R);                    % stable equivalent of inv(R)
    CovB = sigma2 * (invR * invR');  % covariance matrix of coefficients
    
   %% ====== GOODNESS-OF-FIT METRICS ======
   % R2 / adjusted R2 / RMSE
   [R2, R2adj, RMSE] = fit_Rsquare(SSres, SStot, df_res, df_tot);
    
   % AIC / AICc / BIC (SSE-based Gaussian form)
   [AIC, AICc, BIC] = fit_aicbic(SSres, n, k);

   %% ====== PARAMETER-LEVEL INFERENCE ======
   % t-stats, p-values, and (1-alpha) coefficient confidence intervals
    [beta_tval, beta_pval, beta_CI] = fit_betaPval(beta, df_res, CovB, alpha);
    
    %% ====== OVERALL MODEL F-TEST ======
    % Tests whether the model explains significant variance beyond mean-only.
    [model_Fval, model_pval] = fit_modelFtest(SSres, SStot, df_model, df_res);

    %% ====== OUTPUT STATISTICS STRUCTURE ======
    stats.R2 = R2;
    stats.R2adj  = R2adj;
    stats.RMSE = RMSE;
    stats.AIC = AIC;
    stats.BIC = BIC;
    stats.beta_tval = beta_tval;
    stats.beta_pval = beta_pval;
    stats.beta_CI = beta_CI;
    stats.model_Fval = model_Fval;
    stats.model_pval = model_pval;

    %% ====== FITTED CURVE & CONFIDENCE INTERVALS ======
    % Generate fine grid for smooth curve
    xrange = linspace(min(x), max(x), 1000);
    % Calculate fitted values and prediction errors on the grid
    [ypred, delta_fit] = polyval(beta, xrange, S, mu);
    % Store curve data
    Curve.x = xrange; 
    Curve.y = ypred;
    % 95% Prediction confidence intervals for the fitted curve
    Curve.yCI_lower = ypred - 2* delta_fit;
    Curve.yCI_upper = ypred + 2 * delta_fit;
end