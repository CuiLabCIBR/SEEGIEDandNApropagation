function [AIC, AICc, BIC] = fit_aicbic(SSres, n, k)
%F_AICBIC Compute AIC, AICc, and BIC from residual sum of squares (SSE).
%
%   [AIC, AICc, BIC] = f_aicbic(SSres, n, k)
%
%   This helper computes information criteria commonly used for model
%   comparison/selection under the assumption of Gaussian i.i.d. errors.
%   For least-squares regression, the (negative) log-likelihood is proportional
%   to n*log(SSE/n), so AIC/BIC can be computed directly from SSres (SSE).
%
%   INPUTS
%   ------
%   SSres  : Residual sum of squares (SSE) = sum((y - yHat).^2)
%   n      : Number of observations
%   k      : Number of free parameters in the model (including intercept)
%
%   OUTPUTS
%   -------
%   AIC  : Akaike Information Criterion
%          AIC = n*log(SSE/n) + 2k
%
%   AICc : Small-sample corrected AIC (recommended when n is not large vs k)
%          AICc = AIC + (2k(k+1)) / (n - k - 1)
%          Defined only when (n - k - 1) > 0.
%
%   BIC  : Bayesian Information Criterion (a.k.a. Schwarz criterion)
%          BIC = n*log(SSE/n) + k*log(n)
%
%   INTERPRETATION
%   --------------
%   - Lower AIC/AICc/BIC indicates a preferred model (within the same dataset).
%   - AIC tends to favor predictive accuracy; BIC penalizes complexity more strongly.
%   - These criteria are meaningful when comparing models fit to the same response
%     variable and the same observations.
%
%   NOTES / EDGE CASES
%   ------------------
%   - If SSres <= 0 (perfect fit or numerical issue), log(SSres/n) is invalid;
%     the function returns NaN for all criteria.
%   - AICc is undefined when n <= k + 1 (division by zero or negative).
%   - The constant terms in the Gaussian log-likelihood cancel in comparisons,
%     so the SSE-based forms are standard for model selection.

    %% --- Information criteria ---
    % Use the SSE-based Gaussian form: n*log(SSE/n) + penalty(k).
    if (n > 0) && (SSres > 0)
        AIC = n * log(SSres / n) + 2 * k;
        BIC = n * log(SSres / n) + k * log(n);

        % Small-sample correction (AICc)
        if (n - k - 1) > 0
            AICc = AIC + (2 * k * (k + 1)) / (n - k - 1);
        else
            AICc = NaN;
        end
    else
        % Invalid inputs (e.g., SSres <= 0 or n <= 0)
        AIC = NaN;
        AICc = NaN;
        BIC = NaN;
    end
end