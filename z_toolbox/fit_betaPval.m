function [beta_tval, beta_pval,  beta_CI] = fit_betaPval(beta, df_res, CovB, alpha)
%F_BETAPVAL Compute t-stats, two-sided p-values, and confidence intervals for coefficients.
%
%   [beta_pval, beta_tval, beta_CI] = f_betaPval(beta, df_res, CovB, alpha)
%
%   This function performs standard coefficient-level inference for an
%   (approximately) linear regression model:
%     - Standard error:  SE_i = sqrt(Var(beta_i)) = sqrt(diag(CovB))
%     - t-statistic:     t_i  = beta_i / SE_i
%     - two-sided p-val: p_i  = 2 * (1 - tcdf(|t_i|, df_res))
%     - (1-alpha) CI:    beta_i ± t_crit * SE_i, where t_crit = tinv(1-alpha/2, df_res)
%
%   INPUTS
%   ------
%   beta   : Vector of coefficient estimates (any shape; will be vectorized).
%   df_res : Residual degrees of freedom (typically n - k). Must be > 0 to do inference.
%   CovB   : Covariance matrix of beta (k-by-k). Its diagonal entries are Var(beta_i).
%   alpha  : Significance level for confidence intervals (e.g., 0.05 -> 95% CI).
%
%   OUTPUTS
%   -------
%   beta_pval : Two-sided p-values testing H0: beta_i = 0 (k-by-1).
%   beta_tval : t-statistics for each coefficient (k-by-1).
%   beta_CI   : Confidence intervals (k-by-2), columns are [lower, upper].
%
%   NOTES
%   -----
%   - If df_res <= 0, the variance estimate and t-distribution inference are not defined;
%     outputs are set to NaN.
%   - This assumes the usual regression conditions for t-inference (e.g., approximately
%     Gaussian errors / large-sample normality of beta).
%   - CovB must match beta in size: size(CovB,1) should equal numel(beta).

    % (Optional) Default alpha if missing/empty
    if nargin < 4 || isempty(alpha)
        alpha = 0.05;
    end

    % Vectorize beta for consistent output shapes
    beta_col = beta(:);
    k = numel(beta_col);

    % Pre-allocate outputs to safe defaults
    beta_tval = NaN(k, 1);
    beta_pval = NaN(k, 1);
    beta_CI   = NaN(k, 2);

    % Only compute inference when residual degrees of freedom are valid
    if df_res > 0
        % Basic dimension check to prevent silent mismatch errors
        if size(CovB,1) ~= k || size(CovB,2) ~= k
            error('Dimension mismatch: numel(beta)=%d but CovB is %dx%d.', ...
                  k, size(CovB,1), size(CovB,2));
        end

        % Standard errors of coefficients (square root of diagonal variances)
        SE = sqrt(diag(CovB));

        % t-statistics for H0: beta_i = 0
        beta_tval = beta_col ./ SE;

        % Two-sided p-values using Student''s t distribution with df_res degrees of freedom
        beta_pval = 2 * (1 - tcdf(abs(beta_tval), df_res));

        % Critical t value for (1-alpha) confidence interval
        t_crit = tinv(1 - alpha/2, df_res);

        % Confidence interval bounds: beta ± t_crit * SE
        beta_CI_lower = beta_col - t_crit .* SE;
        beta_CI_upper = beta_col + t_crit .* SE;

        % Pack into k-by-2 matrix: [lower, upper]
        beta_CI = [beta_CI_lower, beta_CI_upper];
    end
end