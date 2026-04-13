function [beta, residual, stats, Curve] = fit_lsqcurvefit(x, y, modelfun, beta0, alpha)
%F_LSQCURVEFIT - Nonlinear least squares curve fitting with comprehensive statistics
%   [pHat, ci, residual, stat] = F_LSQCURVEFIT(x, y, modelfun, p0, lb, ub)
%   fits the nonlinear model defined by modelfun to the data (x, y) using
%   constrained nonlinear least squares optimization. The function returns
%   parameter estimates, confidence intervals, residuals, and goodness-of-fit
%   statistics.
%
% INPUTS:
%   x       - Independent variable data (vector or matrix)
%   y       - Dependent variable (response) data (vector)
%   modelfun - Function handle to the nonlinear model: y_hat = modelfun(p, x)
%   p0      - Initial parameter estimates (vector)
%   lb      - Lower bounds for parameters (vector, same size as p0)
%   ub      - Upper bounds for parameters (vector, same size as p0)
%
% OUTPUTS:
%   pHat    - Optimal parameter estimates (vector)
%   ci      - 95% confidence intervals for parameters (k×2 matrix)
%   residual - Residuals: y - modelfun(pHat, x) (vector)
%   stat    - Structure containing comprehensive goodness-of-fit statistics
%
% NOTES:
%   1. The function assumes the model includes an intercept term
%   2. Statistical inference (p-values, confidence intervals) is based on
%      linear approximations and should be interpreted cautiously for
%      highly nonlinear models or small sample sizes
%   3. Requires the Optimization Toolbox and Statistics and Machine Learning Toolbox
%
% See also: lsqcurvefit, nlparci, nlinfit
        x = double(x);
        y = double(y(:));

        %% ===== Degrees of Freedom =====
        n = numel(y);                    % Sample size (number of observations)
        k = numel(beta0);              % Number of estimated parameters
        df_tot = n - 1;                     % Degrees of freedom for total variation
        df_model = k - 1;               % Degrees of freedom for model (assumes intercept)
        df_res = n - k;                     % Degrees of freedom for residuals

        %% Nonlinear least squares fitting
        % Perform constrained nonlinear least squares optimization
        % pHat   - Optimal parameter estimates minimizing sum of squared residuals
        % SSres  - Sum of squared residuals (RSS) at optimal solution
        % residual - Vector of residuals: y - modelfun(pHat, x)
        % exitflag - Convergence indicator (>0: converged, 0: max iterations, <0: failed)
        % output   - Structure with optimization details (iterations, algorithm, etc.)
        % lambda   - Lagrange multipliers for bound constraints
        % J        - Jacobian matrix (n×k) of partial derivatives at solution
        lb = -1*Inf .* ones(size(beta0));
        ub = Inf .* ones(size(beta0));
        options = optimoptions('lsqcurvefit', 'Display', 'off');
        [beta, SSres, residual_neg, exitflag, output, lambda, J] ...
            = lsqcurvefit(modelfun, beta0, x, y, lb, ub, options);
        % Display warning if optimization did not converge properly
        if exitflag <= 0
            warning('Optimization may not have converged. Exit flag: %d', exitflag);
            fprintf('Optimization message: %s\n', output.message);
        end
        % Convert residual sign to y - yhat
        residual = -1 * residual_neg;
        
        % GOODNESS-OF-FIT METRICS
       y_mean = mean(y);
       SStot = sum((y - y_mean).^2);
       [R2, R2adj, RMSE] = fit_Rsquare(SSres, SStot, df_res, df_tot);
       if R2 <= 0
           test = 1;
       end
        
       % AIC / AICc / BIC (SSE-based Gaussian form)
       [AIC, AICc, BIC] = fit_aicbic(SSres, n, k);

       % t-stats, p-values, and (1-alpha) coefficient confidence intervals
        CovB = ( SSres / df_res) * inv(J'*J);
        [beta_tval,  beta_pval, beta_CI] = fit_betaPval(beta, df_res, CovB, alpha);

        % Tests whether the model explains significant variance beyond mean-only.
        [model_Fval, model_pval] = fit_modelFtest(SSres, SStot, df_model, df_res);

        %% ====== OUTPUT STATISTICS STRUCTURE ======
        stats.exitflag = exitflag;
        stats.message  = output.message;
        stats.R2 = R2;
        stats.R2adj = R2adj;
        stats.RMSE = RMSE;
        stats.AIC = AIC;
        stats.BIC = BIC;
        stats.beta_tval = beta_tval;
        stats.beta_pval  = beta_pval;
        stats.beta_CI = beta_CI;
        stats.model_Fval = model_Fval;
        stats.model_pval = model_pval;

        %% ====== FITTED CURVE & CONFIDENCE INTERVALS ======
        % nlpredci returns ypred and delta such that CI_mean = ypred ± delta
        % Provide CovB and MSE to avoid recomputing/ambiguity.
        xrange = [];
        for xn = 1:size(x, 2)
            xrange(:, xn) = linspace(min(x(:, xn)), max(x(:, xn)), 1000);
        end
        try 
            [ypred, delta_mean] = nlpredci(modelfun, xrange, beta, residual, 'Covar', CovB, 'Alpha', alpha);
        catch
            CovB = full(CovB);
            [ypred, delta_mean] = nlpredci(modelfun, xrange, beta, residual, 'Covar', CovB, 'Alpha', alpha);
        end
        Curve.x = xrange;
        Curve.y = ypred;
        Curve.yCI_lower = ypred - delta_mean;   % (1-alpha) CI for mean response
        Curve.yCI_upper = ypred + delta_mean;
end