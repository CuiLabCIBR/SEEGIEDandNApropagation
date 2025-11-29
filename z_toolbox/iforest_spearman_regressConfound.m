function [R, P] = iforest_spearman_regressConfound(x, y, confound)
    % IFOREST_SPEARMAN_REGRESSCONFOUND Computes Spearman correlation after outlier removal and confound adjustment
    %
    % This function performs a comprehensive analysis that:
    % 1. Identifies and removes outliers using Isolation Forest on the combined dataset
    % 2. Adjusts both variables for the effects of a confounding variable using polynomial regression
    % 3. Computes Spearman correlation between the residual values across multiple polynomial orders
    %
    % Inputs:
    %   x       - Primary variable vector (numeric)
    %   y       - Secondary variable vector (numeric, same length as x)
    %   confound - Confounding variable vector (numeric, same length as x and y)
    %             This variable represents a potential source of spurious correlation
    %             between x and y that needs to be statistically controlled
    %
    % Outputs:
    %   R - 1x6 vector of Spearman correlation coefficients for polynomial orders 1-6
    %   P - 1x6 vector of p-values for each correlation coefficient
    %
    % Methodology:
    %   1. Outlier Detection: Uses Isolation Forest algorithm to identify multivariate outliers
    %      in the combined space of x, y, and the confound variable
    %   2. Confound Adjustment: For each polynomial order (1-6), removes the effect of the
    %      confound variable from both x and y using polynomial regression
    %   3. Correlation Analysis: Computes Spearman correlation between the residual values
    %      of x and y (after confound adjustment) to reveal the true relationship
    %
    % Note: This approach helps distinguish genuine relationships from spurious correlations
    %       caused by shared dependence on a confounding variable

    % Apply Isolation Forest algorithm to detect outliers in trivariate space
    % The algorithm identifies anomalies in the combined distribution of x, y, and confound
    [forest, ~, scores] = iforest([x, y, confound], "ContaminationFraction", 0.1);
    
    % Identify inliers: data points with anomaly scores below the threshold
    % These points represent the "normal" data after removing multivariate outliers
    retainIdx = find(scores < forest.ScoreThreshold);
    
    % Extract inlier values for all variables
    x1 = x(retainIdx);          % Primary variable without outliers
    y1 = y(retainIdx);          % Secondary variable without outliers
    confound1 = confound(retainIdx); % Confound variable without outliers
    
    % Initialize output vectors for correlation coefficients and p-values
    R = zeros(1, 6);  % Will store Spearman's rho for polynomial orders 1-6
    P = zeros(1, 6);  % Will store corresponding p-values
    
    % Analyze relationship across multiple polynomial orders (1-6)
    % This allows determining the optimal degree of polynomial needed to model
    % and remove the confound effect from both variables
    for order = 1:6
        % Remove confound effect from x using polynomial regression
        % Fits a polynomial of specified order to model x as a function of confound
        [polyModel, ~, ~] = f_polyfit(confound1, x1, order);
        xresidual = polyModel.residual;  % Residuals represent x after removing confound effect
        
        % Remove confound effect from y using polynomial regression
        % Fits a polynomial of specified order to model y as a function of confound
        [polyModel, ~, ~] = f_polyfit(confound1, y1, order);
        yresidual = polyModel.residual;  % Residuals represent y after removing confound effect
        
        % Calculate Spearman correlation between the confound-adjusted variables
        % This measures the relationship between x and y after removing shared variance
        % with the confound variable at the specified polynomial order
        [R(order), P(order)] = corr(xresidual, yresidual, "type", "Spearman");
    end
end