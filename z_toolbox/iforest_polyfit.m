function adjR2 = iforest_polyfit(x, y)
    % IFOREST_POLYFIT Computes adjusted R-squared values for polynomial fits after outlier removal
    %
    % This function performs the following analysis:
    % 1. Identifies and removes outliers using Isolation Forest algorithm
    % 2. Fits polynomial models of varying orders to the cleaned data
    % 3. Returns adjusted R-squared values for each polynomial order
    %
    % Inputs:
    %   x - Independent variable vector (predictor)
    %   y - Dependent variable vector (response)
    %
    % Output:
    %   adjR2 - 1x6 vector of adjusted R-squared values for polynomial orders 1 through 6
    %
    % Methodology:
    %   1. Outlier Detection: Applies Isolation Forest to identify anomalies in the
    %      bivariate distribution of x and y
    %   2. Data Cleaning: Retains only inliers (non-outlier points) for analysis
    %   3. Polynomial Modeling: Fits polynomial regression models of increasing
    %      complexity (orders 1-6) to the cleaned data
    %   4. Model Evaluation: Computes adjusted R-squared for each model, which
    %      penalizes model complexity to prevent overfitting
    %
    % Note: Adjusted R-squared is preferred over regular R-squared as it accounts
    %       for the number of predictors in the model, providing a more accurate
    %       measure of model fit quality, especially when comparing models with
    %       different numbers of parameters

    % Apply Isolation Forest algorithm to detect outliers in bivariate space
    % The algorithm identifies anomalies in the joint distribution of x and y
    [forest, ~, scores] = iforest([x, y], "ContaminationFraction", 0.1);
    
    % Identify inliers: data points with anomaly scores below the threshold
    % These points represent the "normal" data after removing multivariate outliers
    retainIdx = find(scores < forest.ScoreThreshold);
    
    % Extract inlier values for both variables
    x1 = x(retainIdx); % Independent variable without outliers
    y1 = y(retainIdx); % Dependent variable without outliers
    
    % Initialize output vector for adjusted R-squared values
    adjR2 = zeros(1, 6);
    
    % Fit polynomial models of increasing complexity (orders 1-6)
    for order = 1:6
        % Fit polynomial regression model of specified order
        % Models the relationship: y = f(x) where f is a polynomial of degree 'order'
        [polyModel, ~, ~] = f_polyfit(x1, y1, order);
        
        % Store the adjusted R-squared value for this model
        % Adjusted R-squared accounts for model complexity by penalizing additional terms
        adjR2(order) = polyModel.adjR2;
    end
end
