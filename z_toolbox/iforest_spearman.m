function [R, P] = iforest_spearman(x, y)
    % IFOREST_SPEARMAN Computes Spearman correlation after removing outliers using Isolation Forest
    %
    % This function identifies and removes outliers from bivariate data using
    % Isolation Forest algorithm, then calculates Spearman rank correlation
    % between the two variables on the cleaned dataset.
    %
    % Inputs:
    %   x - First variable vector (numeric)
    %   y - Second variable vector (numeric, same length as x)
    %
    % Outputs:
    %   R - Spearman correlation coefficient between x and y after outlier removal
    %   P - p-value indicating statistical significance of the correlation
    %
    % Method:
    %   1. Combines x and y into a bivariate dataset
    %   2. Applies Isolation Forest algorithm with 10% contamination parameter
    %   3. Identifies inliers based on anomaly score threshold
    %   4. Computes Spearman correlation on the cleaned dataset
    
    % Apply Isolation Forest algorithm to detect outliers in bivariate space
    [forest, ~, scores] = iforest([x, y], "ContaminationFraction", 0.1);
    
    % Identify inliers: points with anomaly scores below the threshold
    retainIdx = find(scores < forest.ScoreThreshold);
    
    % Extract inlier values from both variables
    x1 = x(retainIdx);
    y1 = y(retainIdx);
    
    % Calculate Spearman rank correlation coefficient and p-value
    [R, P] = corr(x1, y1, "type", "Spearman");
end