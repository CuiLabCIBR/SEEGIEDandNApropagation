function [R2, trainR2, testR2, y_pred_cv] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold)
%FIT_LSQCURVEFIT_CV Perform k-fold cross-validation for nonlinear curve fitting.
%
%   [R2, trainR2, testR2, y_pred_cv] = fit_lsqcurvefit_cv(x, y, modelfun, beta0, kfold)
%
%   Inputs
%   ------
%   x        : Predictor matrix of size N-by-P, where N is the number of
%              samples and P is the number of predictors.
%   y        : Response vector of length N.
%   modelfun : Function handle for the model. The function should have the form:
%                  yhat = modelfun(beta, x)
%              where beta is the parameter vector and x is the predictor matrix.
%   beta0    : Initial parameter vector for nonlinear least-squares fitting.
%   kfold    : Number of folds for cross-validation.
%
%   Outputs
%   -------
%   R2       : Overall cross-validated R^2, computed from all out-of-fold predictions.
%   trainR2  : Training-set R^2 for each fold (kfold-by-1 vector).
%   testR2   : Test-set R^2 for each fold (kfold-by-1 vector).
%   y_pred_cv: Cross-validated predicted values for all samples (N-by-1 vector).
%%    
    % Number of samples
    n = size(x, 1);

    % Randomly shuffle sample indices
    idxRandom = randperm(n);

    % Define fold boundaries
    foldStart = floor(linspace(1, n, kfold+1));
    foldStart = foldStart(1:end-1);
    foldEnd = foldStart-1;
    foldEnd = foldEnd(2:end);
    foldEnd = [foldEnd, n];
    
    % Cross-validation loop
    alpha = 0.05;
    y_pred_cv = nan(n, 1);
    y_cv = nan(n, 1);
    for k = 1:kfold
        % Create a logical indicator for the current test fold
        test01 = zeros(n, 1);
        test01(foldStart(k):foldEnd(k)) = 1;

        % Map fold positions back to the shuffled sample indices
        testIdx = idxRandom(test01==1);
        trainIdx = idxRandom(test01==0);

        % Extract training and test data
        xTest = x(testIdx, :);
        yTest = y(testIdx);
        xTrain= x(trainIdx, :);
        yTrain = y(trainIdx);

        % Fit the model on the training set
        [beta, ~, stats, ~] = fit_lsqcurvefit(xTrain, yTrain, modelfun, beta0, alpha);
        % Store training R^2 returned by the fitting function
        trainR2(k) = stats.R2;
        % Predict the test set
        yTestPred = modelfun(beta, xTest);

        % Compute test-set R^2
        testSStot = sum((yTest - mean(yTest)).^2);
        testResidual = yTest - yTestPred;
        testSSres = sum(testResidual.^2);
        testR2(k) = 1- testSSres/testSStot;
        
        % Save out-of-fold predictions
        y_pred_cv(testIdx) = yTestPred;
        y_cv(testIdx) = yTest;
    end
    
    % Compute overall cross-validated R^2 from all out-of-fold predictions
    residual = y_cv - y_pred_cv;
    SSres  = sum(residual.^2);
    SStot = sum((y_cv - mean(y_cv)).^2);
    R2 = 1 - SSres / SStot;
end