
function fit_gam(dataTable, formula, alpha)
    
    rng("default");
    CVMdl  = fitrgam(dataTable, formula, 'CrossVal', 'on', 'KFold', 2);
    yHat = kfoldPredict(CVMdl);
    L = kfoldLoss(CVMdl);
    
end