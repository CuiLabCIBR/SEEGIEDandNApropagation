function order_opt = f_polyfit_crossValidation(x, y, orderList, holdOut, rpN)

%%
    BICmean = zeros(1, length(orderList));
    for m = 1:length(orderList)
        order = orderList(m);
        BIC = zeros(1, rpN);
        for rp = 1:rpN
            cv = cvpartition(length(x), "HoldOut", holdOut);
            x_train = x(training(cv));
            y_train = y(training(cv));
            x_test = x(test(cv));
            y_test = y(test(cv));
            [p_train, s, mu] = polyfit(x_train, y_train, order);
            y_test_pred = polyval(p_train, x_test, s, mu);
            test_SSE = sum((y_test - y_test_pred).^2);
            % BIC
            n = length(y_test);      
            k = length(p_train);
            BIC(rp) = n*log(test_SSE/n) + k*log(n);
        end
        BICmean(m) = mean(BIC);
    end
    [~, index] = min(BICmean);
    order_opt = orderList(index);
end

