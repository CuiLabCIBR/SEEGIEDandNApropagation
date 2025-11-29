function cliffDelta = f_cliff_delta(x, y)
    n1 = length(x); n2 = length(y);
    for ny = 1:length(y)
        Nwin(ny) = sum(x>y(ny));
        Nloss(ny) = sum(x<y(ny));
        Ntie(ny) = sum(x==y(ny));
    end
    Nwin = sum(Nwin); Nloss = sum(Nloss); Ntie = sum(Ntie);
    cliffDelta = (Nwin-Nloss)/(n1*n2-Ntie);
end