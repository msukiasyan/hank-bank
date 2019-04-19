function res = ind2subind(size, ind, subind)
    [i, j, k]       = ind2sub(size, ind);
    switch subind
        case 1
            res = i;
        case 2
            res = j;
        case 3
            res = k;
    end
end