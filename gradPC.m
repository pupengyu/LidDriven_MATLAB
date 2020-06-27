%% This is function to compute the gradient of pressure correction
function out = gradPC(pC, ii, N)
    global L mu h
    cellType = getCellType(ii, N);
    out = zeros(2, 1);
    
    %% dpC/dx
    if contains(cellType, 'L')
        out(1) = (pC(ii + 1) + pC(ii)) / (2 * h); % pC_b = 0
    elseif contains(cellType, 'R')
        out(1) = -(pC(ii) + pC(ii - 1)) / (2 * h);
    else
        out(1) = (pC(ii + 1) - pC(ii - 1)) / (2 * h);
    end
    %% dpC/dy
    if contains(cellType, 'D')
        out(2) = (pC(ii + N) + pC(ii)) / (2 * h);
    elseif contains(cellType, 'U')
        out(2) = -(pC(ii) + pC(ii - N)) / (2 * h);
    else
        out(2) = (pC(ii + N) - pC(ii - N)) / (2 * h);
    end