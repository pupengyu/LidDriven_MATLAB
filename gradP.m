%% This is function to compute the gradient of p
function out = gradP(cellType, p, ii, N)
    global L mu h
%     cellType = getCellType(ii, N);
    out = zeros(2, 1);
    
    %% dp/dx
    if contains(cellType, 'L')
        out(1) = (p(ii + 1) - p(ii)) / (2 * h); % pb = pC
    elseif contains(cellType, 'R')
        out(1) = (p(ii) - p(ii - 1)) / (2 * h);
    else
        out(1) = (p(ii + 1) - p(ii - 1)) / (2 * h);
    end
    %% dp/dy
    if contains(cellType, 'D')
        out(2) = (p(ii + N) - p(ii)) / (2 * h);
    elseif contains(cellType, 'U')
        out(2) = (p(ii) - p(ii - N)) / (2 * h);
    else
        out(2) = (p(ii + N) - p(ii - N)) / (2 * h);
    end