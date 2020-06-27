%% This is function to compute the gradient of U
function out = gradU(U, ii, N)
    global L mu h ub
    cellType = getCellType(ii, N);
    gradu = zeros(2, 1);
    gradv = zeros(2, 1);
    
    %% dU/dx
    if contains(cellType, 'L')
        gradu(1) = (U(ii + 1, 1) + U(ii, 1)) / (2 * h); % u = uw should manifest
        gradv(1) = (U(ii + 1, 2) + U(ii, 2)) / (2 * h);
    elseif contains(cellType, 'R')
        gradu(1) = -(U(ii, 1) + U(ii - 1, 1)) / (2 * h);
        gradv(1) = -(U(ii, 2) + U(ii - 1, 2)) / (2 * h);
    else
        gradu(1) = (U(ii + 1, 1) - U(ii - 1, 1)) / (2 * h);
        gradv(1) = (U(ii + 1, 2) - U(ii - 1, 2)) / (2 * h);
    end
    %% dU/dy
    if contains(cellType, 'D')
        gradu(2) = (U(ii + N, 1) + U(ii, 1)) / (2 * h);
        gradv(2) = (U(ii + N, 2) + U(ii, 2)) / (2 * h);
    elseif contains(cellType, 'U')
        gradu(2) = (2 * ub - (U(ii, 1) + U(ii - N, 1))) / (2 * h);
        gradu(2) = - (U(ii, 1) + U(ii - N, 1)) / (2 * h);
    else
        gradu(2) = (U(ii + N, 1) - U(ii - N, 1)) / (2 * h);
        gradv(2) = (U(ii + N, 1) - U(ii - N, 1)) / (2 * h);
    end
    
    out = [gradu, gradv];