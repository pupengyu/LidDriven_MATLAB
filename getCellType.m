%% Get cell type for boundary condition treatments
function out = getCellType(ii, N)
    %% Distinguish interior cells and boundary cells
    out = '';
    if ii >= 1 && ii <= N
        out = [out, 'D'];
    end
    if ii >= (N - 1) * N + 1 && ii <= N * N
        out = [out, 'U'];
    end
    if mod(ii, N) == 1
        out = [out, 'L'];
    end
    if mod(ii, N) == 0
        out = [out, 'R'];
    end
end