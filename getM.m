%% Get mass flow rate
function out = getM(jj, S_nb, U, ii, N)
    global L mu h
    out = zeros(4, 1);
    UC = U(ii, :);
    switch jj
        case 1
            UF = U(ii - 1, :);
        case 2
            UF = U(ii + 1, :);
        case 3
            UF = U(ii - N, :);
        case 4
            UF = U(ii + N, :);
    end
    out = dot(0.5 * (UF + UC), S_nb(:, jj));
end