% This is a function to predict U %
% by solving momentum equations %

function [UNEW, D] = predictU(U, p, N)
    global L mu h ub alphaU alphaP
    
    %% Define coefficient matrices A and right-hand side vectors B
    A = cell(2, 1);
    A{1} = zeros(N * N, N * N);
    A{2} = zeros(N * N, N * N);
    
    B = zeros(N * N, 2);
    
    %% Neighbor surfaces
    S_nb = [-h, h, 0, 0; 0, 0, -h, h];  % Surface vector for w, e, s, n
    N_nb = length(S_nb(1, :));
        
    for ii = 1:N*N
        cellType = getCellType(ii, N);
        mf = zeros(N_nb, 1);
        aF = zeros(N_nb, 1);
        aC = zeros(1, 2);
        bC = zeros(1, 2);
        nb_id = [ii - 1, ii + 1, ii - N, ii + N];
        
        %% Boundary contributions
        if contains(cellType, 'L')
            nb_id(1) = 0;
            aC(2) = aC(2) + 2 * mu;
        end
        if contains(cellType, 'R')
            nb_id(2) = 0;
            aC(2) = aC(2) + 2 * mu;
        end
        if contains(cellType, 'D')
            nb_id(3) = 0;
            aC(1) = aC(1) + 2 * mu;
        end
        if contains(cellType, 'U')
            nb_id(4) = 0;
            aC(1) = aC(1) + 2 * mu;
            bC(1) = bC(1) + 2 * mu * ub;
        end
        
%         disp(['cellType : ', cellType]);
        
        %% Interior contributions
        for jj = 1:N_nb
            if nb_id(jj) > 0    %% if interior faces
                mf(jj) = getM(jj, S_nb, U, ii, N);
                aF(jj) = -max(-mf(jj), 0) - mu;
                aC = aC + max(mf(jj), 0) + mu;
                bC = bC + (mu * 0.5 * (gradU(U, ii, N) + gradU(U, nb_id(jj), N)) * S_nb(:, jj))';   %% explicit diffusion term
            else
                mf(jj) = 0;
                aF(jj) = 0;
            end
        end
        bC = bC - gradP(p, ii, N)' * h ^ 2;
        

        
        %% Add to matrices
        % diagnal
        A{1}(ii, ii) = aC(1) / alphaU;
        A{2}(ii, ii) = aC(2) / alphaU;
        % neighbors
        for jj = 1:N_nb
            if nb_id(jj) > 0
                A{1}(ii, nb_id(jj)) = aF(jj);
                A{2}(ii, nb_id(jj)) = aF(jj);
            end
        end
        % right-hand side
        B(ii, :) = bC + (1 - alphaU) / alphaU * aC .* U(ii, :);
    end
    
%% solve momentum equations
    A{1} = sparse(A{1});
    A{2} = sparse(A{2});

    UNEW(:, 1) = A{1} \ B(:, 1); % solve u
    UNEW(:, 2) = A{2} \ B(:, 2); % solve v

    
%% store D for pressure correction
    D(:, 1) = h ^ 2 ./ diag(A{1});
    D(:, 2) = h ^ 2 ./ diag(A{2});
end