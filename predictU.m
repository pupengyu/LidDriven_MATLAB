% This is a function to predict U %
% by solving momentum equations %

function [UNEW, D] = predictU(cellType, U, p, N)
    global L mu h ub alphaU alphaP
    
    %% Define coefficient matrices A and right-hand side vectors B
    A = cell(2, 1);
    A{1} = zeros(N * N, N * N);
    A{2} = zeros(N * N, N * N);
    
    B = zeros(N * N, 2);
    
    %% Neighbor surfaces
    S_nb = [-h, h, 0, 0; 0, 0, -h, h];  % Surface vector for w, e, s, n
    N_nb = length(S_nb(1, :));
    
    %% Temporary variables
    tempGradU = cell(N * N, 1);
    tempGradP = cell(N * N, 1);
    tempM = zeros(N * N, 4);
    nb_id = zeros(N * N, 4);
    
    for ii = 1:N*N
        tempGradU{ii} = gradU(cellType{ii}, U, ii, N);
        tempGradP{ii} = gradP(cellType{ii}, p, ii, N);
        nb_id(ii, :) = [ii - 1, ii + 1, ii - N, ii + N];
        
        %% Boundary contributions
        if contains(cellType{ii}, 'L')
            nb_id(ii, 1) = 0;       
        end
        if contains(cellType{ii}, 'R')
            nb_id(ii, 2) = 0;
        end
        if contains(cellType{ii}, 'D')
            nb_id(ii, 3) = 0;
        end
        if contains(cellType{ii}, 'U')
            nb_id(ii, 4) = 0;
        end
        for jj = 1:4
            if nb_id(ii, jj) > 0
                tempM(ii, jj) = getM(jj, S_nb, U, ii, N);
            end
        end
    end
        
    for ii = 1:N*N

        mf = zeros(N_nb, 1);
        aF = zeros(N_nb, 1);
        aC = zeros(1, 2);
        bC = zeros(1, 2);

        
        %% Boundary contributions
        if contains(cellType{ii}, 'L')

            aC(2) = aC(2) + 2 * mu;
        end
        if contains(cellType{ii}, 'R')

            aC(2) = aC(2) + 2 * mu;
        end
        if contains(cellType{ii}, 'D')

            aC(1) = aC(1) + 2 * mu;
        end
        if contains(cellType{ii}, 'U')

            aC(1) = aC(1) + 2 * mu;
            bC(1) = bC(1) + 2 * mu * ub;
        end
        
        
        %% Interior contributions
        for jj = 1:N_nb
            if nb_id(ii, jj) > 0    %% if interior faces
                aF(jj) = -max(-tempM(ii, jj), 0) - mu;
                aC = aC + max(tempM(ii, jj), 0) + mu;
                bC = bC + (mu * 0.5 * (tempGradU{ii} + tempGradU{nb_id(ii, jj)}) * S_nb(:, jj))';   %% explicit diffusion term
            else
                aF(jj) = 0;
            end
        end
        bC = bC - tempGradP{ii}' * h ^ 2;
        

        
        %% Add to matrices
        % diagnal
        A{1}(ii, ii) = aC(1) / alphaU;
        A{2}(ii, ii) = aC(2) / alphaU;
        % neighbors
        for jj = 1:N_nb
            if nb_id(ii, jj) > 0
                A{1}(ii, nb_id(ii, jj)) = aF(jj);
                A{2}(ii, nb_id(ii, jj)) = aF(jj);
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