% This is a function to solve pressure
% correction equation and to correct
% the pressure field and velocity field

function [pNew, UNew] = correctP(cellType, U, p, D, N)
global L mu h ub alphaU alphaP

    %% Matrices and vectors
    A = zeros(N * N, N * N);
    B = zeros(N * N, 1);
    UNew = U;
    
    %% Neighbor surfaces
    S_nb = [-h, h, 0, 0; 0, 0, -h, h];  % Surface vector for w, e, s, n
    N_nb = length(S_nb(1, :));
    
    %% Temporary variables
    tempM = zeros(N * N, 4);
    tempGradP = cell(N * N, 1);
    nb_id = zeros(N * N, 4);
    
    for ii = 1:N*N
        tempGradP{ii} = gradP(cellType{ii}, p, ii, N);
        
        nb_id(ii, :) = [ii - 1, ii + 1, ii - N, ii + N];
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
        
        
    end
    
    for ii = 1:N*N
        for jj = 1:4
            if nb_id(ii, jj) > 0
                Df = 0.5 * (D(ii, :) + D(nb_id(ii, jj), :));
                tempM(ii, jj) = getM(jj, S_nb, U, ii, N) - (2 - alphaU) * ((p(nb_id(ii, jj)) - p(ii)) * Df(ceil(jj / 2)) - 0.5 * (diag(Df) * (tempGradP{ii} + tempGradP{nb_id(ii, jj)}))' * S_nb(:, jj));
            end
        end
    end
    
    
    for ii = 1:N*N
        mf = zeros(N_nb, 1);
        aF = zeros(N_nb, 1);
        aC = 0;
        bC = 0;
        
       %% Interior contributions
        for jj = 1:N_nb
            if nb_id(ii, jj) > 0
                %% Rhie-Chow interpolation
                Df = 0.5 * (D(ii, :) + D(nb_id(ii, jj), :));
%                 mf(jj) = getM(jj, S_nb, U, ii, N) - (2 - alphaU) * ((p(nb_id(ii, jj)) - p(ii)) * Df(ceil(jj / 2)) - 0.5 * (diag(Df) * (gradP(cellType{ii}, p, ii, N) + gradP(cellType{nb_id(ii, jj)}, p, nb_id(ii, jj), N)))' * S_nb(:, jj));
                aF(jj) = -Df(ceil(jj / 2));
                aC = aC - aF(jj);
                bC = bC - tempM(ii, jj); 
            else
%                 mf(jj) = 0;
                aF(jj) = 0;
            end
        end
        
        %% Add to matrices
        % diagnal
        A(ii, ii) = aC;
        % neighbors
        for jj = 1:N_nb
            if nb_id(ii, jj) > 0
                A(ii, nb_id(ii, jj)) = aF(jj);                
            end
        end
        % right-hand side
        B(ii) = bC;
        
    end
    
    %% set reference pressure
    refID = 1;
    A(refID, :) = 0;
    A(refID, refID) = 1;


    
    %% correct p
    A = sparse(A);
    pC = A \ B;
    pC(1) = 0;
    pNew = p + alphaP * pC;
    pNew(1) = 0;
    
    %% correct U
    for ii = 1:N*N     
        UNew(ii, :) = U(ii, :) - (diag(D(ii, :)) * gradPC(cellType{ii}, pC, ii, N))'; 
    end
end