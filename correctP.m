% This is a function to solve pressure
% correction equation and to correct
% the pressure field and velocity field

function [pNew, UNew] = correctP(U, p, D, N)
global L mu h ub alphaU alphaP

    %% Matrices and vectors
    A = zeros(N * N, N * N);
    B = zeros(N * N, 1);
    UNew = U;
    
    %% Neighbor surfaces
    S_nb = [-h, h, 0, 0; 0, 0, -h, h];  % Surface vector for w, e, s, n
    N_nb = length(S_nb(1, :));
    
    for ii = 1:N*N
        cellType = getCellType(ii, N);
        mf = zeros(N_nb, 1);
        aF = zeros(N_nb, 1);
        aC = 0;
        bC = 0;
        nb_id = [ii - 1, ii + 1, ii - N, ii + N];
        
        %% Boundary contributions
        if contains(cellType, 'L')
            nb_id(1) = 0;
        end
        if contains(cellType, 'R')
            nb_id(2) = 0;
        end
        if contains(cellType, 'D')
            nb_id(3) = 0;
        end
        if contains(cellType, 'U')
            nb_id(4) = 0;
        end
        
       %% Interior contributions
        for jj = 1:N_nb
            if nb_id(jj) > 0
                %% Rhie-Chow interpolation
                Df = 0.5 * (D(ii, :) + D(nb_id(jj), :));
                mf(jj) = getM(jj, S_nb, U, ii, N) - (2 - alphaU) * ((p(nb_id(jj)) - p(ii)) * Df(ceil(jj / 2)) - 0.5 * (diag(Df) * (gradP(p, ii, N) + gradP(p, nb_id(jj), N)))' * S_nb(:, jj));
                aF(jj) = -Df(ceil(jj / 2));
                aC = aC - aF(jj);
                bC = bC - mf(jj); 
            else
                mf(jj) = 0;
                aF(jj) = 0;
            end
        end
        
        %% Add to matrices
        % diagnal
        A(ii, ii) = aC;
        % neighbors
        for jj = 1:N_nb
            if nb_id(jj) > 0
                A(ii, nb_id(jj)) = aF(jj);                
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
        cellType = getCellType(ii, N);
        nb_id = [ii - 1, ii + 1, ii - N, ii + N];        
        UNew(ii, :) = U(ii, :) - (diag(D(ii, :)) * gradPC(pC, ii, N))'; 
    end
end