function xjac = jacg(np, alpha, beta)
% JACG Computes np Gauss-Jacobi points
% Based on Fortran implementation from gll_library.f90
%
% Inputs:
%   np - number of points
%   alpha - Jacobi parameter (alpha = beta = 0 -> Legendre points)
%   beta - Jacobi parameter (alpha = beta = -0.5 -> Chebyshev points)
%
% Outputs:
%   xjac - Gauss-Jacobi points

    K_MAX_ITER = 10;
    eps = 1.0e-12;

    xjac = zeros(np, 1);
    n = np - 1;
    dth = 4*atan(1.0)/(2*n+2);

    xlast = 0.0;

    for j = 1:np
        if j == 1
            x = cos((2*(j-1)+1)*dth);
        else
            x1 = cos((2*(j-1)+1)*dth);
            x2 = xlast;
            x = (x1+x2)/2;
        end

        for k = 1:K_MAX_ITER
            [p, pd, ~, ~, ~, ~] = jacobf(np, alpha, beta, x);
            recsum = 0.0;
            jm = j - 1;
            for i = 1:jm
                recsum = recsum + 1.0/(x-xjac(np-i+1));
            end
            delx = -p/(pd-recsum*p);
            x = x + delx;
            if abs(delx) < eps
                break;
            end
        end

        xjac(np-j+1) = x;
        xlast = x;
    end

    % Sort the points
    xjac = sort(xjac);
end
