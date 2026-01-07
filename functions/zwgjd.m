function [z, w] = zwgjd(np, alpha, beta)
% ZWGJD Generate np Gauss-Jacobi points and weights
% associated with Jacobi polynomial of degree n = np-1
%
% Note: Coefficients alpha and beta must be greater than -1
% Based on Fortran implementation from gll_library.f90

    if np <= 0
        error('minimum number of Gauss points is 1');
    end

    if (alpha <= -1) || (beta <= -1)
        error('alpha and beta must be greater than -1');
    end

    n = np - 1;
    apb = alpha + beta;

    if np == 1
        z = (beta-alpha)/(apb+2);
        w = gammaf(alpha+1)*gammaf(beta+1)/gammaf(apb+2) * 2^(apb+1);
        return;
    end

    z = jacg(np, alpha, beta);
    w = zeros(np, 1);

    np1 = n + 1;
    np2 = n + 2;
    dnp1 = double(np1);
    dnp2 = double(np2);
    fac1 = dnp1 + alpha + beta + 1;
    fac2 = fac1 + dnp1;
    fac3 = fac2 + 1;
    fnorm = pnormj(np1, alpha, beta);
    rcoef = (fnorm*fac2*fac3)/(2*fac1*dnp2);

    for i = 1:np
        [p, pd, ~, pdm1, ~, ~] = jacobf(np2, alpha, beta, z(i));
        w(i) = -rcoef/(p*pdm1);
    end
end
