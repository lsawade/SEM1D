function [z, w] = zwgljd(np, alpha, beta)
% ZWGLJD Generate np Gauss-Lobatto-Jacobi points and weights
% associated with Jacobi polynomials of degree n = np-1
%
% Note: alpha and beta coefficients must be greater than -1
%       Legendre polynomials are special case (alpha = beta = 0)
% Based on Fortran implementation from gll_library.f90

    if np <= 1
        error('minimum number of Gauss-Lobatto points is 2');
    end

    % with spectral elements, use at least 3 points
    if np <= 2
        error('minimum number of Gauss-Lobatto points for the SEM is 3');
    end

    if (alpha <= -1) || (beta <= -1)
        error('alpha and beta must be greater than -1');
    end

    n = np - 1;
    nm1 = n - 1;

    z = zeros(np, 1);
    w = zeros(np, 1);

    if nm1 > 0
        alpg = alpha + 1;
        betg = beta + 1;
        [z_inner, w_inner] = zwgjd(nm1, alpg, betg);
        z(2:np-1) = z_inner;
        w(2:np-1) = w_inner;
    end

    z(1) = -1.0;
    z(np) = 1.0;

    for i = 2:(np-1)
        w(i) = w(i)/(1-z(i)^2);
    end

    [p, pd, ~, ~, ~, ~] = jacobf(n, alpha, beta, z(1));
    w(1) = endw1(n, alpha, beta)/(2*pd);

    [p, pd, ~, ~, ~, ~] = jacobf(n, alpha, beta, z(np));
    w(np) = endw2(n, alpha, beta)/(2*pd);
end
