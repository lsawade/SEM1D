function [xigll, wgll, hprime] = define_derivative_matrix(NGLL)
% DEFINE_DERIVATIVE_MATRIX Set up GLL points, weights, and derivative matrix
%
% Input:
%   NGLL - number of GLL points (polynomial degree plus one)
%
% Outputs:
%   xigll - Gauss-Lobatto-Legendre points of integration
%   wgll - weights
%   hprime - array with derivatives of Lagrange polynomials
%            hprime(i,j) = h'_i(xigll_j)
%
% Based on Fortran implementation from define_derivative_matrix.f90

    % Parameters for Gauss-Lobatto-Legendre points
    GAUSSALPHA = 0.0;
    GAUSSBETA = 0.0;

    % Set up coordinates of the Gauss-Lobatto-Legendre points
    [xigll, wgll] = zwgljd(NGLL, GAUSSALPHA, GAUSSBETA);

    % If number of points is odd, the middle abscissa is exactly zero
    if mod(NGLL, 2) ~= 0
        xigll((NGLL-1)/2 + 1) = 0.0;
    end

    % Calculate derivatives of the Lagrange polynomials
    % hprime(i,j) = h'_i(xigll_j) by definition of the derivative matrix
    hprime = zeros(NGLL, NGLL);
    for i1 = 1:NGLL
        for i2 = 1:NGLL
            hprime(i1, i2) = lagrange_deriv_GLL(i1, i2, xigll);
        end
    end
end
