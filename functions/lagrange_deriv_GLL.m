function deriv = lagrange_deriv_GLL(I, j, zgll)
% LAGRANGE_DERIV_GLL Compute derivative of I-th Lagrange interpolant
% at GLL point j
%
% Inputs:
%   I - index of Lagrange polynomial (1-based)
%   j - index of GLL point at which to evaluate (1-based)
%   zgll - GLL points (column vector)
%
% Output:
%   deriv - derivative value
%
% Based on Fortran implementation from lagrange_poly.f90
% Note: Original Fortran uses 0-based indexing

    nz = length(zgll);
    degpoly = nz - 1;

    % Convert to 0-based indexing for comparison with original
    I0 = I - 1;
    j0 = j - 1;

    if (I0 == 0) && (j0 == 0)
        deriv = -degpoly*(degpoly+1.0) / 4.0;
    elseif (I0 == degpoly) && (j0 == degpoly)
        deriv = degpoly*(degpoly+1.0) / 4.0;
    elseif I0 == j0
        deriv = 0.0;
    else
        deriv = pnleg(zgll(j), degpoly) / ...
                (pnleg(zgll(I), degpoly)*(zgll(j)-zgll(I))) ...
                + (1.0-zgll(j)*zgll(j))*pndleg(zgll(j), degpoly) / ...
                (degpoly*(degpoly+1.0)*pnleg(zgll(I), degpoly)* ...
                (zgll(j)-zgll(I))*(zgll(j)-zgll(I)));
    end
end
