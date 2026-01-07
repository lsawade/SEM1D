function [poly, pder, polym1, pderm1, polym2, pderm2] = jacobf(n, alp, bet, x)
% JACOBF Computes the Jacobi polynomial of degree n and its derivative at x
% Based on Fortran implementation from gll_library.f90
%
% Inputs:
%   n - degree of polynomial
%   alp - alpha parameter
%   bet - beta parameter
%   x - point at which to evaluate
%
% Outputs:
%   poly - polynomial value
%   pder - derivative value
%   polym1 - polynomial of degree n-1
%   pderm1 - derivative of polynomial of degree n-1
%   polym2 - polynomial of degree n-2
%   pderm2 - derivative of polynomial of degree n-2

    apb = alp + bet;
    poly = 1.0;
    pder = 0.0;
    psave = 0.0;
    pdsave = 0.0;

    if n == 0
        polym1 = 0;
        pderm1 = 0;
        polym2 = 0;
        pderm2 = 0;
        return;
    end

    polyl = poly;
    pderl = pder;
    poly = (alp-bet+(apb+2)*x)/2;
    pder = (apb+2)/2;

    if n == 1
        polym1 = polyl;
        pderm1 = pderl;
        polym2 = 0;
        pderm2 = 0;
        return;
    end

    for k = 2:n
        dk = double(k);
        a1 = 2*dk*(dk+apb)*(2*dk+apb-2);
        a2 = (2*dk+apb-1)*(alp^2-bet^2);
        b3 = (2*dk+apb-2);
        a3 = b3*(b3+1)*(b3+2);
        a4 = 2*(dk+alp-1)*(dk+bet-1)*(2*dk+apb);
        polyn = ((a2+a3*x)*poly-a4*polyl)/a1;
        pdern = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1;
        psave = polyl;
        pdsave = pderl;
        polyl = poly;
        poly = polyn;
        pderl = pder;
        pder = pdern;
    end

    polym1 = polyl;
    pderm1 = pderl;
    polym2 = psave;
    pderm2 = pdsave;
end
