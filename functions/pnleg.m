function poly = pnleg(z, n)
% PNLEG Compute the value of the Nth order Legendre polynomial at Z
% Based on the recursion formula for the Legendre polynomials
% Based on Fortran implementation from gll_library.f90

    p1 = 1.0;
    p2 = z;
    p3 = p2;

    for k = 1:(n-1)
        fk = double(k);
        p3 = ((2*fk+1)*z*p2 - fk*p1)/(fk+1);
        p1 = p2;
        p2 = p3;
    end

    poly = p3;
end
