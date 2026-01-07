function pder = pndleg(z, n)
% PNDLEG Compute the derivative of the Nth order Legendre polynomial at Z
% Based on the recursion formula for the Legendre polynomials
% Based on Fortran implementation from gll_library.f90

    p1 = 1.0;
    p2 = z;
    p1d = 0.0;
    p2d = 1.0;
    p3d = 1.0;

    for k = 1:(n-1)
        fk = double(k);
        p3 = ((2*fk+1)*z*p2 - fk*p1)/(fk+1);
        p3d = ((2*fk+1)*p2 + (2*fk+1)*z*p2d - fk*p1d)/(fk+1);
        p1 = p2;
        p2 = p3;
        p1d = p2d;
        p2d = p3d;
    end

    pder = p3d;
end
