function norm = pnormj(n, alpha, beta)
% PNORMJ Compute normalization for Jacobi polynomials
% Based on Fortran implementation from gll_library.f90

    dn = double(n);
    const = alpha + beta + 1;

    if n <= 1
        prod = gammaf(dn+alpha)*gammaf(dn+beta);
        prod = prod/(gammaf(dn)*gammaf(dn+alpha+beta));
        norm = prod * 2^const/(2*dn+const);
        return;
    end

    prod = gammaf(alpha+1)*gammaf(beta+1);
    prod = prod/(2*(1+const)*gammaf(const+1));
    prod = prod*(1+alpha)*(2+alpha);
    prod = prod*(1+beta)*(2+beta);

    for i = 3:n
        dindx = double(i);
        frac = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta));
        prod = prod*frac;
    end

    norm = prod * 2^const/(2*dn+const);
end
