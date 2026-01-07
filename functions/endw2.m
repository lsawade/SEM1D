function weight = endw2(n, alpha, beta)
% ENDW2 Compute end weight 2 for Gauss-Lobatto-Jacobi quadrature
% Based on Fortran implementation from gll_library.f90

    apb = alpha + beta;

    if n == 0
        weight = 0.0;
        return;
    end

    f1 = gammaf(alpha+1)*gammaf(beta+2)/gammaf(apb+3);
    f1 = f1*(apb+2)*2^(apb+2)/2;

    if n == 1
        weight = f1;
        return;
    end

    fint1 = gammaf(alpha+1)*gammaf(beta+2)/gammaf(apb+3);
    fint1 = fint1*2^(apb+2);
    fint2 = gammaf(alpha+2)*gammaf(beta+2)/gammaf(apb+4);
    fint2 = fint2*2^(apb+3);
    f2 = (2*(alpha+2)*fint1 - (apb+4)*fint2) * (apb+3)/4;

    if n == 2
        weight = f2;
        return;
    end

    for i = 3:n
        di = double(i-1);
        abn = alpha + beta + di;
        abnn = abn + di;
        a1 = -(2*(di+alpha)*(di+beta))/(abn*abnn*(abnn+1));
        a2 = (2*(alpha-beta))/(abnn*(abnn+2));
        a3 = (2*(abn+1))/((abnn+2)*(abnn+1));
        f3 = -(a2*f2+a1*f1)/a3;
        f1 = f2;
        f2 = f3;
    end

    weight = f3;
end
