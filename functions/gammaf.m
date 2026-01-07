function gamma_val = gammaf(x)
% GAMMAF Gamma function for specific values
% Based on Fortran implementation from gll_library.f90

    gamma_val = 1.0;

    if x == -0.5
        gamma_val = -2*sqrt(pi);
    elseif x == 0.5
        gamma_val = sqrt(pi);
    elseif x == 1.0
        gamma_val = 1.0;
    elseif x == 2.0
        gamma_val = 1.0;
    elseif x == 1.5
        gamma_val = sqrt(pi)/2.0;
    elseif x == 2.5
        gamma_val = 1.5*sqrt(pi)/2.0;
    elseif x == 3.5
        gamma_val = 2.5*1.5*sqrt(pi)/2.0;
    elseif x == 3.0
        gamma_val = 2.0;
    elseif x == 4.0
        gamma_val = 6.0;
    elseif x == 5.0
        gamma_val = 24.0;
    elseif x == 6.0
        gamma_val = 120.0;
    end
end
