function R = residfun(A, xs, xm, keff_old, fission_sum_old);
    
%     D1 = xs.D1;
%     D2 = xs.D2;
%     sigma_21 = xs.sigma_21;
%     sigma_a1 = xs.sigma_a1;
%     sigma_a2 = xs.sigma_a2;
    nusf1 = xs.nusf1;
    nusf2 = xs.nusf2;
%     sigma_r1 = xs.sigma_r1;
%     sigma_r2 = xs.sigma_r2;

    G = 2;
    N = size(A,1)/(5*G);
    matrix_idx = @(coeff,n,g)(coeff + 5*(n-1) + 5*N*(g-1));

    keff = xm(end);
    x = xm(1:5*N*G);

    flux = zeros(N,G);
    for n = 1:N
        for g = 1:G
            flux(n,g) = x(matrix_idx(1,n,g));
        end
    end

    % construct fission source
    fission_source = zeros(N,1);
    for n = 1:N
        fission_source(n) = nusf1*flux(n,1) + nusf2*flux(n,2);
    end
    fission_sum = sum(fission_source);

    b = build_b_vector(xs, x, N);

    R = zeros(5*N*G+1,1);
    R(1:5*N*G) = A*xm(1:5*N*G) - b/keff;
    R(5*N*G+1) = keff - keff_old*fission_sum/fission_sum_old;

end
