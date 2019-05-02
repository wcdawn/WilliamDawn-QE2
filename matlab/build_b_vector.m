function b = build_b_vector(xs, x, N);

%     D1 = xs.D1;
%     D2 = xs.D2;
%     sigma_21 = xs.sigma_21;
%     sigma_a1 = xs.sigma_a1;
%     sigma_a2 = xs.sigma_a2;
    nusf1 = xs.nusf1;
    nusf2 = xs.nusf2;
%     sigma_r1 = xs.sigma_r1;
%     sigma_r2 = xs.sigma_r2;

    matrix_idx = @(coeff,n,g)(coeff + 5*(n-1) + 5*N*(g-1));

    % build b vector
    b = zeros(size(x));
    for n = 1:N
        % flux continuity g = 1
        g = 1;
        b(matrix_idx(1,n,g)) = 0.0;

        % flux continuity g = 2
        g = 2;
        b(matrix_idx(1,n,g)) = 0.0;
        
        % current continuity g = 1
        g = 1;
        b(matrix_idx(2,n,g)) = 0.0;

        % current continuity g = 1
        g = 2;
        b(matrix_idx(2,n,g)) = 0.0;

        % weighed residual w_p = w_0 = 1 (nodal balance) g = 1
        g = 1;
        b(matrix_idx(3,n,g)) = (nusf1*x(matrix_idx(1,n,1)) + ...
            nusf2*x(matrix_idx(1,n,2)));

        % weighed residual w_p = w_0 = 1 (nodal balance) g = 2
        g = 2;
        b(matrix_idx(3,n,g)) = 0.0;

        % weighed residual w_p = w_1 g = 1
        g = 1;
        b(matrix_idx(4,n,g)) = (nusf1*(x(matrix_idx(2,n,1))/12.0 + ...
            x(matrix_idx(4,n,1))/120.0) + nusf2*(x(matrix_idx(2,n,2))/12.0 + ...
            x(matrix_idx(4,n,2))/120.0));

        % weighed residual w_p = w_1 g = 2
        g = 2;
        b(matrix_idx(4,n,g)) = 0.0;

        % weighed residual w_p = w_2 g = 1
        g = 1;
        b(matrix_idx(5,n,g)) = (nusf1*(x(matrix_idx(3,n,1))*0.05 + ...
            x(matrix_idx(5,n,1))/700.0) + nusf2*(x(matrix_idx(3,n,2))*0.05 + ...
            x(matrix_idx(5,n,2))/700.0));

        % weighed residual w_p = w_2 g = 2
        g = 2;
        b(matrix_idx(5,n,g)) = 0.0;
    end % N

end
