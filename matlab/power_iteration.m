function [flux, keff, x, k_conv_plot, x_conv_plot, time_plot] = ...
    power_iteration(A, xs, maxout, epspow, epsk)

    [L, U] = lu(A);

%     D1 = xs.D1;
%     D2 = xs.D2;
%     sigma_21 = xs.sigma_21;
%     sigma_a1 = xs.sigma_a1;
%     sigma_a2 = xs.sigma_a2;
    nusf1 = xs.nusf1;
    nusf2 = xs.nusf2;
%     sigma_r1 = xs.sigma_r1;
%     sigma_r2 = xs.sigma_r2;

    % array dimensions
    G = 2;
    N = size(A,1)/(5*G);
    matrix_idx = @(coeff,n,g)(coeff + 5*(n-1) + 5*N*(g-1));

    % initialization 
    flux = ones(N,G);
    keff = 1.0;

    tic;
    x = ones(5*N*G,1);
    keff_old = 1.1;
    fission_sum_old = 1.0;
    k_conv_plot = [];
    x_conv_plot = [];
    time_plot = [];

    for it = 1:maxout

        b = build_b_vector(xs, x, N);

        b = b / keff;

        x = U\(L\b);

        % parse solution for phi
        for n = 1:N
            for g = 1:G
                flux(n,g) = x(matrix_idx(1,n,g));
            end % G
        end % N

        fission_source = zeros(N,1);
        for n = 1:N
            fission_source(n) = nusf1*flux(n,1) + nusf2*flux(n,2);
        end
        fission_sum = sum(fission_source);

        keff = keff_old * fission_sum/fission_sum_old;

        tolx = norm(fission_sum(:) - fission_sum_old(:),inf);
        tolk = abs(keff - keff_old);
        keff_old = keff;
        fission_sum_old = fission_sum;

        fprintf('it=%5d, keff=%.6f, tolx=%.6e, tolk=%.6e\n',it,keff,tolx,tolk);
        k_conv_plot = [k_conv_plot, tolk];
%         x_conv_plot = [x_conv_plot, tolx];
        x_conv_plot = [x_conv_plot, ...
            norm(residfun(A, xs, [x;keff], keff_old, fission_sum_old),2)];

        power_it_time = toc;
        time_plot = [time_plot, power_it_time];

        if ((tolx <= epspow) && (tolk <= epsk)) 
            fprintf('Terminated Successfully\n');
            break
        end

    end % main iteration loop

    toc;

    fprintf('final residual = %.6e\n', ...
        norm(residfun(A, xs, [x;keff], keff_old, fission_sum_old)));

    
end
