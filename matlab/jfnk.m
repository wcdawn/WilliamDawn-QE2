function [flux, keff, x, k_conv_plot, x_conv_plot] = ...
    jfnk(A, xs, tol, maxit)

%     D1 = xs.D1;
%     D2 = xs.D2;
%     sigma_21 = xs.sigma_21;
%     sigma_a1 = xs.sigma_a1;
%     sigma_a2 = xs.sigma_a2;
    nusf1 = xs.nusf1;
    nusf2 = xs.nusf2;
%     sigma_r1 = xs.sigma_r1;
%     sigma_r2 = xs.sigma_r2;

    eta = 0.1;
    max_inner = 300;

    % array dimensions
    G = 2;
    N = size(A,1)/(5*G);
    matrix_idx = @(coeff,n,g)(coeff + 5*(n-1) + 5*N*(g-1));

    flux = ones(N,G);
    keff = 1.0;
    keff_old = 1.1;

    fission_sum_old = 1.0;

    x = ones(5*N*G,1);
    xm = [x; keff];

    k_conv_plot = [];
    x_conv_plot = [];

    tic;

    for m = 1:maxit

        R = residfun(A, xs, xm, keff_old, fission_sum_old);
        myresid = @(xm)(residfun(A, xs, xm, keff_old, fission_sum_old));
        matvec = @(v)(dirder(myresid, xm, v));
        % deltax = my_bicgstab(matvec, -R, eta, max_inner);
        deltax = gmres(matvec, -R, max_inner, eta, max_inner);
        xm = xm + deltax;

        % parse solution
        x = xm(1:5*N*G);
        kconv = abs(keff_old-keff);
        keff_old = keff;
        keff = xm(end);

        % parse solution for phi
        for n = 1:N
            for g = 1:G
                flux(n,g) = x(matrix_idx(1,n,g));
            end % G
        end % N
        % construct fission source
        fission_source = zeros(N,1);
        for n = 1:N
            fission_source(n) = nusf1*flux(n,1) + nusf2*flux(n,2);
        end
        fission_sum_old = sum(fission_source);


        R = residfun(A, xs, xm, keff_old, fission_sum_old);
        conv = norm(R,2);
        fprintf('it=%5d, keff=%.6f, kconv=%.6e, conv=%.6e\n',m,keff,kconv,conv);
        k_conv_plot = [k_conv_plot, kconv];
        x_conv_plot = [x_conv_plot, conv];

        if (conv < tol)
            break
        end
    end

    toc;

    if (m == maxit)
        warning('in jfnk() m == maxit');
    end
    fprintf('jfnk m = %d\n',m);

end

function dd = dirder(this_resid, x, v)

    if (norm(v) == 0d0)
        dd = zeros(size(v));
        return
    end

    fd_step = sqrt(eps);
    xs = (x.'*v)/norm(v,2);
    if (xs ~= 0d0)
        fd_step = fd_step*max([abs(xs),1d0])*sign(xs);
    end
    fd_step = fd_step/norm(v,2);
    dd = (this_resid(x+fd_step*v) - this_resid(x)) / fd_step;

end
