function [x] = my_bicgstab(matvec, b, tol, maxit)

    % algorithm from Kelley "Iterative Methods for Linear and Nonlinear
    % Equations" (red book), Algorithm 3.6.3, p. 50.

    % I've written this to work for general matvec().
    % If we have an explicit form of A,
    % matvec = @(x)(A*x);

    % lots of initialization

    % use zero vector as initial iterate
    N = length(b);
    x = zeros(N,1);
    r = b;
    % otherwise
    % if (x0 is present)
    %     x = x0;
    %     r = -matvec(x) - b;
    % end

    r0_hat = r;

    rho_0 = 1.0;
    alpha = rho_0;
    omega = rho_0;

    v = zeros(N,1);
    p = v;

    rho_1 = r0_hat'*r;

    conv = tol*norm(b,2);

    for it = 1:maxit
        if (omega == 0) 
            disp(it);
            disp(rnorm);
            error('BiCGSTAB breakdown, omega=0');
        end

        % (b)
        beta = (rho_1/rho_0) * (alpha/omega);
        % (c)
        p = r + beta*(p - omega*v);
        % (d)
        v = matvec(p);
        tau = r0_hat'*v;
        if (tau == 0)
            error('BiCGSTAB breakdown, tau=0');
        end
        % (e)
        alpha = rho_1 / tau;
        % (f)
        s = r - alpha*v;
        t = matvec(s);
        tau = t'*t;
        if (tau == 0)
            error('BiCGSTAB breakdown, t=0');
        end
        % (g)
        omega = (t'*s) / tau;
        rho_0 = rho_1;
        rho_1 = -omega * (r0_hat'*t);
        % (h)
        x = x + alpha*p + omega*s;
        % (i)
        r = s - omega*t;
        rnorm = norm(r,2);

        if (rnorm < conv)
            break
        end
    end

    fprintf('it = %d\n',it);
    if (it == maxit)
        warning('in bicgstab() it == maxit');
    end
end
