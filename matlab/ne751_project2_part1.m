clear
close all

% geometry input
Lx = 200.0; % [cm]
N = 50; % number of nodes
G = 2; % ngroup
albedo = 0.0;

% plot properties
FN = 'Times New Roman';
FS = 12;
LW = 2;

% tolerance/iteration control
maxout = 200;
epsk = 1d-10;
epspow = 1d-9;

% material B
D1 = 1.4291; % [cm]
D2 = 0.4453; % [cm]
sigma_21 = 0.0168; % [1/cm]
sigma_a1 = 0.0095; % [1/cm]
sigma_a2 = 0.0745; % [1/cm]
nusf1 = 0.0058; % [neut/cm]
nusf2 = 0.1033; % [neut/cm]

% calcualte removal xs
sigma_r1 = sigma_a1 + sigma_21;
sigma_r2 = sigma_a2;

% build a mesh
h = Lx/N;
node = h*0.5:h:Lx;

% calculate analytic results
buckling = (pi/Lx)^2; % note this is actually B^2
c = [1.0, sigma_21/(D2*buckling+sigma_r2)]; % my thesis Eq. (A.89)
phi_analytic = @(x,g)(c(g)*sin(sqrt(buckling)*x));
keff_analytic = (nusf1 + nusf2*c(2))/(D1*buckling+sigma_r1);
cinf = sigma_21/sigma_r2;
kinf_analytic = (nusf1 + nusf2*cinf)/(sigma_r1);
fprintf('Analytic Results\n');
fprintf('k2/k1 = %.16f\n',c(2));
fprintf('keff  = %.16f\n',keff_analytic);
fprintf('cinf  = %.10f\n',cinf);
fprintf('kinf  = %.10f\n',kinf_analytic);

% build A matrix
% A = zeros(5*N*G); % square, dense matrix
A = sparse(5*N*G); % sparse is actually faster
matrix_idx = @(coeff,n,g)(coeff + 5*(n-1) + 5*N*(g-1));
for n = 1:N
    % flux continuity g = 1
    g = 1;
    row = matrix_idx(1,n,g);
    if (n == 1)
        A(row,matrix_idx(1,n,g)) = 1.0;
        A(row,matrix_idx(2,n,g)) = 0.5;
        A(row,matrix_idx(3,n,g)) = 0.5;
        A(row,matrix_idx(1,n+1,g)) = -1.0;
        A(row,matrix_idx(2,n+1,g)) = +0.5;
        A(row,matrix_idx(3,n+1,g)) = -0.5;
    elseif (n == N)
        % zero flux at x=0
        A(row,matrix_idx(1,1,g)) = 1.0;
        A(row,matrix_idx(2,1,g)) = -0.5;
        A(row,matrix_idx(3,1,g)) = 0.5;
    else
        A(row,matrix_idx(1,n,g)) = 1.0;
        A(row,matrix_idx(2,n,g)) = 0.5;
        A(row,matrix_idx(3,n,g)) = 0.5;
        A(row,matrix_idx(1,n+1,g)) = -1.0;
        A(row,matrix_idx(2,n+1,g)) = +0.5;
        A(row,matrix_idx(3,n+1,g)) = -0.5;
    end

    % flux continuity g = 2
    g = 2;
    row = matrix_idx(1,n,g);
    if (n == 1)
        A(row,matrix_idx(1,n,g)) = 1.0;
        A(row,matrix_idx(2,n,g)) = 0.5;
        A(row,matrix_idx(3,n,g)) = 0.5;
        A(row,matrix_idx(1,n+1,g)) = -1.0;
        A(row,matrix_idx(2,n+1,g)) = +0.5;
        A(row,matrix_idx(3,n+1,g)) = -0.5;
    elseif (n == N)
        % zero flux at x=0
        A(row,matrix_idx(1,1,g)) = 1.0;
        A(row,matrix_idx(2,1,g)) = -0.5;
        A(row,matrix_idx(3,1,g)) = 0.5;
    else
        A(row,matrix_idx(1,n,g)) = 1.0;
        A(row,matrix_idx(2,n,g)) = 0.5;
        A(row,matrix_idx(3,n,g)) = 0.5;
        A(row,matrix_idx(1,n+1,g)) = -1.0;
        A(row,matrix_idx(2,n+1,g)) = +0.5;
        A(row,matrix_idx(3,n+1,g)) = -0.5;
    end

    % current continuity g = 1
    g = 1;
    % diffusion coefficients are uniform in this problem but not in general
    Dnh  = D1/h; % D^{n}/h
    Dn1h = D1/h; % D^{n+1}/h
    row = matrix_idx(2,n,g);
    if (n == 1)
        A(row,matrix_idx(2,n,g)) = Dnh;
        A(row,matrix_idx(3,n,g)) = 3.0*Dnh;
        A(row,matrix_idx(4,n,g)) = -0.5*Dnh;
        A(row,matrix_idx(5,n,g)) = -0.2*Dnh;
        A(row,matrix_idx(2,n+1,g)) = -Dn1h;
        A(row,matrix_idx(3,n+1,g)) = 3.0*Dn1h;
        A(row,matrix_idx(4,n+1,g)) = 0.5*Dn1h;
        A(row,matrix_idx(5,n+1,g)) = -0.2*Dn1h;
    elseif (n == N)
        % zero flux at x=Lx
        A(row,matrix_idx(1,N,g)) = 1.0;
        A(row,matrix_idx(2,N,g)) = 0.5;
        A(row,matrix_idx(3,N,g)) = 0.5;
    else
        A(row,matrix_idx(2,n,g)) = Dnh;
        A(row,matrix_idx(3,n,g)) = 3.0*Dnh;
        A(row,matrix_idx(4,n,g)) = -0.5*Dnh;
        A(row,matrix_idx(5,n,g)) = -0.2*Dnh;
        A(row,matrix_idx(2,n+1,g)) = -Dn1h;
        A(row,matrix_idx(3,n+1,g)) = 3.0*Dn1h;
        A(row,matrix_idx(4,n+1,g)) = 0.5*Dn1h;
        A(row,matrix_idx(5,n+1,g)) = -0.2*Dn1h;
    end

    % current continuity g = 2
    g = 2;
    % diffusion coefficients are uniform in this problem but not in general
    Dnh = D2/h; % D^{n}/h
    Dn1h = D2/h; % D^{n+1}/h
    row = matrix_idx(2,n,g);
    if (n == 1)
        A(row,matrix_idx(2,n,g)) = Dnh;
        A(row,matrix_idx(3,n,g)) = 3.0*Dnh;
        A(row,matrix_idx(4,n,g)) = -0.5*Dnh;
        A(row,matrix_idx(5,n,g)) = -0.2*Dnh;
        A(row,matrix_idx(2,n+1,g)) = -Dn1h;
        A(row,matrix_idx(3,n+1,g)) = 3.0*Dn1h;
        A(row,matrix_idx(4,n+1,g)) = 0.5*Dn1h;
        A(row,matrix_idx(5,n+1,g)) = -0.2*Dn1h;
    elseif (n == N)
        % zero flux at x=Lx
        A(row,matrix_idx(1,N,g)) = 1.0;
        A(row,matrix_idx(2,N,g)) = 0.5;
        A(row,matrix_idx(3,N,g)) = 0.5;
    else
        A(row,matrix_idx(2,n,g)) = Dnh;
        A(row,matrix_idx(3,n,g)) = 3.0*Dnh;
        A(row,matrix_idx(4,n,g)) = -0.5*Dnh;
        A(row,matrix_idx(5,n,g)) = -0.2*Dnh;
        A(row,matrix_idx(2,n+1,g)) = -Dn1h;
        A(row,matrix_idx(3,n+1,g)) = 3.0*Dn1h;
        A(row,matrix_idx(4,n+1,g)) = 0.5*Dn1h;
        A(row,matrix_idx(5,n+1,g)) = -0.2*Dn1h;
    end

    % weighed residual w_p = w_0 = 1 (nodal balance) g = 1
    g = 1;
    Dn = D1;
    row = matrix_idx(3,n,g);
    A(row,matrix_idx(1,n,g)) = sigma_r1;
    A(row,matrix_idx(3,n,g)) = -Dn/h^2 * 6.0;
    A(row,matrix_idx(5,n,g)) = Dn/h^2 * 0.4;

    % weighed residual w_p = w_0 = 1 (nodal balance) g = 2
    g = 2;
    Dn = D2;
    row = matrix_idx(3,n,g);
    A(row,matrix_idx(1,n,g)) = sigma_r2;
    A(row,matrix_idx(3,n,g)) = -Dn/h^2 * 6.0;
    A(row,matrix_idx(5,n,g)) = Dn/h^2 * 0.4;
    A(row,matrix_idx(1,n,1)) = -sigma_21;

    % weighed residual w_p = w_1 g = 1
    g = 1;
    Dn = D1;
    row = matrix_idx(4,n,g);
    A(row,matrix_idx(2,n,g)) = sigma_r1/12.0;
    A(row,matrix_idx(4,n,g)) = 0.5*Dn/h^2 + sigma_r1/120.0;

    % weighed residual w_p = w_1 g = 2
    g = 2;
    Dn = D2;
    row = matrix_idx(4,n,g);
    A(row,matrix_idx(2,n,g)) = sigma_r2/12.0;
    A(row,matrix_idx(4,n,g)) = 0.5*Dn/h^2 + sigma_r2/120.0;
    A(row,matrix_idx(2,n,1)) = -sigma_21/12.0;
    A(row,matrix_idx(4,n,1)) = -sigma_21/120.0;

    % weighed residual w_p = w_2 g = 1
    g = 1;
    Dn = D1;
    row = matrix_idx(5,n,g);
    A(row,matrix_idx(3,n,g)) = sigma_r1*0.05;
    A(row,matrix_idx(5,n,g)) = 0.2*Dn/h^2 + sigma_r2/700.0;

    % weighed residual w_p = w_2 g = 2
    g = 2;
    Dn = D2;
    row = matrix_idx(5,n,g);
    A(row,matrix_idx(3,n,g)) = sigma_r2*0.05;
    A(row,matrix_idx(5,n,g)) = 0.2*Dn/h^2 + sigma_r2/700.0;
    A(row,matrix_idx(3,n,1)) = -sigma_21*0.05;
    A(row,matrix_idx(5,n,1)) = -sigma_21/700.0;

end

figure
spy(A)
title('Sparsity Pattern of A 5NG \times 5NG Matrix');
set(gca,'FontName',FN,'FontSize',FS);

% setup data structure
xs.D1 = D1;
xs.D2 = D2;
xs.sigma_21 = sigma_21;
xs.sigma_a1 = sigma_a1;
xs.sigma_a2 = sigma_a2;
xs.nusf1 = nusf1;
xs.nusf2 = nusf2;
xs.sigma_r1 = sigma_r1;
xs.sigma_r2 = sigma_r2;

[flux, keff, x, k_conv_plot, pi_conv] = power_iteration(A, xs, maxout, epspow, epsk);
[flux, keff, x, k_conv_plot, jfnk_conv] = jfnk(A, xs, epspow, maxout);

figure
semilogy(pi_conv,'-*');
hold on
semilogy(jfnk_conv,'-*');
hold off
legend({'Power Iteration','JFNK'});
xlabel('Iteration');
ylabel('Convergence Criteria');

if (length(k_conv_plot) == maxout)
    fprintf('MAX OUTERS EXCEEDED!\n');
end

flux_norm_const = max(abs(flux(:,1)));
flux(:,:) = flux(:,:) / flux_norm_const;

fprintf('keff_computed  = %.6f\n',keff);
fprintf('kerr           = %.6e [pcm]\n',(keff_analytic-keff)*1e5);
fprintf('ratio_computed = %.6f\n',max(flux(:,2)));
fprintf('ratio_err      = %.6e\n',c(2)-max(flux(:,2)));
fprintf('inf_norm_1     = %.6e\n',norm(flux(:,1)-phi_analytic(node,1)',inf));
fprintf('inf_norm_2     = %.6e\n',norm(flux(:,2)-phi_analytic(node,2)',inf));

rt_edges = h:h:Lx;
current = zeros(N,G);
for n = 1:N
    for g = 1:G
        if (g == 1)
            Dnh = D1/h;
        else
            Dnh = D2/h;
        end
        current(n,g) = -Dnh * (x(matrix_idx(2,n,g)) + 3*x(matrix_idx(3,n,g)) - ...
            0.5*x(matrix_idx(4,n,g)) - 0.2*x(matrix_idx(5,n,g)));
    end
end

% figure
% semilogy(k_conv_plot);
% hold on
% semilogy(x_conv_plot);
% hold off
% xlabel('Iteration');
% ylabel('Convergence Crieteria');
% legend({'tolk','tolx'});

% figure
% h = plot(node,phi_analytic(node,1),'-o','LineWidth',LW);
% set(h,'MarkerFaceColor',get(h,'Color'));
% hold on
% h = plot(node,phi_analytic(node,2),'-^','LineWidth',LW);
% set(h,'MarkerFaceColor',get(h,'Color'));
% hold off
% xlim([0,Lx])
% xlabel('x [cm]');
% ylabel('\phi_g');
% legend({'\phi_1','\phi_2'});
% title('Analytic Solution');
% set(gca,'FontName',FN,'FontSize',FS);

if (N < 30)
    label_circle = '-o';
    label_triangle = '-^';
else
    label_circle = '-';
    label_triangle = '--';
end
figure
h = plot(node,flux(:,1),label_circle,'LineWidth',LW);
set(h,'MarkerFaceColor',get(h,'Color'));
hold on
h = plot(node,flux(:,2),label_triangle,'LineWidth',LW);
set(h,'MarkerFaceColor',get(h,'Color'));
hold off
xlim([node(1),node(end)])
xlabel('x [cm]');
ylabel('\phi_g');
legend({'\phi_1','\phi_2'});
title('Numeric/Computed Flux Solution');
set(gca,'FontName',FN,'FontSize',FS);
print(gcf,'./figs/flux_solution.eps','-depsc2');
%close(gcf);

figure
h = plot(node,current(:,1),label_circle,'LineWidth',LW);
set(h,'MarkerFaceColor',get(h,'Color'));
hold on
h = plot(node,current(:,2),label_triangle,'LineWidth',LW);
set(h,'MarkerFaceColor',get(h,'Color'));
hold off
xlim([node(1),node(end)])
xlabel('x [cm]');
ylabel('J_g');
legend({'J_1','J_2'});
title('Numeric/Computed Current Solution');
set(gca,'FontName',FN,'FontSize',FS);
print(gcf,'./figs/current_solution.eps','-depsc2');
%close(gcf);

figure
yyaxis left
plot(node,flux(:,1),'LineWidth',LW);
hold on
plot(node,phi_analytic(node,1)','LineWidth',LW);
ylabel('\phi(x)')
yyaxis right
plot(node,flux(:,1)-phi_analytic(node,1)','LineWidth',LW)
ylabel('Error');
xlabel('x [cm]');
title('\phi_1 Error');
legend({'NEM','Exact','Absolute Error'})
set(gca,'FontName',FN,'FontSize',FS);
print(gcf,'./figs/phi1_error.eps','-depsc2');
%close(gcf);

figure
yyaxis left
plot(node,flux(:,2),'LineWidth',LW);
hold on
plot(node,phi_analytic(node,2)','LineWidth',LW);
ylabel('\phi(x)')
yyaxis right
plot(node,(flux(:,2)-phi_analytic(node,2)'),'LineWidth',LW)
ylabel('Error');
xlabel('x [cm]');
title('\phi_2 Error');
legend({'NEM','Exact','Absolute Error'})
set(gca,'FontName',FN,'FontSize',FS);
print(gcf,'./figs/phi2_error.eps','-depsc2');
%close(gcf)
