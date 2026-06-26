%% Simplified SDCIB Eigenvalue and Load-Step Study
% Computes the equilibrium, small-signal eigenvalues, participation factors,
% and one load-step response for the reduced SDCIB model.
clear; clc; close all;

%% Parameters
p = struct();

% Base / grid (SDCIB)
p.fb      = 60;
p.omega_b = 2*pi*p.fb;      % rad/s
p.omega_s = 1.0;            % pu synchronous speed

% Infinite bus in ri-frame
p.Vinf = 1.0;
p.Rinf = 0.01;
p.Xinf = 0.5;

% AFE
p.omega_lp    = 2*pi*100;   % rad/s
p.l_afe       = 0.05;
p.r_afe       = 0.003;
p.c_dc        = 2.0;
p.vdc_ups_ref = 1.0;

% VSI
p.l_vsi       = 0.05;
p.r_vsi       = 0.003;
p.c_vsi       = 0.2;
p.vu_vsi_ref  = 1.0;
p.omega_vsi   = p.omega_s;

% PSU reduced equivalent
p.c_psu    = 2.0;
p.r_psu    = 0.005;
p.vpsu_ref = 1.0;

% Downstream DC-DC/load reduced equivalent
p.c_eq    = 0.2;
p.veq_ref = 0.5;

% Load input
p.p_load0 = 0.5;      % equilibrium load before the step (pu)
p.p_load1 = 0.6;      % load after the step (pu)

%% Bandwidth-based tuning targets (fbw, zeta)
tg = struct();

tg.pll.fbw      = 20.0;    tg.pll.zeta    = 0.707;
tg.dc_afe.fbw   = 10.0;     tg.dc_afe.zeta = 1;
tg.c_afe.fbw    = 400.0;   tg.c_afe.zeta  = 0.707;
tg.v_vsi.fbw    = 120;     tg.v_vsi.zeta  = 1.0;
tg.c_vsi.fbw    = 600.0;   tg.c_vsi.zeta  = 1.0;
tg.v_psu.fbw    = 39.0;    tg.v_psu.zeta  = 1.0;
tg.v_eq.fbw     = 60.0;    tg.v_eq.zeta   = 1.0;

%% Compute PI gains from (fbw, zeta)
[p.kp_dc_afe, p.ki_dc_afe] = pi_bw_tuning('voltage', tg.dc_afe.fbw, tg.dc_afe.zeta, p.c_dc,  p.omega_b);
[p.kp_c_afe,  p.ki_c_afe]  = pi_bw_tuning('current', tg.c_afe.fbw,  tg.c_afe.zeta,  p.l_afe, p.omega_b, p.r_afe);

[p.kp_v_vsi, p.ki_v_vsi] = pi_bw_tuning('voltage', tg.v_vsi.fbw, tg.v_vsi.zeta, p.c_vsi, p.omega_b);
[p.kp_c_vsi, p.ki_c_vsi] = pi_bw_tuning('current', tg.c_vsi.fbw, tg.c_vsi.zeta, p.l_vsi, p.omega_b, p.r_vsi);

[p.kp_v_psu, p.ki_v_psu] = pi_bw_tuning('voltage', tg.v_psu.fbw, tg.v_psu.zeta, p.c_psu, p.omega_b);
[p.kp_v_eq,  p.ki_v_eq]  = pi_bw_tuning('voltage', tg.v_eq.fbw,  tg.v_eq.zeta,  p.c_eq,  p.omega_b);

[p.kp_pll, p.ki_pll] = pll_bw_tuning_simple(tg.pll.fbw, tg.pll.zeta);

fprintf('\n=== Gains from bandwidth-based tuning ===\n');
fprintf('PLL:      kp=%.4g, ki=%.4g\n', p.kp_pll,    p.ki_pll);
fprintf('AFE-DC:   kp=%.4g, ki=%.4g\n', p.kp_dc_afe, p.ki_dc_afe);
fprintf('AFE-I:    kp=%.4g, ki=%.4g\n', p.kp_c_afe,  p.ki_c_afe);
fprintf('VSI-V:    kp=%.4g, ki=%.4g\n', p.kp_v_vsi,  p.ki_v_vsi);
fprintf('VSI-I:    kp=%.4g, ki=%.4g\n', p.kp_c_vsi,  p.ki_c_vsi);
fprintf('PSU-V:    kp=%.4g, ki=%.4g\n', p.kp_v_psu,  p.ki_v_psu);
fprintf('EQ-V:     kp=%.4g, ki=%.4g\n', p.kp_v_eq,   p.ki_v_eq);

%% State ordering
state_names = { ...
    'theta_pll','eps_pll','vq_pll','id_afe','iq_afe','xi_dc_afe','gamd_afe','gamq_afe', ...
    'vdc_ups', 'iU_cv','iV_cv','vU_vsi','vV_vsi','xiU_vsi','xiV_vsi','gamU_vsi','gamV_vsi', ...
    'v_psu','xi_psu','v_eq','xi_eq' };

n   = numel(state_names);
idx = cell2struct(num2cell(1:n), state_names, 2);

%% Initial guess
x = zeros(n,1);

x(idx.vdc_ups) = p.vdc_ups_ref;
x(idx.vU_vsi)  = p.vu_vsi_ref;
x(idx.vV_vsi)  = 0;
x(idx.v_psu)   = p.vpsu_ref;
x(idx.v_eq)    = p.veq_ref;

g_load0 = p.p_load0/(3*p.veq_ref^2);
i_eq0   = g_load0*p.veq_ref;
i_psu0  = (p.veq_ref/p.vpsu_ref)*i_eq0;

g_eq0 = i_psu0 / max(p.vu_vsi_ref,1e-6);
x(idx.xi_psu) = g_eq0 / p.ki_v_psu;
x(idx.xi_eq)  = i_eq0 / p.ki_v_eq;

x(idx.iU_cv) = g_eq0 * p.vu_vsi_ref;
x(idx.iV_cv) = 0;

x(idx.id_afe) = p.p_load0;
x(idx.iq_afe) = 0;

%% Newton equilibrium
fprintf('\n=== Newton equilibrium (xdot=0) at p_load=%.3f pu ===\n', p.p_load0);

rhs_u = @(xx,uu) odefun_rhs_reduced(0, xx, p, @(t) uu, idx);
x_eq = run_newton_solver(x, rhs_u, p.p_load0, 80, 1e-10, 1e-10, true);

fprintf('Equilibrium done. ||F||=%.3e\n', norm(rhs_u(x_eq,p.p_load0),2));

%% Small-signal linearization
fprintf('\n=== Linearization around equilibrium ===\n');

rhs_eq = @(xx) odefun_rhs_reduced(0, xx, p, @(tt) p.p_load0, idx);
A      = numeric_jacobian(rhs_eq, x_eq);

du = 1e-6*(1+abs(p.p_load0));
b  = (rhs_u(x_eq, p.p_load0+du) - rhs_u(x_eq, p.p_load0)) / du;

yfun = @(xx,uu) output_ppcc(xx,p,uu,idx);
c    = numeric_jacobian(@(xx) yfun(xx,p.p_load0), x_eq);

print_eigs_participation(A, state_names, b, c);

%% Load-step response: 0.5 pu to 0.6 pu
t_step = 0.20;      % s
t_end  = 2.00;      % s

fprintf('\n=== Time-domain load step response ===\n');
fprintf('p_load: %.3f pu -> %.3f pu at t = %.3f s\n', p.p_load0, p.p_load1, t_step);

% Simulate only after the step.  Before t_step the system stays at the
% 0.5 pu equilibrium, so the pre-step trajectory is filled analytically.
u_fun = @(t) p.p_load0 + (t >= t_step).*(p.p_load1 - p.p_load0);

ode_opts = odeset('RelTol',1e-4, 'AbsTol',1e-6, ...
                  'MaxStep',2e-3, 'InitialStep',1e-5);
sol = ode15s(@(t,xx) odefun_rhs_reduced(t, xx, p, @(tt) p.p_load1, idx), ...
             [t_step t_end], x_eq, ode_opts);

t_pre  = linspace(0, t_step, 300).';
X_pre  = repmat(x_eq.', numel(t_pre), 1);
t_post = linspace(t_step, t_end, 3500).';
X_post = deval(sol, t_post).';

t = [t_pre; t_post(2:end)];
X = [X_pre; X_post(2:end,:)];
S = collect_signals_reduced(t, X, p, u_fun, idx);

p_pcc_eq0 = output_ppcc(x_eq, p, p.p_load0, idx);
dp_load   = S.p_load_cmd - p.p_load0;
dp_pcc    = S.p_pcc      - p_pcc_eq0;

figure('Name','Load step response: 0.5 pu to 0.6 pu');
plot(t, S.p_load_cmd, '--', 'LineWidth', 1.2); grid on; hold on;
plot(t, S.p_pcc, 'LineWidth', 1.2);
xline(t_step, '--', 'HandleVisibility','off');
xlabel('Time (s)');
ylabel('Power (pu)');
legend('p_{load}(t)', 'p_{pcc}(t)', 'Location','best');
title('Load step response: 0.5 pu \rightarrow 0.6 pu');
format_ieee_singlecol(gcf, 2, 4);

fprintf('\nStep-response summary:\n');
fprintf('  max|Delta p_load| = %.4g pu\n', max(abs(dp_load)));
fprintf('  max|Delta p_pcc | = %.4g pu\n', max(abs(dp_pcc)));
fprintf('  final p_pcc      = %.6f pu\n', S.p_pcc(end));


%% Local functions
function [x, nF] = run_newton_solver(x, rhs_u, p_load, maxIt, tolF, tolX, do_print)
    for it = 1:maxIt
        F  = rhs_u(x, p_load);
        nF = norm(F,2);
        if do_print
            fprintf('it=%2d  ||F||=%.3e\n', it, nF);
        end
        if nF < tolF, break; end

        J  = numeric_jacobian(@(xx) rhs_u(xx, p_load), x);
        dx = -J \ F;
        if norm(dx,2) < tolX*(1+norm(x,2)), break; end

        alpha = 1.0; c1 = 1e-4;
        while alpha > 1e-6
            xt = x + alpha*dx;
            Ft = rhs_u(xt, p_load);
            if norm(Ft,2) <= (1-c1*alpha)*nF
                x = xt; break;
            end
            alpha = 0.5*alpha;
        end
        if alpha <= 1e-6
            if do_print
                warning('Newton line-search stalled; applying damped step.');
            end
            x = x + 1e-3*dx;
        end
    end
end

function [xdot, alg] = odefun_rhs_reduced(t, x, p, pload_fun, idx)
    theta_pll = x(idx.theta_pll); eps_pll = x(idx.eps_pll); vq_pll_f = x(idx.vq_pll);
    id_afe = x(idx.id_afe); iq_afe = x(idx.iq_afe); xi_dc_afe = x(idx.xi_dc_afe);
    gamd_afe = x(idx.gamd_afe); gamq_afe = x(idx.gamq_afe); vdc_ups = x(idx.vdc_ups);
    iU_cv = x(idx.iU_cv); iV_cv = x(idx.iV_cv); vU_vsi = x(idx.vU_vsi); vV_vsi = x(idx.vV_vsi);
    xiU_vsi = x(idx.xiU_vsi); xiV_vsi = x(idx.xiV_vsi); gamU_vsi = x(idx.gamU_vsi);
    gamV_vsi = x(idx.gamV_vsi); v_psu = x(idx.v_psu); xi_psu = x(idx.xi_psu);
    v_eq = x(idx.v_eq); xi_eq = x(idx.xi_eq);

    p_load = pload_fun(t);
    g_load = p_load / (3*max(p.veq_ref,1e-9)^2);
    i_eq   = p.kp_v_eq*(p.veq_ref - v_eq) + p.ki_v_eq*xi_eq;
    i_psu  = safe_div(v_eq, v_psu) * i_eq;
    g_eq   = p.kp_v_psu*(p.vpsu_ref - v_psu) + p.ki_v_psu*xi_psu;

    iU_vsi = g_eq * vU_vsi;
    iV_vsi = g_eq * vV_vsi;
    vuv_sq = vU_vsi^2 + vV_vsi^2;
    psu_injection_term = ((g_eq - p.r_psu*g_eq^2) * vuv_sq) / (3*max(v_psu,1e-9));

    omega_vsi  = p.omega_vsi; vU_ref_vsi = p.vu_vsi_ref; vV_ref_vsi = 0;
    iU_ref_cv = p.kp_v_vsi*(vU_ref_vsi - vU_vsi) + p.ki_v_vsi*xiU_vsi - omega_vsi*p.c_vsi*vV_vsi;
    iV_ref_cv = p.kp_v_vsi*(vV_ref_vsi - vV_vsi) + p.ki_v_vsi*xiV_vsi + omega_vsi*p.c_vsi*vU_vsi;
    vU_ref_cv = p.kp_c_vsi*(iU_ref_cv - iU_cv) + p.ki_c_vsi*gamU_vsi - omega_vsi*p.l_vsi*iV_cv;
    vV_ref_cv = p.kp_c_vsi*(iV_ref_cv - iV_cv) + p.ki_c_vsi*gamV_vsi + omega_vsi*p.l_vsi*iU_cv;
    mU = safe_div(vU_ref_cv, vdc_ups); mV = safe_div(vV_ref_cv, vdc_ups);

    omega_pll = p.omega_s + p.kp_pll*vq_pll_f + p.ki_pll*eps_pll;
    s_ang = sin(theta_pll + pi/2); c_ang = cos(theta_pll + pi/2);

    ir_pcc =  s_ang*id_afe + c_ang*iq_afe; ii_pcc = -c_ang*id_afe + s_ang*iq_afe;
    vr_pcc = p.Vinf - p.Rinf*ir_pcc + p.Xinf*ii_pcc; vi_pcc = 0 - p.Rinf*ii_pcc - p.Xinf*ir_pcc;
    vd_pcc = s_ang*vr_pcc - c_ang*vi_pcc; vq_pcc = c_ang*vr_pcc + s_ang*vi_pcc;

    id_ref_afe = p.kp_dc_afe*(p.vdc_ups_ref - vdc_ups) + p.ki_dc_afe*xi_dc_afe; iq_ref_afe = 0;
    vd_ref_afe = p.kp_c_afe*(id_afe - id_ref_afe) + p.ki_c_afe*gamd_afe + omega_pll*p.l_afe*iq_afe;
    vq_ref_afe = p.kp_c_afe*(iq_afe - iq_ref_afe) + p.ki_c_afe*gamq_afe - omega_pll*p.l_afe*id_afe;

    md = safe_div(vd_ref_afe, vdc_ups); mq = safe_div(vq_ref_afe, vdc_ups);
    i_dc_in  = md*id_afe + mq*iq_afe; i_dc_out = mU*iU_cv + mV*iV_cv;

    xdot = zeros(size(x));
    xdot(idx.theta_pll) = p.omega_b*(omega_pll - p.omega_s);
    xdot(idx.eps_pll)   = vq_pll_f;
    xdot(idx.vq_pll)    = p.omega_lp*(vq_pcc - vq_pll_f);
    xdot(idx.id_afe)    = (p.omega_b/p.l_afe) * (vd_pcc - md*vdc_ups - p.r_afe*id_afe + omega_pll*p.l_afe*iq_afe);
    xdot(idx.iq_afe)    = (p.omega_b/p.l_afe) * (vq_pcc - mq*vdc_ups - p.r_afe*iq_afe - omega_pll*p.l_afe*id_afe);
    xdot(idx.xi_dc_afe) = (p.vdc_ups_ref - vdc_ups);
    xdot(idx.gamd_afe)  = (id_afe - id_ref_afe);
    xdot(idx.gamq_afe)  = (iq_afe - iq_ref_afe);
    xdot(idx.vdc_ups)   = (p.omega_b/p.c_dc) * (i_dc_in - i_dc_out);
    xdot(idx.iU_cv)     = (p.omega_b/p.l_vsi) * (mU*vdc_ups - vU_vsi - p.r_vsi*iU_cv + omega_vsi*p.l_vsi*iV_cv);
    xdot(idx.iV_cv)     = (p.omega_b/p.l_vsi) * (mV*vdc_ups - vV_vsi - p.r_vsi*iV_cv - omega_vsi*p.l_vsi*iU_cv);
    xdot(idx.vU_vsi)    = (p.omega_b/p.c_vsi) * (iU_cv - iU_vsi + omega_vsi*p.c_vsi*vV_vsi);
    xdot(idx.vV_vsi)    = (p.omega_b/p.c_vsi) * (iV_cv - iV_vsi - omega_vsi*p.c_vsi*vU_vsi);
    xdot(idx.xiU_vsi)   = (vU_ref_vsi - vU_vsi);
    xdot(idx.xiV_vsi)   = (vV_ref_vsi - vV_vsi);
    xdot(idx.gamU_vsi)  = (iU_ref_cv - iU_cv);
    xdot(idx.gamV_vsi)  = (iV_ref_cv - iV_cv);
    xdot(idx.v_psu)     = (p.omega_b/p.c_psu) * (psu_injection_term - i_psu);
    xdot(idx.xi_psu)    = (p.vpsu_ref - v_psu);
    xdot(idx.v_eq)      = (p.omega_b/p.c_eq) * (i_eq - g_load*v_eq);
    xdot(idx.xi_eq)     = (p.veq_ref - v_eq);

    alg.p_load = p_load; alg.g_load = g_load; alg.i_eq = i_eq; alg.i_psu = i_psu; alg.g_eq = g_eq;
    alg.vdc_ups = vdc_ups; alg.v_psu = v_psu; alg.v_eq = v_eq;
    alg.p_pcc = vr_pcc*ir_pcc + vi_pcc*ii_pcc;
    alg.p_dc_in = vdc_ups*(md*id_afe + mq*iq_afe);
    alg.p_dc_out = vdc_ups*(mU*iU_cv + mV*iV_cv);
end

function z = output_ppcc(x, p, p_load_scalar, idx)
    [~, alg] = odefun_rhs_reduced(0, x, p, @(t) p_load_scalar, idx);
    z = alg.p_pcc;
end

function S = collect_signals_reduced(t, X, p, u_fun, idx)
    nt = numel(t);

    fields = {'p_load_cmd','p_pcc','p_dc_in','p_dc_out','vdc_ups','v_psu','v_eq','g_eq','g_load','i_eq','i_psu'};
    for kf = 1:numel(fields), S.(fields{kf}) = zeros(nt,1); end

    for k = 1:nt
        [~, alg] = odefun_rhs_reduced(0, X(k,:).', p, @(tt) u_fun(t(k)), idx);

        S.p_load_cmd(k) = alg.p_load; S.p_pcc(k) = alg.p_pcc; S.p_dc_in(k) = alg.p_dc_in;
        S.p_dc_out(k) = alg.p_dc_out; S.vdc_ups(k) = alg.vdc_ups; S.v_psu(k) = alg.v_psu;
        S.v_eq(k) = alg.v_eq; S.g_eq(k) = alg.g_eq; S.g_load(k) = alg.g_load;
        S.i_eq(k) = alg.i_eq; S.i_psu(k) = alg.i_psu;
    end
end

function J = numeric_jacobian(fun, x0)
    f0 = fun(x0);
    n = numel(x0); m = numel(f0);
    J = zeros(m,n);
    for i = 1:n
        x1 = x0; hi = 1e-7*(1+abs(x0(i))); x1(i) = x1(i) + hi;
        J(:,i) = (fun(x1) - f0)/hi;
    end
end

function y = safe_div(a,b)
    if abs(b) < 1e-9, y = a/(signnz(b)*1e-9); else, y = a/b; end
end

function s = signnz(x)
    if x >= 0, s = 1; else, s = -1; end
end

function print_eigs_participation(A, state_names, b, c_out)
    fprintf('\nEigenvalues and participation factors (A-matrix)\n');

    % Right/left eigenvectors
    [V,D] = eig(A);
    lam = diag(D);
    [Wraw,Draw] = eig(A.');
    lamL = diag(Draw);

    % Match left eigenvectors to right eigenvalues
    W = zeros(size(Wraw));
    used = false(numel(lamL),1);
    for k = 1:numel(lam)
        cost = abs(lamL - lam(k));
        cost(used) = cost(used) + 1e3;
        [~,jj] = min(cost);
        used(jj) = true;
        W(:,k) = Wraw(:,jj);
    end
    W = W.';  % rows are left eigenvectors

    % Bi-orthonormalize so that W(k,:)*V(:,k)=1
    for k = 1:numel(lam)
        s = W(k,:)*V(:,k);
        if abs(s) > 1e-12
            W(k,:) = W(k,:)/s;
        end
    end

    % Participation factors
    P = V .* (W');
    Pabs = abs(P);
    Pnorm = Pabs ./ (sum(Pabs,1)+1e-15);

    % Residues
    R = zeros(numel(lam),1);
    for k = 1:numel(lam)
        R(k) = (c_out * V(:,k)) * (W(k,:) * b);
    end

    % Sort by real part descending
    [~,ks] = sort(real(lam),'descend');
    lam_s = lam(ks);
    Pnorm_s = Pnorm(:,ks);
    R_s = R(ks);
    Rabs_s = abs(R_s);
    sigma = real(lam_s);
    omg   = imag(lam_s);
    fhz   = abs(omg)/(2*pi);

    % Damping ratio
    zeta = nan(size(lam_s));
    for k = 1:numel(lam_s)
        den = sqrt(sigma(k)^2 + omg(k)^2);
        if den > 1e-12
            zeta(k) = -sigma(k)/den;
        end
    end

    % Color by log10(|R|)
    Rcolor = log10(Rabs_s + 1e-16);
    cmin = floor(min(Rcolor));
    cmax = ceil(max(Rcolor));
    if cmin == cmax
        cmin = cmin-1;
        cmax = cmax+1;
    end

    % Side-by-side eigenmap layout
    figure('Name','Eigenvalues in complex plane');
    % Set global colormap
    colormap(gcf, parula);

    % Left panel: global eigenmap
    % Position: [left bottom width height] -> Width is 0.50
    ax_main = axes('Position', [0.08 0.20 0.50 0.70]);
    scatter(ax_main, real(lam_s), imag(lam_s), 15, Rcolor, 'filled');
    grid(ax_main,'on'); hold(ax_main,'on');
    xline(ax_main, 0,'--k', 'LineWidth', 0.8, 'Alpha', 0.5);
    yline(ax_main, 0,'--k', 'LineWidth', 0.8, 'Alpha', 0.5);

    % Global Y limits and ticks (spaced by 2000)
    ylim(ax_main, [-6000, 6000]);
    yticks(ax_main, -6000:2000:6000);
    xlabel(ax_main,'$\Re(\lambda)$','Interpreter','latex');
    ylabel(ax_main,'$\Im(\lambda)$','Interpreter','latex');
    caxis(ax_main, [cmin cmax]);

    % Assign indices to each conjugate root
    lam_temp = lam_s;
    [~, sort_idx] = sort(abs(real(lam_temp)));
    lam_temp = lam_temp(sort_idx);
    mode_idx = 1;
    labels = {};
    coords  = complex([]);   % store eigenvalues to label (complex)
    while ~isempty(lam_temp)
        curr = lam_temp(1);
        if abs(imag(curr)) > 1e-5
            match_idx = find(abs(lam_temp - conj(curr)) < 1e-5, 1);
            if imag(curr) > 0
                upper_root = curr;
                lower_root = conj(curr);
            else
                upper_root = conj(curr);
                lower_root = curr;
            end
            coords(end+1) = upper_root;
            labels{end+1} = sprintf(' %d', mode_idx);
            coords(end+1) = lower_root;
            labels{end+1} = sprintf(' %d', mode_idx + 1);
            mode_idx = mode_idx + 2;
            remove_idx = 1;
            if ~isempty(match_idx)
                remove_idx = [1, match_idx];
            end
            lam_temp(remove_idx) = [];
        else
            coords(end+1) = curr;
            labels{end+1} = sprintf(' %d', mode_idx);
            mode_idx = mode_idx + 1;
            lam_temp(1) = [];
        end
    end

    % Right panel: zoomed-in eigenmap
    % Make width exactly half of the left panel (0.25 vs 0.50)
    ax_in = axes('Position', [0.66 0.20 0.25 0.70]);
    box(ax_in,'on'); grid(ax_in,'on'); hold(ax_in,'on');
    scatter(ax_in, real(lam_s), imag(lam_s), 15, Rcolor, 'filled');
    xline(ax_in, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.5);
    yline(ax_in, 0, '--k', 'LineWidth', 0.8, 'Alpha', 0.5);
    xlim(ax_in, [-400, 20]);
    roi = (real(lam_s) >= -500) & (real(lam_s) <= 20);
    if any(roi)
        yim = imag(lam_s(roi));
        yr  = max(yim) - min(yim);
        ypad = 0.10 * max(1, yr);
        ylim(ax_in, [min(yim)-ypad, max(yim)+ypad]);
    else
        ylim(ax_in, [-1000, 1000]);
    end
    xlabel(ax_in,'$\Re(\lambda)$','Interpreter','latex');
    ylabel(ax_in,'$\Im(\lambda)$','Interpreter','latex');
    caxis(ax_in, [cmin cmax]);

    % Label roots on the main and zoomed axes
    for m = 1:numel(coords)
        r = real(coords(m));
        i = imag(coords(m));

        % Extract integer mode number from the label string
        mnum = sscanf(labels{m}, '%d');

        % Main labels: skip labeling the first 13 modes to avoid clutter
        if mnum > 13
            if i >= -1e-5
                text(ax_main, r, i, labels{m}, 'Interpreter','latex', ...
                    'FontSize',7, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
            else
                text(ax_main, r, i, labels{m}, 'Interpreter','latex', ...
                    'FontSize',7, 'VerticalAlignment','top', 'HorizontalAlignment','left');
            end
        end

        % Inset labels only within ROI
        if (r >= -500) && (r <= 20)
            if i >= -1e-5
                text(ax_in, r, i, labels{m}, 'Interpreter','latex', ...
                    'FontSize',7, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
            else
                text(ax_in, r, i, labels{m}, 'Interpreter','latex', ...
                    'FontSize',7, 'VerticalAlignment','top', 'HorizontalAlignment','left');
            end
        end
    end

    % Add zoom rectangle and connecting lines
    % Force update to make sure limits and positions are current before drawing lines
    drawnow;

    xin_lim = xlim(ax_in);
    yin_lim = ylim(ax_in);

    % Draw rectangle on main axis matching the limits of zoomed axis
    rect_width = xin_lim(2) - xin_lim(1);
    rect_height = yin_lim(2) - yin_lim(1);
    rectangle(ax_main, 'Position', [xin_lim(1), yin_lim(1), rect_width, rect_height], ...
              'EdgeColor', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.2);

    % Map data coordinates to normalized figure coordinates for the connecting lines
    pos_main = get(ax_main, 'Position');
    xlim_main = xlim(ax_main);
    ylim_main = ylim(ax_main);

    pos_in = get(ax_in, 'Position');

    % Right edge of the rectangle in normalized units
    norm_rect_x = pos_main(1) + (xin_lim(2) - xlim_main(1)) / (xlim_main(2) - xlim_main(1)) * pos_main(3);
    norm_rect_y_top = pos_main(2) + (yin_lim(2) - ylim_main(1)) / (ylim_main(2) - ylim_main(1)) * pos_main(4);
    norm_rect_y_bot = pos_main(2) + (yin_lim(1) - ylim_main(1)) / (ylim_main(2) - ylim_main(1)) * pos_main(4);

    % Left edge of the zoomed-in axis in normalized units
    norm_ax2_x = pos_in(1);
    norm_ax2_y_top = pos_in(2) + pos_in(4);
    norm_ax2_y_bot = pos_in(2);

    % Draw connecting lines
    annotation('line', [norm_rect_x, norm_ax2_x], [norm_rect_y_top, norm_ax2_y_top], ...
               'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1);
    annotation('line', [norm_rect_x, norm_ax2_x], [norm_rect_y_bot, norm_ax2_y_bot], ...
               'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1);

    % Shared colorbar
    cb = colorbar(ax_in);
    cb.Position = [0.94 0.20 0.015 0.70]; % Explicitly place at the far right
    cb.Ticks = cmin:cmax;
    cb.TickLabels = arrayfun(@(e) sprintf('$10^{%d}$', e), cb.Ticks, 'UniformOutput', false);
    cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$|R_k|$';
    cb.Label.Interpreter = 'latex';

    % Apply (2, 9) sizing globally for this figure
    format_ieee_singlecol(gcf, 2, 10);

    % Print tables
    fprintf('\n--- Eigenvalues (sorted by Re) ---\n');
    fprintf('%4s  %14s  %14s  %10s  %10s  %10s  %14s\n', ...
        'k', 'Re(lambda)', 'Im(lambda)', 'f(Hz)', 'zeta', 'stable?', '|Residue|');
    for kk = 1:numel(lam_s)
        st = "stable";
        if sigma(kk) > 1e-8
            st = "UNSTABLE";
        end
        fprintf('%4d  %14.6e  %14.6e  %10.3f  %10.3f  %10s  %14.6e\n', ...
            kk, sigma(kk), omg(kk), fhz(kk), zeta(kk), st, Rabs_s(kk));
    end
    topN = 6;
    fprintf('\n--- Dominant participation per mode (Top-%d states) ---\n', topN);
    for kk = 1:numel(lam_s)
        [vals,ii] = sort(Pnorm_s(:,kk),'descend');
        fprintf('\nMode %d: lambda = %.6e %+.6ej, f=%.3f Hz, |R|=%.3e\n', ...
            kk, sigma(kk), omg(kk), fhz(kk), Rabs_s(kk));
        for j = 1:min(topN,numel(ii))
            fprintf('  %2d) %-12s  PF=%.4f\n', j, state_names{ii(j)}, vals(j));
        end
    end

    eig_tbl = table((1:numel(lam_s))', sigma, omg, fhz, zeta, Rabs_s, real(R_s), imag(R_s), ...
        'VariableNames', {'mode_sorted','Re','Im','freq_Hz','zeta','residue_abs','residue_re','residue_im'});
    writetable(eig_tbl, 'eigvals_sorted.csv');

    Ptab = array2table(Pnorm_s, 'VariableNames', cellstr(compose('mode%02d', 1:numel(lam_s))));
    Ptab.State = string(state_names(:));
    writetable(movevars(Ptab, 'State', 'Before', 1), 'participation_abs_norm.csv');

    fprintf('\nSaved: eigvals_sorted.csv, participation_abs_norm.csv\n');
end
function [kp, ki] = pi_bw_tuning(loopType, fbw, zeta, param1, omega_b, param2)
    if nargin < 6, param2 = 0; end
    wn = 2*pi*fbw;
    switch lower(loopType)
        case 'voltage', kp = 2*zeta*wn * (param1/omega_b); ki = (wn^2) * (param1/omega_b);
        case 'current', kp = 2*zeta*wn * (param1/omega_b) - param2; ki = (wn^2) * (param1/omega_b);
    end
end

function [kp, ki] = pll_bw_tuning_simple(fbw, zeta)
    wn = 2*pi*fbw; kp = (2*zeta*wn) / (2*pi*60); ki = (wn^2) / (2*pi*60);
end

function format_ieee_singlecol(fig, h_in, w_in)
    if nargin < 2 || isempty(h_in), h_in = 2.0; end
    if nargin < 3 || isempty(w_in), w_in = 4.0; end

    set(fig, 'Color', 'none', 'Units', 'inches', 'Position', [1 1 w_in h_in], 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in], 'PaperSize', [w_in h_in]);
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k), 'FontName','Times New Roman', 'FontSize',8, 'LineWidth',0.8, 'Box','on', 'TickDir','out');
        grid(ax(k),'on'); ax(k).GridAlpha = 0.18; ax(k).MinorGridAlpha = 0.12;
    end
    ln = findall(fig, 'Type', 'line');
    for k = 1:numel(ln), if ln(k).LineWidth < 1.0, ln(k).LineWidth = 1.0; end, end
    tx = findall(fig, '-property', 'FontSize');
    for k = 1:numel(tx)
        try
            if isprop(tx(k),'FontName'), tx(k).FontName = 'Times New Roman'; end
            if isprop(tx(k),'FontSize')
                if isa(tx(k), 'matlab.graphics.illustration.Legend'), tx(k).FontSize = 7; else, tx(k).FontSize = 8; end
            end
        catch
        end
    end
end
