%% GPU Load Response and FFT Analysis
% Simulates the single-data-center infinite-bus model under a measured GPU
% power trace using implicit trapezoidal integration, then extracts the
% dominant load and PCC spectral components.
%
% Required input:
%   GPU_data.csv with columns t_seconds and P_gpu_W.
clear; clc; close all;

% 1. Parameters
wb = 2*pi*60;
omega_s = 1.0;
omega_lp = 2*pi*100;

Vinf = 1.0; Rinf = 0.02; Xinf = 0.19;
l_afe = 0.05; r_afe = 0.003; c_dc = 2.0; vdc_ups_ref = 1.0;
l_vsi = 0.05; r_vsi = 0.003; c_vsi = 0.2; vu_vsi_ref = 1.0;
l_psu = 0.05; c_psu = 2.0; r_psu = 0.005; v_psu_ref = 1.0;
l_eq  = 0.05; c_eq  = 0.2; v_eq_ref = 0.5;

% 2. PI gains
[kp_pll, ki_pll]       = pll_bw_tuning(20.0, 0.707);
[kp_dc_afe, ki_dc_afe] = pi_bw_tuning('v', 5.0,   1.0,   c_dc,  wb, 0);
[kp_c_afe,  ki_c_afe]  = pi_bw_tuning('c', 200.0, 0.707, l_afe, wb, r_afe);
[kp_v_vsi,  ki_v_vsi]  = pi_bw_tuning('v', 100.0, 1.0,   c_vsi, wb, 0);
[kp_c_vsi,  ki_c_vsi]  = pi_bw_tuning('c', 400.0, 1.0,   l_vsi, wb, r_vsi);
[kp_v_psu,  ki_v_psu]  = pi_bw_tuning('v', 10.0,  1.0,   c_psu, wb, 0);
[kp_v_eq,   ki_v_eq]   = pi_bw_tuning('v', 100.0, 1.0,   c_eq,  wb, 0);

params = {wb, omega_s, omega_lp, Vinf, Rinf, Xinf, ...
          l_afe, r_afe, c_dc, vdc_ups_ref, ...
          l_vsi, r_vsi, c_vsi, vu_vsi_ref, ...
          l_psu, c_psu, r_psu, v_psu_ref, ...
          l_eq, c_eq, v_eq_ref, ...
          kp_pll, ki_pll, kp_dc_afe, ki_dc_afe, kp_c_afe, ki_c_afe, ...
          kp_v_vsi, ki_v_vsi, kp_c_vsi, ki_c_vsi, ...
          kp_v_psu, ki_v_psu, kp_v_eq, ki_v_eq};

% 3. Load real GPU trace and build workload on [100, 200] s
fprintf('Loading GPU_data.csv ...\n');

gpu_tbl = readtable('GPU_data.csv');

t_gpu_raw = gpu_tbl.t_seconds;
P_gpu_raw = gpu_tbl.P_gpu_W;

valid = isfinite(t_gpu_raw) & isfinite(P_gpu_raw);
t_gpu_raw = t_gpu_raw(valid);
P_gpu_raw = P_gpu_raw(valid);

if numel(t_gpu_raw) < 2
    error('GPU_data.csv does not contain enough valid samples.');
end

t_gpu_raw = t_gpu_raw - t_gpu_raw(1);
[t_gpu_raw, uniq_idx] = unique(t_gpu_raw, 'stable');
P_gpu_raw = P_gpu_raw(uniq_idx);

dt_load = median(diff(t_gpu_raw));
if dt_load <= 0
    error('Invalid time stamps in GPU_data.csv.');
end

t_start = 1030.0;
t_end   = 1130.0;
t_load  = (t_start:dt_load:t_end).';

% Linear mapping: min -> 0.3 pu, max -> 0.8 pu
Pmin = min(P_gpu_raw);
Pmax = max(P_gpu_raw);

if abs(Pmax - Pmin) < 1e-12
    error('GPU power trace is nearly constant; cannot map min/max to different pu values.');
end

P_gpu_pu_raw = 0.3 + (0.8 - 0.3) * (P_gpu_raw - Pmin) / (Pmax - Pmin);

% Build [100, 200] s trace by repeating the original trace if needed
P_load_samples = tile_signal_to_time(t_gpu_raw, P_gpu_pu_raw, t_load);

fprintf('GPU trace loaded.\n');
fprintf('Original GPU power (W): min = %.3f, mean = %.3f, max = %.3f\n', ...
    min(P_gpu_raw), mean(P_gpu_raw), max(P_gpu_raw));
fprintf('Mapped per-unit load on [%.1f, %.1f] s: min = %.3f, mean = %.3f, max = %.3f\n', ...
    t_start, t_end, min(P_load_samples), mean(P_load_samples), max(P_load_samples));
fprintf('Load sample interval = %.4f s\n', dt_load);

% 4. Newton initialization at first p_load sample of the window
p_load_init = P_load_samples(1);

fprintf('Running Newton initialization at first p_load sample ...\n');
fprintf('Initial p_load = %.4f pu\n', p_load_init);

g_load0 = p_load_init / (3*v_eq_ref^2);

X_guess = zeros(21,1);
X_guess(4)  = 0.5;
X_guess(6)  = 0.5 / ki_dc_afe;
X_guess(9)  = 1.0;
X_guess(10) = 0.5;
X_guess(12) = 1.0;
X_guess(14) = 0.5 / ki_v_vsi;

i_eq0  = g_load0 * v_eq_ref;
i_psu0 = (v_eq_ref / v_psu_ref) * i_eq0;

disc0 = 1 - 12*r_psu*i_psu0;
disc0 = max(disc0, 1e-8);
g_eq0 = (1 - sqrt(disc0)) / (2*r_psu);

X_guess(18) = 1.0;
X_guess(19) = g_eq0 / max(ki_v_psu, 1e-9);
X_guess(20) = 0.5;
X_guess(21) = i_eq0 / max(ki_v_eq, 1e-9);

[X, init_info] = newton_initialize_state(X_guess, p_load_init, params);
fprintf('Initialization finished: iter = %d, ||xdot(x0)||_inf = %.3e\n', ...
    init_info.iter, init_info.res_inf);

if ~init_info.converged
    warning('Newton initialization did not fully converge.');
end

% 5. Implicit trapezoidal time-domain simulation
n_sub = 40;

fprintf('Running implicit trapezoidal simulation from %.1f s to %.1f s ...\n', ...
    t_start, t_end);
fprintf('Using n_sub = %d substeps per load interval.\n', n_sub);

n_seg = numel(t_load) - 1;
N_out = n_seg * n_sub + 1;

t_arr      = zeros(N_out,1);
P_load_res = zeros(N_out,1);
P_pcc_res  = zeros(N_out,1);

idx_out = 1;
t_arr(idx_out) = t_load(1);
P_load_res(idx_out) = P_load_samples(1);
[~, P_pcc0] = get_sys_B_derivatives(X, P_load_samples(1), params);
P_pcc_res(idx_out) = P_pcc0;

for k = 1:n_seg
    t0 = t_load(k);
    t1 = t_load(k+1);
    p_load_k = P_load_samples(k);

    h = (t1 - t0) / n_sub;

    for j = 1:n_sub
        X = implicit_trap_step(X, p_load_k, h, params);

        idx_out = idx_out + 1;
        t_arr(idx_out) = t0 + j*h;
        P_load_res(idx_out) = p_load_k;

        [~, P_pcc_now] = get_sys_B_derivatives(X, p_load_k, params);
        P_pcc_res(idx_out) = P_pcc_now;
    end
end

% 6. Time-domain plot
t_plot = t_arr - t_arr(1);   % make time start from 0

fig1 = figure;
plot(t_plot, P_load_res, 'r', 'LineWidth', 0.5, 'DisplayName', '$p_{\mathrm{load}}$'); hold on;
plot(t_plot, P_pcc_res,  'b', 'LineWidth', 0.5, 'DisplayName', '$p_{\mathrm{pcc}}$');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Power (p.u.)', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'northwest');
grid on;
xlim([0, t_plot(end)]);
ylim([0, 1.0]);
format_fig(fig1, 2.0, 5);

% 7. FFT spectrum extraction with selectable time window
fft_t_start = t_start+50.0;   % FFT window start time (s)
fft_t_end   = t_end;   % FFT window end time   (s)

idx_fft = (t_arr >= fft_t_start) & (t_arr <= fft_t_end);

if nnz(idx_fft) < 2
    error('FFT time window is too short or outside t_arr range.');
end

t_fft_use      = t_arr(idx_fft);
P_load_fft_use = P_load_res(idx_fft);
P_pcc_fft_use  = P_pcc_res(idx_fft);

dt_fft = median(diff(t_fft_use));
Fs = 1/dt_fft;

[f_load, A_load] = extract_single_sided_spectrum(P_load_fft_use, Fs);
[f_pcc,  A_pcc ] = extract_single_sided_spectrum(P_pcc_fft_use, Fs);

fmin = 0.01;
fmax = 10.0;
idx_load = (f_load >= fmin) & (f_load <= fmax);
idx_pcc  = (f_pcc  >= fmin) & (f_pcc  <= fmax);

[fpk_load, Apk_load] = find_dominant_peak(f_load(idx_load), A_load(idx_load));
[fpk_pcc,  Apk_pcc ] = find_dominant_peak(f_pcc(idx_pcc),   A_pcc(idx_pcc));

fprintf('\n=== FFT Time Window ===\n');
fprintf('Window: [%.2f s, %.2f s]\n', fft_t_start, fft_t_end);

fprintf('\n=== Dominant Spectral Components ===\n');
fprintf('Load: f_peak = %.4f Hz, amplitude = %.6e\n', fpk_load, Apk_load);
fprintf('PCC : f_peak = %.4f Hz, amplitude = %.6e\n', fpk_pcc,  Apk_pcc);

% 8. Spectrum plot
fig2 = figure;
plot(f_load(idx_load), A_load(idx_load), 'r', 'LineWidth', 0.5, ...
    'DisplayName', '$p_{\mathrm{load}}$'); hold on;
plot(f_pcc(idx_pcc), A_pcc(idx_pcc), 'b', 'LineWidth', 0.5, ...
    'DisplayName', '$p_{\mathrm{pcc}}$');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel('Amplitude (p.u.)', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'best');
grid on;
xlim([fmin, fmax]);
format_fig(fig2, 2.0, 5);


% Local functions


function P_target = tile_signal_to_time(t_base, P_base, t_target)
    T_base = t_base(end);
    if T_base <= 0
        error('Base signal duration must be positive.');
    end

    dt_base = median(diff(t_base));
    if dt_base <= 0
        error('Base signal sample interval must be positive.');
    end

    T_period = T_base + dt_base;
    P_target = zeros(size(t_target));

    for k = 1:numel(t_target)
        tau = mod(t_target(k), T_period);
        if tau > t_base(end)
            tau = t_base(end);
        end
        P_target(k) = interp1(t_base, P_base, tau, 'previous', 'extrap');
    end
end

function Xnext = implicit_trap_step(Xnow, p_load, h, params)
    % Implicit trapezoidal:
    % X_{n+1} = X_n + h/2 * (f(X_n) + f(X_{n+1}))
    %
    % Residual:
    % R(Xnext) = Xnext - Xnow - h/2*(fnow + f(Xnext)) = 0

    newton_tol = 1e-8;
    newton_max_iter = 20;
    fd_eps = 1e-7;
    reg = 1e-9;

    [fnow, ~] = get_sys_B_derivatives(Xnow, p_load, params);

    % Predictor: explicit Euler
    Xnext = Xnow + h * fnow;

    % Clamp key positive states
    Xnext(9)  = max(Xnext(9),  1e-4);
    Xnext(18) = max(Xnext(18), 1e-4);
    Xnext(20) = max(Xnext(20), 1e-4);

    for iter = 1:newton_max_iter
        R = trap_residual(Xnext, Xnow, fnow, p_load, h, params);
        if norm(R, inf) < newton_tol
            return;
        end

        J = numerical_jacobian(@(x) trap_residual(x, Xnow, fnow, p_load, h, params), ...
                               Xnext, R, fd_eps);

        dX = -(J + reg*eye(numel(Xnext))) \ R;

        alpha = 1.0;
        accepted = false;
        res0 = norm(R, inf);

        for ls = 1:10
            Xtry = Xnext + alpha*dX;
            Xtry(9)  = max(Xtry(9),  1e-4);
            Xtry(18) = max(Xtry(18), 1e-4);
            Xtry(20) = max(Xtry(20), 1e-4);

            Rtry = trap_residual(Xtry, Xnow, fnow, p_load, h, params);
            if norm(Rtry, inf) < res0
                Xnext = Xtry;
                accepted = true;
                break;
            end
            alpha = 0.5 * alpha;
        end

        if ~accepted
            Xnext = Xnext + 0.1*dX;
            Xnext(9)  = max(Xnext(9),  1e-4);
            Xnext(18) = max(Xnext(18), 1e-4);
            Xnext(20) = max(Xnext(20), 1e-4);
        end
    end

    warning('implicit_trap_step: Newton did not fully converge at one time step.');
end

function R = trap_residual(Xnext, Xnow, fnow, p_load, h, params)
    [fnext, ~] = get_sys_B_derivatives(Xnext, p_load, params);
    R = Xnext - Xnow - 0.5*h*(fnow + fnext);
end

function [dX, P_pcc] = get_sys_B_derivatives(X, p_load, params)
    [wb, omega_s, omega_lp, Vinf, Rinf, Xinf, ...
     l_afe, r_afe, c_dc, vdc_ups_ref, ...
     l_vsi, r_vsi, c_vsi, vu_vsi_ref, ...
     l_psu, c_psu, r_psu, v_psu_ref, ...
     l_eq, c_eq, v_eq_ref, ...
     kp_pll, ki_pll, kp_dc_afe, ki_dc_afe, kp_c_afe, ki_c_afe, ...
     kp_v_vsi, ki_v_vsi, kp_c_vsi, ki_c_vsi, ...
     kp_v_psu, ki_v_psu, kp_v_eq, ki_v_eq] = params{:};

    dX = zeros(21,1);
    xs = X(1:17);

    g_load = p_load / (3*v_eq_ref^2);
    v_psu  = X(18);
    xi_psu = X(19);
    v_eq   = X(20);
    xi_eq  = X(21);

    i_eq_ref = kp_v_eq*(v_eq_ref - v_eq) + ki_v_eq*xi_eq;
    dX(20) = (wb/c_eq) * (i_eq_ref - g_load*v_eq);
    dX(21) = v_eq_ref - v_eq;

    i_psu_ld = (v_eq / max(0.1, v_psu)) * i_eq_ref;
    g_eq = kp_v_psu*(v_psu_ref - v_psu) + ki_v_psu*xi_psu;
    v_vsi_sq = X(12)^2 + X(13)^2;

    dX(18) = (wb/c_psu) * ((g_eq - r_psu*g_eq^2) * v_vsi_sq / (3*max(v_psu,1e-9)) - i_psu_ld);
    dX(19) = v_psu_ref - v_psu;

    iU_ld = g_eq * X(12);
    iV_ld = g_eq * X(13);

    iU_ref = kp_v_vsi*(vu_vsi_ref - xs(12)) + ki_v_vsi*xs(14) - omega_s*c_vsi*xs(13);
    iV_ref = kp_v_vsi*(0 - xs(13))          + ki_v_vsi*xs(15) + omega_s*c_vsi*xs(12);

    vU_ref = kp_c_vsi*(iU_ref - xs(10)) + ki_c_vsi*xs(16) - omega_s*l_vsi*xs(11);
    vV_ref = kp_c_vsi*(iV_ref - xs(11)) + ki_c_vsi*xs(17) + omega_s*l_vsi*xs(10);

    mU = vU_ref / max(xs(9), 1e-9);
    mV = vV_ref / max(xs(9), 1e-9);

    omega_pll = omega_s + kp_pll*xs(3) + ki_pll*xs(2);
    s = sin(xs(1) + pi/2);
    c = cos(xs(1) + pi/2);

    ir =  s*xs(4) + c*xs(5);
    ii = -c*xs(4) + s*xs(5);

    vr = Vinf - Rinf*ir + Xinf*ii;
    vi = -Rinf*ii - Xinf*ir;

    vd = s*vr - c*vi;
    vq = c*vr + s*vi;

    P_pcc = vd*xs(4) + vq*xs(5);

    id_ref = kp_dc_afe*(vdc_ups_ref - xs(9)) + ki_dc_afe*xs(6);
    vd_ref = kp_c_afe*(xs(4) - id_ref) + ki_c_afe*xs(7) + omega_pll*l_afe*xs(5);
    vq_ref = kp_c_afe*xs(5)            + ki_c_afe*xs(8) - omega_pll*l_afe*xs(4);

    md = vd_ref / max(xs(9), 1e-9);
    mq = vq_ref / max(xs(9), 1e-9);

    dX(1)  = wb*(omega_pll - omega_s);
    dX(2)  = xs(3);
    dX(3)  = omega_lp*(vq - xs(3));
    dX(4)  = (wb/l_afe)*(vd - md*xs(9) - r_afe*xs(4) + omega_pll*l_afe*xs(5));
    dX(5)  = (wb/l_afe)*(vq - mq*xs(9) - r_afe*xs(5) - omega_pll*l_afe*xs(4));
    dX(6)  = vdc_ups_ref - xs(9);
    dX(7)  = xs(4) - id_ref;
    dX(8)  = xs(5);
    dX(9)  = (wb/c_dc)*((md*xs(4) + mq*xs(5)) - (mU*xs(10) + mV*xs(11)));
    dX(10) = (wb/l_vsi)*(mU*xs(9) - xs(12) - r_vsi*xs(10) + omega_s*l_vsi*xs(11));
    dX(11) = (wb/l_vsi)*(mV*xs(9) - xs(13) - r_vsi*xs(11) - omega_s*l_vsi*xs(10));
    dX(12) = (wb/c_vsi)*(xs(10) - iU_ld + omega_s*c_vsi*xs(13));
    dX(13) = (wb/c_vsi)*(xs(11) - iV_ld - omega_s*c_vsi*xs(12));
    dX(14) = vu_vsi_ref - xs(12);
    dX(15) = -xs(13);
    dX(16) = iU_ref - xs(10);
    dX(17) = iV_ref - xs(11);
end

function [X, info] = newton_initialize_state(X0, p_load, params)
    X = X0;
    tol = 1e-9;
    max_iter = 30;
    fd_eps = 1e-7;
    reg = 1e-8;

    for iter = 1:max_iter
        F = steady_state_residual(X, p_load, params);
        res_inf = norm(F, inf);

        if res_inf < tol
            info.iter = iter - 1;
            info.res_inf = res_inf;
            info.converged = true;
            return;
        end

        J = numerical_jacobian(@(x) steady_state_residual(x, p_load, params), X, F, fd_eps);
        dx = -(J + reg*eye(numel(X))) \ F;

        alpha = 1.0;
        accepted = false;
        for ls = 1:12
            X_try = X + alpha*dx;
            X_try(9)  = max(X_try(9),  0.05);
            X_try(18) = max(X_try(18), 0.05);
            X_try(20) = max(X_try(20), 0.05);

            F_try = steady_state_residual(X_try, p_load, params);
            if norm(F_try, inf) < res_inf
                X = X_try;
                accepted = true;
                break;
            end
            alpha = 0.5 * alpha;
        end

        if ~accepted
            X = X + 0.1*dx;
            X(9)  = max(X(9),  0.05);
            X(18) = max(X(18), 0.05);
            X(20) = max(X(20), 0.05);
        end
    end

    F = steady_state_residual(X, p_load, params);
    info.iter = max_iter;
    info.res_inf = norm(F, inf);
    info.converged = false;
end

function F = steady_state_residual(X, p_load, params)
    [F, ~] = get_sys_B_derivatives(X, p_load, params);
end

function J = numerical_jacobian(fun, x, Fx, eps_fd)
    n = numel(x);
    J = zeros(n, n);
    for k = 1:n
        xk = x;
        hk = eps_fd * max(1, abs(x(k)));
        xk(k) = xk(k) + hk;
        Fk = fun(xk);
        J(:, k) = (Fk - Fx) / hk;
    end
end

function [kp, ki] = pi_bw_tuning(type, fbw, zeta, p1, wb, p2)
    wn = 2*pi*fbw;
    if type == 'v'
        kp = 2*zeta*wn*(p1/wb);
    else
        kp = 2*zeta*wn*(p1/wb) - p2;
    end
    ki = wn^2*(p1/wb);
end

function [kp, ki] = pll_bw_tuning(fbw, zeta)
    wn = 2*pi*fbw;
    kp = 2*zeta*wn/(2*pi*60);
    ki = wn^2/(2*pi*60);
end

function [f, A] = extract_single_sided_spectrum(x, Fs)
    x = x(:) - mean(x);
    N = numel(x);
    w = hann(N);
    X = fft(x .* w);
    N_half = floor(N/2);
    f = (0:N_half)' * Fs / N;
    A = abs(X(1:N_half+1)) / (N*mean(w));
    if N > 2
        A(2:end-1) = 2*A(2:end-1);
    end
end

function [fpk, Apk] = find_dominant_peak(f, A)
    idx = f > 0.01;
    f = f(idx);
    A = A(idx);
    if isempty(f)
        fpk = NaN;
        Apk = NaN;
        return;
    end
    [Apk, k] = max(A);
    fpk = f(k);
end

function format_fig(fig, h, w)
    set(fig, 'Units', 'inches', 'Position', [1 1 w h]);
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k), 'FontSize', 9, 'FontName', 'Times New Roman', 'LineWidth', 0.8);
    end
    lgd = findall(fig, 'Type', 'legend');
    for k = 1:numel(lgd)
        set(lgd(k), 'FontSize', 9, 'FontName', 'Times New Roman');
    end
end
