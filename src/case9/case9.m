%% IEEE 9-Bus Data-Center Stability Study
% Supplementary MATLAB script for the IEEE 9-bus case.
% The workflow solves the power flow, performs Kron reduction, initializes
% SM/GFM/GFL/data-center port models, evaluates small-signal modes and POA
% gain, and runs a GPU-load time-domain simulation.
%
% Required input for the time-domain portion:
%   GPU_data.csv with columns t_seconds and P_gpu_W.
clear; clc; close all;

%% 1. System parameters and Ybus
baseMVA = 100; fb = 60; scale_lev = 0.5;

% Bus data: bus 8 is the dynamic data-center load; other loads are constant PQ.
% Format: [bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin]
bus = [
  1 3 0   0   0 0 1 1.040 0 345 1 1.1 0.9;
  2 2 0   0   0 0 1 1.025 0 345 1 1.1 0.9;
  3 1 0   0   0 0 1 1.025 0 345 1 1.1 0.9;
  4 1 0   0   0 0 1 1.000 0 345 1 1.1 0.9;
  5 1 125 50  0 0 1 1.000 0 345 1 1.1 0.9; % Constant-impedance load
  6 1 90  30  0 0 1 1.000 0 345 1 1.1 0.9; % Constant-impedance load
  7 1 0   0   0 0 1 1.000 0 345 1 1.1 0.9;
  8 1 100 0   0 0 1 1.000 0 345 1 1.1 0.9; % Dynamic data-center bus
  9 1 0   0   0 0 1 1.000 0 345 1 1.1 0.9;
];
bus(:,3:4) = bus(:,3:4) * scale_lev;

% Generator power-flow data
gen_pf = [1 0    0   300 -300 1.04  baseMVA 1 250 10;
          2 163  0   300 -300 1.025 baseMVA 1 300 10];
gen_pf(:,2:3) = gen_pf(:,2:3) * scale_lev;

% Branch data
branch = [
  1 4 0.0000 0.0576 0.0000 250 250 250 0 0 1 -360 360;
  4 5 0.0170 0.0920 0.1580 250 250 250 0 0 1 -360 360;
  5 6 0.0390 0.1700 0.3580 150 150 150 0 0 1 -360 360;
  3 6 0.0000 0.0586 0.0000 300 300 300 0 0 1 -360 360;
  6 7 0.0119 0.1008 0.2090 150 150 150 0 0 1 -360 360;
  7 8 0.0085 0.0720 0.1490 250 250 250 0 0 1 -360 360;
  8 2 0.0000 0.0625 0.0000 250 250 250 0 0 1 -360 360;
  8 9 0.0320 0.1610 0.3060 250 250 250 0 0 1 -360 360;
  9 4 0.0100 0.0850 0.1760 250 250 250 0 0 1 -360 360;
];

% Construct Ybus matrix
nb = size(bus,1);
Ybus = zeros(nb,nb);
for k = 1:size(branch,1)
    f = branch(k,1);
    t = branch(k,2);
    y = 1 / complex(branch(k,3), branch(k,4));
    ysh = 1j * branch(k,5) / 2;
    Ybus(f,f) = Ybus(f,f) + y + ysh;
    Ybus(t,t) = Ybus(t,t) + y + ysh;
    Ybus(f,t) = Ybus(f,t) - y;
    Ybus(t,f) = Ybus(t,f) - y;
end
Ybus = Ybus + diag(complex(bus(:,5), bus(:,6)) / baseMVA);

%% 2. Power flow (Newton-Raphson)
Pd = bus(:,3) / baseMVA;
Qd = bus(:,4) / baseMVA;
Pg = zeros(nb,1);
Qg0 = zeros(nb,1);
Pg(gen_pf(:,1)) = gen_pf(:,2) / baseMVA;
Qg0(gen_pf(:,1)) = gen_pf(:,3) / baseMVA;

% GFL inverter injection at Bus 3
S3_pu = (85.00 + 1j * -10.94) * scale_lev / baseMVA;
Pg(3) = real(S3_pu);
Qg0(3) = imag(S3_pu);

V = bus(:,8) .* exp(1j * bus(:,9) * pi / 180);
pvpq = find(bus(:,2) ~= 3);
pq = find(bus(:,2) == 1);
pf_fun = @(x) pf_mismatch(x, Ybus, V, pvpq, pq, Pg - Pd, Qg0 - Qd);
opts_pf = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12);
x_pf = fsolve(pf_fun, [angle(V(pvpq)); abs(V(pq))], opts_pf);

V(pvpq) = abs(V(pvpq)) .* exp(1j * x_pf(1:numel(pvpq)));
V(pq)   = x_pf(numel(pvpq)+1:end) .* exp(1j * angle(V(pq)));
I_pf = Ybus * V;
S_pf = V .* conj(I_pf);

%% 3. Kron reduction
% Dynamic ports: SM(1), GFM(2), GFL(3), and DC(8)
idx_p = [1, 2, 3, 8];
idx_l = [4, 5, 6, 7, 9];

Pd_net = Pd;
Qd_net = Qd;
Pd_net(8) = 0;
Qd_net(8) = 0; % Bus 8 load handled by the dynamic DC model

Ynet = Ybus + diag(conj(Pd_net + 1j * Qd_net) ./ (abs(V).^2));
Yeq = Ynet(idx_p, idx_p) - Ynet(idx_p, idx_l) * (Ynet(idx_l, idx_l) \ Ynet(idx_l, idx_p));
Vp0 = V(idx_p);

Ip0_gen = [conj((S_pf(1) + Pd(1) + 1j * Qd(1)) / Vp0(1));
           conj((S_pf(2) + Pd(2) + 1j * Qd(2)) / Vp0(2));
           conj(S3_pu / Vp0(3))];
Ip0_dc = conj((-Pd(8) - 1j * Qd(8)) / Vp0(4));
Ip0 = [Ip0_gen; Ip0_dc];

%% 4. Device initialization
sm = default_sm_params(baseMVA, fb);
gfm = default_gfm_params(fb, 1.0);
gfl = default_gfl_params(baseMVA, fb);
dc_param = default_datacenter_params(fb);

fprintf('\nInitializing Port Devices (Single DC at Bus 8)...\n');
[x0_sm, sm, refs_sm] = sm_init_from_pf_fsolve(Vp0(1), Ip0(1), sm);
[x0_gfm, refs_gfm] = gfm_init_from_pf_qss_newton(Vp0(2), Ip0(2), gfm);
[x0_gfl, gfl, refs_gfl] = gfl_init_from_pf_qss_newton(Vp0(3), Ip0(3), gfl);
[x0_dc1, refs_dc1] = datacenter_init_from_pf_qss(Vp0(4), Ip0(4), dc_param);

nx_sm  = numel(x0_sm);
nx_gfm = numel(x0_gfm);
nx_gfl = numel(x0_gfl);
nx_dc  = numel(x0_dc1);
nx = nx_sm + nx_gfm + nx_gfl + nx_dc;

v0_all = [real(Vp0(1)); imag(Vp0(1)); ...
          real(Vp0(2)); imag(Vp0(2)); ...
          real(Vp0(3)); imag(Vp0(3)); ...
          real(Vp0(4)); imag(Vp0(4))];
y0_final = [x0_sm(:); x0_gfm(:); x0_gfl(:); x0_dc1(:); v0_all(:)];

%% 5. Small-signal stability analysis
fprintf('\nSmall-signal stability analysis\n');
epsilon = 1e-7;
N_var = length(y0_final);
J = zeros(N_var, N_var);
net_base = struct('Yeq', Yeq, 'use_gpu_profile', false, 'exc_active', false);

f0 = dae_rhs(0, y0_final, nx_sm, nx_gfm, nx_gfl, nx_dc, ...
    sm, gfm, gfl, dc_param, refs_sm, refs_gfm, refs_gfl, refs_dc1, net_base);

for i = 1:N_var
    yp = y0_final;
    yp(i) = yp(i) + epsilon;
    f_p = dae_rhs(0, yp, nx_sm, nx_gfm, nx_gfl, nx_dc, ...
        sm, gfm, gfl, dc_param, refs_sm, refs_gfm, refs_gfl, refs_dc1, net_base);
    J(:, i) = (f_p - f0) / epsilon;
end

Am = J(1:nx, 1:nx);
Bm = J(1:nx, nx+1:end);
Cm = J(nx+1:end, 1:nx);
Dm = J(nx+1:end, nx+1:end);
Asys = Am - Bm * (Dm \ Cm);
Asys_red = Asys(2:end, 2:end); % remove SM delta state
[V_eig, D_eig] = eig(Asys_red);
lambda = diag(D_eig);
W_eig = inv(V_eig);

P_mat = abs(V_eig .* W_eig.');
for i = 1:size(P_mat,2)
    P_mat(:,i) = P_mat(:,i) / sum(P_mat(:,i));
end

% Reduced-model state indices
idx_sm_red  = 1 : (nx_sm - 1);
idx_gfm_red = nx_sm : (nx_sm + nx_gfm - 1);
idx_gfl_red = (nx_sm + nx_gfm) : (nx_sm + nx_gfm + nx_gfl - 1);
idx_dc_red  = (nx_sm + nx_gfm + nx_gfl) : size(P_mat, 1);

%% 6. POA gain analysis
fprintf('\nCalculating POA small-signal gain (0.1-1000 Hz)\n');
refs_pert = refs_dc1;
refs_pert.p_load = refs_dc1.p_load + epsilon;
f_u = dae_rhs(0, y0_final, nx_sm, nx_gfm, nx_gfl, nx_dc, ...
    sm, gfm, gfl, dc_param, refs_sm, refs_gfm, refs_gfl, refs_pert, net_base);
Em = (f_u(1:nx) - f0(1:nx)) / epsilon;
Fm = -(f_u(nx+1:end) - f0(nx+1:end)) / epsilon;
Hv = zeros(4, 8);
v_eq = y0_final(nx+1:end);
Vp_0 = v_eq(1:2:end) + 1j * v_eq(2:2:end);
P_base = zeros(4, 1);
P_base(1) = real(Vp_0(1) * conj(Yeq(1,:) * Vp_0));
P_base(2) = real(Vp_0(2) * conj(Yeq(2,:) * Vp_0));
P_base(3) = real(Vp_0(3) * conj(Yeq(3,:) * Vp_0));
P_base(4) = real(Vp_0(4) * conj(Yeq(4,:) * Vp_0));
for j = 1:8
    vp = v_eq;
    vp(j) = vp(j) + epsilon;
    Vp_p = vp(1:2:end) + 1j * vp(2:2:end);
    Hv(1, j) = (real(Vp_p(1) * conj(Yeq(1,:) * Vp_p)) - P_base(1)) / epsilon;
    Hv(2, j) = (real(Vp_p(2) * conj(Yeq(2,:) * Vp_p)) - P_base(2)) / epsilon;
    Hv(3, j) = (real(Vp_p(3) * conj(Yeq(3,:) * Vp_p)) - P_base(3)) / epsilon;
    Hv(4, j) = (real(Vp_p(4) * conj(Yeq(4,:) * Vp_p)) - P_base(4)) / epsilon;
end
Bsys = Em - Bm * (Dm \ Fm);
Csys = -Hv * (Dm \ Cm);
Dsys = -Hv * (Dm \ Fm);
f_scan = logspace(-1, log10(1000), 1000);
w_scan = 2 * pi * f_scan;
POA_mag = zeros(4, length(w_scan));
I_nx = eye(nx);
for k = 1:length(w_scan)
    resp = Csys * ((1j * w_scan(k) * I_nx - Asys) \ Bsys) + Dsys;
    POA_mag(:, k) = abs(resp);
end

fig_poa = figure('Name', 'POA Analysis - Multi-Port Gain');
semilogx(f_scan, POA_mag(1, :), 'LineWidth', 2); grid on; hold on;
semilogx(f_scan, POA_mag(2, :), 'LineWidth', 2);
semilogx(f_scan, POA_mag(3, :), 'LineWidth', 2);
semilogx(f_scan, POA_mag(4, :), 'k--', 'LineWidth', 2);
xlim([0.1, 1000]);
xlabel('Frequency (Hz)');
ylabel('POA');
legend('SM', 'GFM', 'GFL', 'data center', 'Location', 'best');
format_fig(fig_poa, 2, 4);

%% 7. Time-domain simulation with GPU profile
fprintf('\nLoading GPU data for time-domain simulation\n');
% Load and clean GPU trace
gpu_tbl = readtable('GPU_data.csv');
t_gpu_raw = gpu_tbl.t_seconds;
P_gpu_raw = gpu_tbl.P_gpu_W;
t_gpu_raw = t_gpu_raw(1:11301);
P_gpu_raw = P_gpu_raw(1:11301);

valid = isfinite(t_gpu_raw) & isfinite(P_gpu_raw);
t_gpu_raw = t_gpu_raw(valid);
P_gpu_raw = P_gpu_raw(valid);

t_gpu_raw = t_gpu_raw - t_gpu_raw(1);
[t_gpu_raw, uniq_idx] = unique(t_gpu_raw, 'stable');
P_gpu_raw = P_gpu_raw(uniq_idx);

% Define the simulation window
dt_load = median(diff(t_gpu_raw));
t_start = 1020;
t_end   = 1130.0;
t_load  = (t_start:dt_load:t_end).';

% 1. Extract the selected time window
idx_window = (t_gpu_raw >= 1020) & (t_gpu_raw <= 1130.0);
t_local = t_gpu_raw(idx_window);
P_local_raw = P_gpu_raw(idx_window);

% Shift time to start at 0 for the ODE solver
t_local = t_local - t_local(1);

% 2. Compute min/max within the selected window
P_local_min = min(P_local_raw);
P_local_max = max(P_local_raw);
fprintf('Local window limits: min = %.2f W, max = %.2f W\n', P_local_min, P_local_max);

% 3. Map the selected window to [0.3, 0.8] pu
P_load_samples = 0.3 + (0.8 - 0.3) * (P_local_raw - P_local_min) / (P_local_max - P_local_min);
t_load = t_local; % Use the time points from the selected window


fprintf('Mapped per-unit load on [%.1f, %.1f] s: min = %.3f, mean = %.3f, max = %.3f\n', ...
    t_start, t_end, min(P_load_samples), mean(P_load_samples), max(P_load_samples));

% Configure the ODE solver
M = blkdiag(speye(nx), sparse(8, 8));
opts_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'Mass', M);
net_sim = struct('Yeq', Yeq, 'use_gpu_profile', true, ...
                 't_load', t_load, 'P_load', P_load_samples);

fprintf('\nStarting Dynamic Simulation with real GPU profile (%.1f s to %.1f s)...\n', t_start, t_end);
[t, y] = ode15s(@(t,y) dae_rhs(t, y, nx_sm, nx_gfm, nx_gfl, nx_dc, ...
    sm, gfm, gfl, dc_param, refs_sm, refs_gfm, refs_gfl, refs_dc1, net_sim), ...
    t_load, y0_final, opts_ode);
%% 8. Port active-power plot
Nt = length(t);
t_stable = 90;
P_ports = zeros(Nt, 4);
P_plot  = zeros(Nt, 4);
for k = 1:Nt
    v_k = y(k, nx+1:end).';
    Vp_k = v_k(1:2:end) + 1j * v_k(2:2:end);
    I_net_k = Yeq * Vp_k;
    for m = 1:4
        P_ports(k, m) = real(Vp_k(m) * conj(I_net_k(m)));
    end
end
P_plot(:,1:3) = P_ports(:,1:3);
P_plot(:,4)   = -P_ports(:,4);

% Shift time so the plotted window starts at zero
t_plot = t - t_stable;
t_load_plot = t_load - t_stable;

fig_portP = figure('Name', 'Port Active Powers');
plot(t_plot, P_plot(:,1), 'LineWidth', 1.0); hold on; grid on;
plot(t_plot, P_plot(:,2), 'LineWidth', 1.0);
plot(t_plot, P_plot(:,3), 'LineWidth', 1.0);
plot(t_plot, P_plot(:,4), 'LineWidth', 1.0);
plot(t_load_plot, P_load_samples, 'k--', 'LineWidth', 0.8);

xlabel('Time (s)');
ylabel('PCC Active Power (p.u.)');
legend('SM', 'GFM', 'GFL', 'data center', 'Server Load', ...
       'NumColumns', 3, 'Location', 'north');

% Set x-axis limits for the shifted time window
xlim([0, (t_end - t_start) - t_stable]);
ylim([0, 1.5]);
format_fig(fig_portP, 2, 5);

%% 9. Eigenvalue map
fig4 = figure('Name', 'Eigenvalue Map');
N_modes = length(lambda);
LineWidths = zeros(N_modes, 1);
EdgeColors = cell(N_modes, 1);
P_DC_all  = zeros(N_modes, 1);
P_SM_all  = zeros(N_modes, 1);
P_GFM_all = zeros(N_modes, 1);
P_GFL_all = zeros(N_modes, 1);

lw_min = 0.5;
lw_max = 2.0;
dc_edge_th = 0.01;
mach_eps   = 1e-4;
for k = 1:N_modes
    p_sm_raw  = sum(P_mat(idx_sm_red,  k));
    p_gfm_raw = sum(P_mat(idx_gfm_red, k));
    p_gfl_raw = sum(P_mat(idx_gfl_red, k));
    p_dc      = sum(P_mat(idx_dc_red,  k));
    P_SM_all(k)  = p_sm_raw;
    P_GFM_all(k) = p_gfm_raw;
    P_GFL_all(k) = p_gfl_raw;
    P_DC_all(k)  = p_dc;
    if p_dc > dc_edge_th
        EdgeColors{k} = 'k';
        LineWidths(k) = lw_min + (lw_max - lw_min) * max(0, min(1, p_dc));
    else
        EdgeColors{k} = 'none';
        LineWidths(k) = 0.5;
    end
end

format_fig(fig4, 2, 10);
ax_main = axes('Position', [0.08 0.18 0.38 0.72]);
hold(ax_main, 'on'); grid(ax_main, 'on'); box(ax_main, 'on');
ax_in   = axes('Position', [0.56 0.18 0.38 0.72]);
hold(ax_in, 'on'); grid(ax_in, 'on'); box(ax_in, 'on');

xline(ax_main, 0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.5);
yline(ax_main, 0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.5);
xlabel(ax_main, '$\Re(\lambda)$', 'Interpreter', 'latex');
ylabel(ax_main, '$\Im(\lambda)$', 'Interpreter', 'latex');

xline(ax_in, 0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.5);
yline(ax_in, 0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.5);
xlabel(ax_in, '$\Re(\lambda)$', 'Interpreter', 'latex');
ylabel(ax_in, '$\Im(\lambda)$', 'Interpreter', 'latex');

xlim(ax_main, [min(real(lambda))-20, max(real(lambda))+20]);
y_margin = 0.08 * max(1, max(abs(imag(lambda))));
ylim(ax_main, [min(imag(lambda))-y_margin, max(imag(lambda))+y_margin]);

x_main = [-3000, 500];
x_zoom = [-500, 50];
xlim(ax_in, x_zoom);
xlim(ax_main, x_main);
roi = (real(lambda) >= x_zoom(1)) & (real(lambda) <= x_zoom(2));
if any(roi)
    yim_roi = imag(lambda(roi));
    yim_abs = sort(abs(yim_roi(:)));
    idx90 = max(1, ceil(0.90 * numel(yim_abs)));
    y90 = yim_abs(idx90);
    y_lim = max(20, 1.15 * y90);
    if numel(yim_abs) <= 4
        y_lim = max(20, 1.15 * max(yim_abs));
    end
    ylim(ax_in, [-y_lim, y_lim]);
else
    ylim(ax_in, [-100, 100]);
end

h_sm  = plot(ax_main, NaN, NaN, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'SM');
h_gfm = plot(ax_main, NaN, NaN, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'GFM');
h_gfl = plot(ax_main, NaN, NaN, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [0,1,0], 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'GFL');
h_dc  = plot(ax_main, NaN, NaN, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', 'k', ...
    'LineWidth', lw_max, 'DisplayName', 'data center');
legend(ax_main, [h_sm, h_gfm, h_gfl, h_dc], ...
    {'SM', 'GFM', 'GFL', 'data center'}, ...
    'Location', 'north', 'Interpreter', 'latex');

xin_lim = xlim(ax_in);
yin_lim = ylim(ax_in);
rectangle(ax_main, ...
    'Position', [xin_lim(1), yin_lim(1), xin_lim(2)-xin_lim(1), yin_lim(2)-yin_lim(1)], ...
    'EdgeColor', [0.4 0.4 0.4], ...
    'LineStyle', '-.', ...
    'LineWidth', 1.2);

drawnow;
pos_main  = get(ax_main, 'Position');
pos_in    = get(ax_in,   'Position');
xlim_main = xlim(ax_main);
ylim_main = ylim(ax_main);

norm_rect_x = pos_main(1) + ...
    (xin_lim(2) - xlim_main(1)) / (xlim_main(2) - xlim_main(1)) * pos_main(3);
norm_rect_y_top = pos_main(2) + ...
    (yin_lim(2) - ylim_main(1)) / (ylim_main(2) - ylim_main(1)) * pos_main(4);
norm_rect_y_bot = pos_main(2) + ...
    (yin_lim(1) - ylim_main(1)) / (ylim_main(2) - ylim_main(1)) * pos_main(4);

norm_ax2_x     = pos_in(1);
norm_ax2_y_top = pos_in(2) + pos_in(4);
norm_ax2_y_bot = pos_in(2);

annotation(fig4, 'line', [norm_rect_x, norm_ax2_x], [norm_rect_y_top, norm_ax2_y_top], ...
    'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.2);
annotation(fig4, 'line', [norm_rect_x, norm_ax2_x], [norm_rect_y_bot, norm_ax2_y_bot], ...
    'Color', [0.4 0.4 0.4], 'LineStyle', '-.', 'LineWidth', 1.2);

r_main = 30.0;
r_in   = 6;
render_modes_sector(ax_main, lambda, P_SM_all, P_GFM_all, P_GFL_all, ...
    EdgeColors, LineWidths, r_main, mach_eps);
render_modes_sector(ax_in, lambda, P_SM_all, P_GFM_all, P_GFL_all, ...
    EdgeColors, LineWidths, r_in, mach_eps);

uistack(findall(ax_main,'Type','Line'),'top');
uistack(findall(ax_in,'Type','Line'),'top');

%% Local functions
function F = pf_mismatch(x, Ybus, V, pvpq, pq, Psp, Qsp)
    Va = angle(V);
    Vm = abs(V);
    npvpq = numel(pvpq);
    Va(pvpq) = x(1:npvpq);
    Vm(pq) = x(npvpq+1:end);
    V_new = Vm .* exp(1j * Va);
    S = V_new .* conj(Ybus * V_new);
    F = [Psp(pvpq) - real(S(pvpq));
         Qsp(pq)   - imag(S(pq))];
end

function F = dae_rhs(t, y, nx_sm, nx_gfm, nx_gfl, nx_dc, ...
    sm, gfm, gfl, dc, refs_sm, refs_gfm, refs_gfl, r1, net)

    nx_total = nx_sm + nx_gfm + nx_gfl + nx_dc;
    x = y(1:nx_total);
    v = y(nx_total+1:end);
    Vp = v(1:2:end) + 1j * v(2:2:end);
    omega_sm = x(2);

    rt = r1;

    if isfield(net, 'use_gpu_profile') && net.use_gpu_profile
        rt.p_load = interp1(net.t_load, net.P_load, t, 'previous', 'extrap');
    elseif isfield(net, 'exc_active') && net.exc_active
        rt.p_load = r1.p_load * (1 + net.A_exc * sin(2*pi*net.f_exc * (t - net.t_start)));
    end

    [dx_sm,  I1] = sm_port_model(t, x(1:nx_sm), Vp(1), sm, refs_sm.pref);
    dx_sm(1) = 0;

    [dx_gfm, I2] = gfm_port_model(t, x(nx_sm+1:nx_sm+nx_gfm), ...
        Vp(2), gfm, refs_gfm.pref0, refs_gfm.qref0, refs_gfm.ome_ref0, refs_gfm.v_ref0, omega_sm);

    [dx_gfl, I3] = gfl_port_model(t, x(nx_sm+nx_gfm+1:nx_sm+nx_gfm+nx_gfl), ...
        Vp(3), gfl, refs_gfl.pref0, refs_gfl.qref0, omega_sm);

    [dx_dc,  I4] = datacenter_port_model(t, x(nx_sm+nx_gfm+nx_gfl+1:end), ...
        Vp(4), dc, rt, omega_sm);

    g = [I1; I2; I3; I4] - net.Yeq * Vp;
    F = [dx_sm; dx_gfm; dx_gfl; dx_dc; -real(g); -imag(g)];
end

function P_out = local_tile_signal(t_in, P_in, t_out)
    t_mod = mod(t_out - t_out(1), max(t_in));
    P_out = interp1(t_in, P_in, t_mod, 'previous', 'extrap');
end

function sm = default_sm_params(baseMVA, fb)
    sm = struct( ...
        'omega_b', 2*pi*fb, 'ws', 1.0, ...
        'xd', 0.1460, 'xq', 0.0969, 'xdp', 0.0608, 'xqp', 0.0969, ...
        'Td0p', 8.96, 'Tq0p', 0.31, 'H', 3.0, 'D', 0.0, ...
        'ka', 5.0, 'Ta', 0.2, 'Te', 0.314, 'kf', 0.063, 'Tf', 0.35, ...
        'ke', 1.0, 'ae', 0.0039, 'be', 1.555, ...
        'r', 0.15, 'Tsv', 0.1, 'Tch', 0.5);
end

function gfm = default_gfm_params(fb, ome_s)
    ome_b = 2*pi*fb;
    wc_v = 2*pi*80;
    wc_i = 2*pi*200;
    cf = 0.074;
    lf = 0.08;
    gfm = struct( ...
        'base', struct('ome_sys', ome_s, 'ome_b', ome_b), ...
        'filt', struct('lf', lf, 'rf', 0.003, 'cf', cf, 'lg', 0.2, 'rg', 0.01), ...
        'droop', struct('Kp', 0.02, 'ome_z', 20.0, 'Kq', 0.05, 'ome_f', 50.0), ...
        'virt', struct('rv', 0.0, 'lv', 0.2), ...
        'inner', struct('kpv', 4.0*wc_v*cf/ome_b, 'kiv', (wc_v^2)*cf/ome_b, ...
                        'kpc', 1.414*wc_i*lf/ome_b, 'kic', (wc_i^2)*lf/ome_b), ...
        'dc', struct('Vdc', 1.0));
end

function gfl = default_gfl_params(baseMVA, fb)
    base = struct('baseMVA', baseMVA, 'fb_hz', fb, 'ome_b', 2*pi*fb, 'ome_sys', 1.0);
    gfl = struct( ...
        'base', base, ...
        'filt', struct('lf', 0.08, 'rf', 0.003, 'cf', 0.074, 'lg', 0.1, 'rg', 0.01), ...
        'pll', struct('kp', 0.05, 'ki', 1.42, 'ome_lp', 376.99), ...
        'outer', struct('Kp_p', 0.01, 'Ki_p', 0.12, 'Kp_q', 0.01, 'Ki_q', 0.12, ...
                        'ome_z', 41.47, 'ome_f', 41.47), ...
        'inner', struct('kpc', 0.15, 'kic', 0.267, 'kffv', 0, 'lf', 0.08), ...
        'dc', struct('Vdc', 1.0));
end

function dc = default_datacenter_params(fb)
    dc = struct( ...
        'fb', fb, 'omega_b', 2*pi*fb, 'omega_s', 1.0, 'omega_lp', 2*pi*100, ...
        'l_afe', 0.05, 'r_afe', 0.003, 'c_dc', 2.0, 'vdc_ups_ref', 1.0, ...
        'l_vsi', 0.05, 'r_vsi', 0.003, 'c_vsi', 0.2, 'vu_vsi_ref', 1.0, 'omega_vsi', 1.0, ...
        'c_psu', 2, 'r_psu', 0.005, 'vpsu_ref', 1.0, ...
        'c_eq', 0.2, 'veq_ref', 0.5, ...
        'kp_pll', 0.471, 'ki_pll', 41.89);
    wn = 2*pi*5;
    dc.kp_dc_afe = 2*1.0*wn*(dc.c_dc/dc.omega_b);
    dc.ki_dc_afe = (wn^2)*(dc.c_dc/dc.omega_b);
    wn = 2*pi*200;
    dc.kp_c_afe = 2*0.707*wn*(dc.l_afe/dc.omega_b) - dc.r_afe;
    dc.ki_c_afe = (wn^2)*(dc.l_afe/dc.omega_b);
    wn = 2*pi*100;
    dc.kp_v_vsi = 2*1.0*wn*(dc.c_vsi/dc.omega_b);
    dc.ki_v_vsi = (wn^2)*(dc.c_vsi/dc.omega_b);
    wn = 2*pi*400;
    dc.kp_c_vsi = 2*1.0*wn*(dc.l_vsi/dc.omega_b) - dc.r_vsi;
    dc.ki_c_vsi = (wn^2)*(dc.l_vsi/dc.omega_b);
    wn = 2*pi*10;
    dc.kp_v_psu = 2*1.0*wn*(dc.c_psu/dc.omega_b);
    dc.ki_v_psu = (wn^2)*(dc.c_psu/dc.omega_b);
    wn = 2*pi*100;
    dc.kp_v_eq = 2*1.0*wn*(dc.c_eq/dc.omega_b);
    dc.ki_v_eq = (wn^2)*(dc.c_eq/dc.omega_b);
end

function render_modes_sector(ax, lambda, P_SM_all, P_GFM_all, P_GFL_all, ...
    EdgeColors, LineWidths, r_x, mach_eps)
    col_sm  = [0, 0, 1];
    col_gfm = [1, 0, 0];
    col_gfl = [0, 1, 0];
    r_y = get_visual_circle_yradius(ax, r_x);
    for k = 1:length(lambda)
        x0 = real(lambda(k));
        y0 = imag(lambda(k));
        parts = [P_SM_all(k), P_GFM_all(k), P_GFL_all(k)];
        sum_dyn = sum(parts);
        if sum_dyn <= mach_eps
            draw_white_marker(ax, x0, y0, r_x, r_y, EdgeColors{k}, LineWidths(k));
        else
            weights = parts / sum_dyn;
            draw_sector_marker(ax, x0, y0, r_x, r_y, weights, ...
                {col_sm, col_gfm, col_gfl}, EdgeColors{k}, LineWidths(k));
        end
    end
end

function r_y = get_visual_circle_yradius(ax, r_x)
    drawnow;
    old_units = get(ax, 'Units');
    set(ax, 'Units', 'pixels');
    ax_pos = get(ax, 'Position');
    set(ax, 'Units', old_units);
    xl = xlim(ax);
    yl = ylim(ax);
    x_per_pix = diff(xl) / ax_pos(3);
    y_per_pix = diff(yl) / ax_pos(4);
    r_pix = r_x / x_per_pix;
    r_y = r_pix * y_per_pix;
end

function draw_white_marker(ax, x0, y0, r_x, r_y, edgeColor, lineWidth)
    th = linspace(0, 2*pi, 100);
    xp = x0 + r_x*cos(th);
    yp = y0 + r_y*sin(th);
    patch(ax, xp, yp, [1 1 1], ...
        'EdgeColor', edgeColor, ...
        'LineWidth', lineWidth, ...
        'Clipping', 'on', ...
        'HandleVisibility', 'off');
end

function draw_sector_marker(ax, x0, y0, r_x, r_y, weights, colors, edgeColor, lineWidth)
    weights = max(weights(:).', 0);
    s = sum(weights);
    if s <= 1e-12
        draw_white_marker(ax, x0, y0, r_x, r_y, edgeColor, lineWidth);
        return;
    end
    weights = weights / s;
    theta_start = pi/2;
    for i = 1:numel(weights)
        wi = weights(i);
        if wi <= 1e-8
            continue;
        end
        theta_end = theta_start + 2*pi*wi;
        th = linspace(theta_start, theta_end, 40);
        xp = [x0, x0 + r_x*cos(th), x0];
        yp = [y0, y0 + r_y*sin(th), y0];
        patch(ax, xp, yp, colors{i}, ...
            'EdgeColor', 'none', ...
            'LineWidth', 0.5, ...
            'Clipping', 'on', ...
            'HandleVisibility', 'off');
        theta_start = theta_end;
    end
    if ~isequal(edgeColor, 'none')
        th = linspace(0, 2*pi, 120);
        xp = x0 + r_x*cos(th);
        yp = y0 + r_y*sin(th);
        patch(ax, xp, yp, 'w', ...
            'FaceColor', 'none', ...
            'EdgeColor', edgeColor, ...
            'LineWidth', lineWidth, ...
            'Clipping', 'on', ...
            'HandleVisibility', 'off');
    end
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
