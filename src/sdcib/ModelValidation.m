%% Model Validation: Full abc Model versus QSS uv Model
% Compares the full cascaded system and reduced quasi-steady-state
% representation using a fixed-step trapezoidal simulation.
clear; clc; close all;

% 1. System Parameters
wb = 2 * pi * 60;
omega_s = 1.0;
omega_lp = 2 * pi * 100;
Vinf = 1.0; Rinf = 0.02; Xinf = 0.19;

l_afe = 0.05; r_afe = 0.003; c_dc = 2.0; vdc_ups_ref = 1.0;
l_vsi = 0.05; r_vsi = 0.003; c_vsi = 0.2; vu_vsi_ref = 1.0;

l_psu = 0.05; c_psu = 2.0; r_psu = 0.005; v_psu_ref = 1.0;
l_eq  = 0.05; c_eq  = 0.2; v_eq_ref  = 0.5;

% 2. Controller PI Tuning
[kp_pll, ki_pll]       = pll_bw_tuning(20.0, 0.707);
[kp_dc_afe, ki_dc_afe] = pi_bw_tuning('v', 5.0, 1.0, c_dc, wb, 0);
[kp_c_afe, ki_c_afe]   = pi_bw_tuning('c', 200.0, 0.707, l_afe, wb, r_afe);
[kp_v_vsi, ki_v_vsi]   = pi_bw_tuning('v', 100.0, 1.0, c_vsi, wb, 0);
[kp_c_vsi, ki_c_vsi]   = pi_bw_tuning('c', 400.0, 1.0, l_vsi, wb, r_vsi);
[kp_v_psu, ki_v_psu]   = pi_bw_tuning('v', 10.0, 1.0, c_psu, wb, 0);
[kp_c_psu, ki_c_psu]   = pi_bw_tuning('c', 1000.0, 1.0, l_psu, wb, r_psu)
[kp_v_eq, ki_v_eq]     = pi_bw_tuning('v', 100.0, 1.0, c_eq, wb, 0);
[kp_c_eq, ki_c_eq]     = pi_bw_tuning('c', 1000.0, 1.0, l_eq, wb, 0)

% Pack parameters
params = {wb, omega_s, omega_lp, Vinf, Rinf, Xinf, l_afe, r_afe, c_dc, vdc_ups_ref, ...
          l_vsi, r_vsi, c_vsi, vu_vsi_ref, l_psu, c_psu, r_psu, v_psu_ref, l_eq, c_eq, v_eq_ref, ...
          kp_pll, ki_pll, kp_dc_afe, ki_dc_afe, kp_c_afe, ki_c_afe, kp_v_vsi, ki_v_vsi, ...
          kp_c_vsi, ki_c_vsi, kp_v_psu, ki_v_psu, kp_c_psu, ki_c_psu, kp_v_eq, ki_v_eq, kp_c_eq, ki_c_eq};

% 3. Initial Conditions
p_load_init = 0.5;
g_load0 = p_load_init / (3 * v_eq_ref^2);

X = zeros(62, 1);
for offset = [0, 41]
    X(offset+4)  = p_load_init;
    X(offset+6)  = p_load_init / ki_dc_afe;
    X(offset+9)  = vdc_ups_ref;
    X(offset+10) = p_load_init / vu_vsi_ref;
    X(offset+12) = vu_vsi_ref;
    X(offset+14) = p_load_init / vu_vsi_ref / ki_v_vsi;
end

i_eq0 = g_load0 * v_eq_ref;
i_psu0 = (v_eq_ref / v_psu_ref) * i_eq0;
g_eq0 = (1 - sqrt(1 - 4 * r_psu * (i_psu0 * 3))) / (2 * r_psu);

% System A (Full ABC Load)
X(18:20) = g_eq0 * sqrt((sqrt(2/3)*vu_vsi_ref*[1; -0.5; -0.5]).^2 + 1e-6);
X(21:23) = v_psu_ref;
X(24:26) = g_eq0 / ki_v_psu;
X(30:32) = i_eq0;
X(33:35) = v_eq_ref;
X(36:38) = i_eq0 / ki_v_eq;

% System B (QSS UV Load)
X(59) = v_psu_ref; X(60) = g_eq0 / ki_v_psu; X(61) = v_eq_ref; X(62) = i_eq0 / ki_v_eq;

% 4. Fixed-Step Solver Setup
dt = 5e-6;
t_end = 1.55;
t_arr = 0:dt:t_end;
N = length(t_arr);

V_psu_abc_res = zeros(N, 3); V_psu_uv_res  = zeros(N, 1);
V_eq_abc_res  = zeros(N, 3); V_eq_uv_res   = zeros(N, 1);
P_vsi_A_res   = zeros(N, 1); P_vsi_B_res   = zeros(N, 1);
P_pcc_A_res   = zeros(N, 1); P_pcc_B_res   = zeros(N, 1);

fprintf('Running Trapezoidal Simulation (with VSI Power Logging)...\n');

for k = 1:N-1
    t = t_arr(k);

    % Log state variables
    V_psu_abc_res(k, :) = X(21:23)';
    V_psu_uv_res(k)     = X(59);
    V_eq_abc_res(k, :)  = X(33:35)';
    V_eq_uv_res(k)      = X(61);

    if t < 1.025, p_load = 0.50; else, p_load = 0.60; end

    % Step 1: Predictor
    [dX1, PA1, PB1, PPCCA1, PPCCB1] = get_sys_derivatives(t, X, p_load, params);

    % Log Power (At the start of the step)
    P_vsi_A_res(k)  = PA1;
    P_vsi_B_res(k)  = PB1;
    P_pcc_A_res(k)  = PPCCA1;
    P_pcc_B_res(k)  = PPCCB1;

    X_pred = X + dt * dX1;
    X_pred(18:20) = max(0, X_pred(18:20));
    X_pred(30:32) = max(0, X_pred(30:32));

    % Step 2: Corrector
    t_next = t + dt;
    if t_next < 1.025, p_load_next = 0.50; else, p_load_next = 0.60; end
    [dX2, ~, ~, ~, ~] = get_sys_derivatives(t_next, X_pred, p_load_next, params);

    % Step 3: Trapezoidal Update
    X = X + (dt / 2) * (dX1 + dX2);
    X(18:20) = max(0, X(18:20));
    X(30:32) = max(0, X(30:32));
end

% Log final point
V_psu_abc_res(N, :) = X(21:23)'; V_psu_uv_res(N) = X(59);
V_eq_abc_res(N, :)  = X(33:35)'; V_eq_uv_res(N)  = X(61);
[~, PA_end, PB_end, PPCCA_end, PPCCB_end] = get_sys_derivatives(t_end, X, p_load, params);
P_vsi_A_res(N)  = PA_end;
P_vsi_B_res(N)  = PB_end;
P_pcc_A_res(N)  = PPCCA_end;
P_pcc_B_res(N)  = PPCCB_end;

fprintf('Simulation Complete! Plotting with IEEE Format...\n');

% 5. Plotting (IEEE Single Column Format)
plot_idx = t_arr >= 1.0;
t_plot = t_arr(plot_idx);

% Plot 1: Upstream PSU DC Voltage
fig1 = figure('Name', 'PSU DC Voltage');
plot(t_plot, V_psu_abc_res(plot_idx, 1), 'r', 'DisplayName', 'Phase A (abc)'); hold on;
plot(t_plot, V_psu_abc_res(plot_idx, 2), 'g', 'DisplayName', 'Phase B (abc)');
plot(t_plot, V_psu_abc_res(plot_idx, 3), 'b', 'DisplayName', 'Phase C (abc)');
plot(t_plot, V_psu_uv_res(plot_idx), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Aggregated (uv)');
xlabel('Time (s)', 'Interpreter', 'latex'); ylabel('$v^{\mathrm{psu}}$ (p.u.)', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
xlim([1.0, 1.2]); ylim([0, 1.5]);
format_ieee_singlecol(fig1, 2, 4);

% Plot 2: Downstream Load DC Voltage
fig2 = figure('Name', 'Load DC Voltage');
plot(t_plot, V_eq_abc_res(plot_idx, 1), 'r', 'DisplayName', 'Phase A (abc)'); hold on;
plot(t_plot, V_eq_abc_res(plot_idx, 2), 'g', 'DisplayName', 'Phase B (abc)');
plot(t_plot, V_eq_abc_res(plot_idx, 3), 'b', 'DisplayName', 'Phase C (abc)');
plot(t_plot, V_eq_uv_res(plot_idx), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Aggregated (uv)');
xlabel('Time (s)', 'Interpreter', 'latex'); ylabel('$v^{\mathrm{eq}}$ (p.u.)', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
xlim([1.0, 1.5]); ylim([0.3, 0.6]);
format_ieee_singlecol(fig2, 2, 4);

% Plot 3: VSI Total Extracted Power
fig3 = figure('Name', 'VSI Extracted Power');
plot(t_plot, P_vsi_A_res(plot_idx), 'r', 'LineWidth', 1.0, 'DisplayName', 'Full Model ($abc$)'); hold on;
plot(t_plot, P_vsi_B_res(plot_idx), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reduced Model ($uv$)');
xlabel('Time (s)', 'Interpreter', 'latex'); ylabel('$p_{\mathrm{vsi}}$ (p.u.)', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
xlim([1.0, 1.5]); ylim([0.45, 0.7]);
format_ieee_singlecol(fig3, 2, 4);

% Plot 4: PCC Total Extracted Power
fig4 = figure('Name', 'PCC Extracted Power');
plot(t_plot, P_pcc_A_res(plot_idx), 'r', 'LineWidth', 1.0, 'DisplayName', 'Full Model ($abc$)'); hold on;
plot(t_plot, P_pcc_B_res(plot_idx), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reduced Model ($uv$)');
xlabel('Time (s)', 'Interpreter', 'latex'); ylabel('$p_{\mathrm{pcc}}$ (p.u.)', 'Interpreter', 'latex');
legend('Location', 'southeast', 'Interpreter', 'latex');
xlim([1.0, 1.5]); ylim([0.45, 0.7]);
format_ieee_singlecol(fig4, 2, 4);


% LOCAL FUNCTIONS

function [kp, ki] = pi_bw_tuning(type, fbw, zeta, p1, wb, p2)
    wn = 2*pi*fbw;
    if type == 'v', kp = 2*zeta*wn*(p1/wb); ki = (wn^2)*(p1/wb);
    else, kp = 2*zeta*wn*(p1/wb) - p2; ki = (wn^2)*(p1/wb); end
end

function [kp, ki] = pll_bw_tuning(fbw, zeta)
    wn = 2*pi*fbw; kp = (2*zeta*wn)/(2*pi*60); ki = (wn^2)/(2*pi*60);
end

function y = safe_div(a,b)
    if abs(b) < 1e-9, y = a/(sign(b+1e-16)*1e-9); else, y = a/b; end
end

% IEEE figure formatting
function format_ieee_singlecol(fig, h_in, w_in)
    if nargin < 2 || isempty(h_in), h_in = 2.0; end
    if nargin < 3 || isempty(w_in), w_in = 3.5; end % Standard IEEE single column width
    set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 w_in h_in], 'PaperUnits', 'inches', 'PaperPosition', [0 0 w_in h_in], 'PaperSize', [w_in h_in]);
    ax = findall(fig, 'Type', 'axes');
    for k = 1:numel(ax)
        set(ax(k), 'FontName','Times New Roman', 'FontSize',8, 'LineWidth',0.8, 'Box','on', 'TickDir','in');
        grid(ax(k),'on'); ax(k).GridAlpha = 0.25; ax(k).MinorGridAlpha = 0.15;
    end
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

% Unified 62-state DAE solver with instantaneous power outputs
function [dX, P_A, P_B, P_pcc_A, P_pcc_B] = get_sys_derivatives(t, X, p_load, params)
    [wb, omega_s, omega_lp, Vinf, Rinf, Xinf, l_afe, r_afe, c_dc, vdc_ups_ref, ...
     l_vsi, r_vsi, c_vsi, vu_vsi_ref, l_psu, c_psu, r_psu, v_psu_ref, l_eq, c_eq, v_eq_ref, ...
     kp_pll, ki_pll, kp_dc_afe, ki_dc_afe, kp_c_afe, ki_c_afe, kp_v_vsi, ki_v_vsi, ...
     kp_c_vsi, ki_c_vsi, kp_v_psu, ki_v_psu, kp_c_psu, ki_c_psu, kp_v_eq, ki_v_eq, kp_c_eq, ki_c_eq] = params{:};

    dX = zeros(62, 1);
    theta_vsi = wb * omega_s * t;
    g_load = p_load / (3 * v_eq_ref^2);

    sqrt23 = sqrt(2/3);
    c0 = cos(theta_vsi); s0 = sin(theta_vsi);
    c1 = cos(theta_vsi - 2*pi/3); s1 = sin(theta_vsi - 2*pi/3);
    c2 = cos(theta_vsi + 2*pi/3); s2 = sin(theta_vsi + 2*pi/3);
    T_inv = sqrt23 * [c0, -s0; c1, -s1; c2, -s2];
    T_fwd = sqrt23 * [c0, c1, c2; -s0, -s1, -s2];

    % Sub-system A: Full ABC Model

    xA = X(1:41); dXA = zeros(41,1);
    i_rec_abc=xA(18:20); v_psu_abc=xA(21:23); xi_psu_abc=xA(24:26); gam_psu=xA(27:29);
    i_eq_abc=xA(30:32); v_eq_abc=xA(33:35); xi_eq_abc=xA(36:38); gam_eq=xA(39:41);

    v_abc = T_inv * [xA(12); xA(13)];
    v_rec = sqrt(v_abc.^2 + 1e-6);

    i_eq_ref = kp_v_eq*(v_eq_ref - v_eq_abc) + ki_v_eq*xi_eq_abc;
    d_eq_unl = kp_c_eq*(i_eq_ref - i_eq_abc) + ki_c_eq*gam_eq;
    d_eq     = max(0, min(1.0, d_eq_unl));
    dgam_eq  = i_eq_ref - i_eq_abc;
    dgam_eq((d_eq_unl >= 1.0 & dgam_eq > 0) | (d_eq_unl <= 0 & dgam_eq < 0)) = 0;

    dXA(30:32) = (wb/l_eq)*(d_eq.*v_psu_abc - v_eq_abc);
    dXA(33:35) = (wb/c_eq)*(i_eq_abc - g_load*v_eq_abc);
    dXA(36:38) = v_eq_ref - v_eq_abc;
    dXA(39:41) = dgam_eq;

    i_psu_ld = d_eq .* i_eq_abc;
    g_eq_abc = kp_v_psu*(v_psu_ref - v_psu_abc) + ki_v_psu*xi_psu_abc;
    i_ref_psu= g_eq_abc .* v_rec;
    d_psu_unl= kp_c_psu*(i_ref_psu - i_rec_abc) + ki_c_psu*gam_psu;
    d_psu    = max(0, min(0.99, d_psu_unl));
    dgam_psu = i_ref_psu - i_rec_abc;
    dgam_psu((d_psu_unl >= 0.99 & dgam_psu > 0) | (d_psu_unl <= 0 & dgam_psu < 0)) = 0;

    di_rec = (wb/l_psu)*(v_rec - (1-d_psu).*v_psu_abc - r_psu*i_rec_abc);
    di_rec((i_rec_abc <= 0) & (di_rec < 0)) = 0;
    dXA(18:20) = di_rec;
    dXA(21:23) = (wb/c_psu)*((1-d_psu).*i_rec_abc - i_psu_ld);
    dXA(24:26) = v_psu_ref - v_psu_abc;
    dXA(27:29) = dgam_psu;

    iU_vsi_ld_A = T_fwd(1,:) * ((v_abc ./ v_rec) .* i_rec_abc);
    iV_vsi_ld_A = T_fwd(2,:) * ((v_abc ./ v_rec) .* i_rec_abc);

    % Total instantaneous power extracted from VSI AC bus (System A)
    P_A = xA(12) * iU_vsi_ld_A + xA(13) * iV_vsi_ld_A;

    % Sub-system B: QSS UV Model

    xB = X(42:62); dXB = zeros(21,1);
    v_psu_uv=xB(18); xi_psu_uv=xB(19); v_eq_uv=xB(20); xi_eq_uv=xB(21);

    i_eq_qss = kp_v_eq*(v_eq_ref - v_eq_uv) + ki_v_eq*xi_eq_uv;
    dXB(20) = (wb/c_eq)*(i_eq_qss - g_load*v_eq_uv);
    dXB(21) = v_eq_ref - v_eq_uv;

    i_psu_ld_uv = (v_eq_uv / max(0.1, v_psu_uv)) * i_eq_qss;
    g_eq_uv = kp_v_psu*(v_psu_ref - v_psu_uv) + ki_v_psu*xi_psu_uv;
    norm_v_sq = xB(12)^2 + xB(13)^2;
    dXB(18) = (wb/c_psu)*( (g_eq_uv - r_psu*g_eq_uv^2)*norm_v_sq/(3*max(v_psu_uv,1e-9)) - i_psu_ld_uv );
    dXB(19) = v_psu_ref - v_psu_uv;

    iU_vsi_ld_B = g_eq_uv * xB(12);
    iV_vsi_ld_B = g_eq_uv * xB(13);

    % Total instantaneous power extracted from VSI AC bus (System B)
    P_B = xB(12) * iU_vsi_ld_B + xB(13) * iV_vsi_ld_B;

    % Shared Grid/VSI Dynamics

    iU_ld_arr = [iU_vsi_ld_A, iU_vsi_ld_B];
    iV_ld_arr = [iV_vsi_ld_A, iV_vsi_ld_B];

    P_pcc_A = 0;
    P_pcc_B = 0;

    for sys = 1:2
        if sys == 1, xs = xA; iU_ld = iU_ld_arr(1); iV_ld = iV_ld_arr(1);
        else,        xs = xB; iU_ld = iU_ld_arr(2); iV_ld = iV_ld_arr(2); end

        iU_ref_cv = kp_v_vsi*(vu_vsi_ref - xs(12)) + ki_v_vsi*xs(14) - omega_s*c_vsi*xs(13);
        iV_ref_cv = kp_v_vsi*(0 - xs(13)) + ki_v_vsi*xs(15) + omega_s*c_vsi*xs(12);
        vU_ref_cv = kp_c_vsi*(iU_ref_cv - xs(10)) + ki_c_vsi*xs(16) - omega_s*l_vsi*xs(11);
        vV_ref_cv = kp_c_vsi*(iV_ref_cv - xs(11)) + ki_c_vsi*xs(17) + omega_s*l_vsi*xs(10);
        mU = safe_div(vU_ref_cv, xs(9)); mV = safe_div(vV_ref_cv, xs(9));

        omega_pll = omega_s + kp_pll*xs(3) + ki_pll*xs(2);
        s_ang = sin(xs(1) + pi/2); c_ang = cos(xs(1) + pi/2);
        ir_pcc = s_ang*xs(4) + c_ang*xs(5); ii_pcc = -c_ang*xs(4) + s_ang*xs(5);
        vr_pcc = Vinf - Rinf*ir_pcc + Xinf*ii_pcc; vi_pcc = 0 - Rinf*ii_pcc - Xinf*ir_pcc;
        vd_pcc = s_ang*vr_pcc - c_ang*vi_pcc; vq_pcc = c_ang*vr_pcc + s_ang*vi_pcc;

        P_pcc = vd_pcc * xs(4) + vq_pcc * xs(5);

        id_ref_afe = kp_dc_afe*(vdc_ups_ref - xs(9)) + ki_dc_afe*xs(6); iq_ref_afe = 0;
        vd_ref_afe = kp_c_afe*(xs(4) - id_ref_afe) + ki_c_afe*xs(7) + omega_pll*l_afe*xs(5);
        vq_ref_afe = kp_c_afe*(xs(5) - iq_ref_afe) + ki_c_afe*xs(8) - omega_pll*l_afe*xs(4);
        md = safe_div(vd_ref_afe, xs(9)); mq = safe_div(vq_ref_afe, xs(9));

        dx_grid = zeros(17,1);
        dx_grid(1) = wb*(omega_pll - omega_s); dx_grid(2) = xs(3); dx_grid(3) = omega_lp*(vq_pcc - xs(3));
        dx_grid(4) = (wb/l_afe)*(vd_pcc - md*xs(9) - r_afe*xs(4) + omega_pll*l_afe*xs(5));
        dx_grid(5) = (wb/l_afe)*(vq_pcc - mq*xs(9) - r_afe*xs(5) - omega_pll*l_afe*xs(4));
        dx_grid(6) = vdc_ups_ref - xs(9); dx_grid(7) = xs(4) - id_ref_afe; dx_grid(8) = xs(5) - iq_ref_afe;
        dx_grid(9) = (wb/c_dc)*((md*xs(4) + mq*xs(5)) - (mU*xs(10) + mV*xs(11)));
        dx_grid(10) = (wb/l_vsi)*(mU*xs(9) - xs(12) - r_vsi*xs(10) + omega_s*l_vsi*xs(11));
        dx_grid(11) = (wb/l_vsi)*(mV*xs(9) - xs(13) - r_vsi*xs(11) - omega_s*l_vsi*xs(10));
        dx_grid(12) = (wb/c_vsi)*(xs(10) - iU_ld + omega_s*c_vsi*xs(13));
        dx_grid(13) = (wb/c_vsi)*(xs(11) - iV_ld - omega_s*c_vsi*xs(12));
        dx_grid(14) = vu_vsi_ref - xs(12); dx_grid(15) = 0 - xs(13);
        dx_grid(16) = iU_ref_cv - xs(10); dx_grid(17) = iV_ref_cv - xs(11);

        if sys == 1
            dXA(1:17) = dx_grid;
            P_pcc_A = P_pcc;
        else
            dXB(1:17) = dx_grid;
            P_pcc_B = P_pcc;
        end
    end

    dX(1:41)  = dXA;
    dX(42:62) = dXB;
end
