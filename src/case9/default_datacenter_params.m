function dc = default_datacenter_params(fb)
%DEFAULT_DATACENTER_PARAMS Return parameters and controller gains for the data-center model.

    dc = struct();
    dc.fb = fb;
    dc.omega_b = 2 * pi * fb;
    dc.omega_s = 1.0;

    % AFE
    dc.omega_lp    = 2 * pi * 100;
    dc.l_afe       = 0.05;
    dc.r_afe       = 0.003;
    dc.c_dc        = 1.0;
    dc.vdc_ups_ref = 1.0;

    % VSI
    dc.l_vsi       = 0.05;
    dc.r_vsi       = 0.003;
    dc.c_vsi       = 0.2;
    dc.vu_vsi_ref  = 1.0;
    dc.omega_vsi   = 1.0;

    % PSU and EQ
    dc.c_psu    = 0.3;
    dc.r_psu    = 0.005;
    dc.vpsu_ref = 1.0;
    dc.c_eq     = 0.2;
    dc.veq_ref  = 0.5;

    % PI Tuning via Bandwidth
    % PLL: fbw=20, zeta=0.707
    wn_pll = 2 * pi * 20.0;
    dc.kp_pll = (2 * 0.707 * wn_pll) / dc.omega_b;
    dc.ki_pll = (wn_pll^2) / dc.omega_b;

    % Voltage loops helper
    tune_v = @(fbw, zeta, C) struct('kp', 2*zeta*(2*pi*fbw)*(C/dc.omega_b), 'ki', ((2*pi*fbw)^2)*(C/dc.omega_b));
    % Current loops helper
    tune_c = @(fbw, zeta, L, R) struct('kp', 2*zeta*(2*pi*fbw)*(L/dc.omega_b) - R, 'ki', ((2*pi*fbw)^2)*(L/dc.omega_b));

    % Apply tuning
    dc_afe_pi = tune_v(10.0, 1.0, dc.c_dc);
    dc.kp_dc_afe = dc_afe_pi.kp; dc.ki_dc_afe = dc_afe_pi.ki;

    c_afe_pi = tune_c(200.0, 0.707, dc.l_afe, dc.r_afe);
    dc.kp_c_afe = c_afe_pi.kp; dc.ki_c_afe = c_afe_pi.ki;

    v_vsi_pi = tune_v(100.0, 1.0, dc.c_vsi);
    dc.kp_v_vsi = v_vsi_pi.kp; dc.ki_v_vsi = v_vsi_pi.ki;

    c_vsi_pi = tune_c(400.0, 1.0, dc.l_vsi, dc.r_vsi);
    dc.kp_c_vsi = c_vsi_pi.kp; dc.ki_c_vsi = c_vsi_pi.ki;

    v_psu_pi = tune_v(20.0, 1.0, dc.c_psu);
    dc.kp_v_psu = v_psu_pi.kp; dc.ki_v_psu = v_psu_pi.ki;

    v_eq_pi = tune_v(100.0, 1.0, dc.c_eq);
    dc.kp_v_eq = v_eq_pi.kp; dc.ki_v_eq = v_eq_pi.ki;
end
