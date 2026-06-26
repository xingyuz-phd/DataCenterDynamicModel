function [xdot, Iport, aux] = datacenter_port_model(t, x, Vport, dc, refs, omega_sys)
%DATACENTER_PORT_MODEL Dynamic port model for the data-center load.
% Current injection follows the network convention used by the reduced
% admittance matrix: positive current is injected into the grid.

    % Handle reference frequency fallback for initialization
    if nargin < 6 || isempty(omega_sys)
        omega_sys = 1.0;
    end

    % Allow dynamic load profiles if p_load is provided as a function handle
    if isa(refs.p_load, 'function_handle')
        p_load_val = refs.p_load(t);
    else
        p_load_val = refs.p_load;
    end

    % Core dynamics
    [xdot, I_afe_in, p_pcc] = datacenter_core_equations(x, Vport, dc, p_load_val, omega_sys);

    % Shunt and injection calculation
    % Shunt current absorbed from grid
    I_shunt_in = 1j * refs.B_shunt * Vport;

    % Total current absorbed by the data-center node
    I_total_in = I_afe_in + I_shunt_in;

    % Current injected into the grid network under the Ybus convention
    Iport = -I_total_in;

    aux = struct('p_load', p_load_val, 'p_pcc', p_pcc, 'I_afe', I_afe_in);
end

function [xdot, I_afe_in, p_pcc] = datacenter_core_equations(x, Vport, dc, p_load, omega_sys)
    % Unpack 21-state vector
    theta_pll=x(1); eps_pll=x(2); vq_pll_f=x(3); id_afe=x(4); iq_afe=x(5);
    xi_dc_afe=x(6); gamd_afe=x(7); gamq_afe=x(8); vdc_ups=x(9); iU_cv=x(10);
    iV_cv=x(11); vU_vsi=x(12); vV_vsi=x(13); xiU_vsi=x(14); xiV_vsi=x(15);
    gamU_vsi=x(16); gamV_vsi=x(17); v_psu=x(18); xi_psu=x(19); v_eq=x(20); xi_eq=x(21);

    % Safe divisions to prevent singularity during initial transients
    vdc_ups_safe = max(vdc_ups, 1e-6);
    v_psu_safe   = max(v_psu, 1e-6);
    veq_ref_safe = max(dc.veq_ref, 1e-6);

    % PCC voltage in the global RI frame
    vr_pcc = real(Vport);
    vi_pcc = imag(Vport);

    % Downstream DC-DC/load equivalent
    g_load = p_load / (3 * veq_ref_safe^2);
    i_eq   = dc.kp_v_eq*(dc.veq_ref - v_eq) + dc.ki_v_eq*xi_eq;
    i_psu  = (v_eq / v_psu_safe) * i_eq;

    % PSU reduced equivalent
    g_eq  = dc.kp_v_psu*(dc.vpsu_ref - v_psu) + dc.ki_v_psu*xi_psu;
    iU_vsi = g_eq * vU_vsi;
    iV_vsi = g_eq * vV_vsi;
    vuv_sq = vU_vsi^2 + vV_vsi^2;
    psu_injection_term = ((g_eq - dc.r_psu*g_eq^2) * vuv_sq) / (3 * v_psu_safe);

    % VSI controls
    vU_ref_vsi = dc.vu_vsi_ref;
    vV_ref_vsi = 0;
    iU_ref_cv = dc.kp_v_vsi*(vU_ref_vsi - vU_vsi) + dc.ki_v_vsi*xiU_vsi - dc.omega_vsi*dc.c_vsi*vV_vsi;
    iV_ref_cv = dc.kp_v_vsi*(vV_ref_vsi - vV_vsi) + dc.ki_v_vsi*xiV_vsi + dc.omega_vsi*dc.c_vsi*vU_vsi;
    vU_ref_cv = dc.kp_c_vsi*(iU_ref_cv - iU_cv) + dc.ki_c_vsi*gamU_vsi - dc.omega_vsi*dc.l_vsi*iV_cv;
    vV_ref_cv = dc.kp_c_vsi*(iV_ref_cv - iV_cv) + dc.ki_c_vsi*gamV_vsi + dc.omega_vsi*dc.l_vsi*iU_cv;
    mU = vU_ref_cv / vdc_ups_safe;
    mV = vV_ref_cv / vdc_ups_safe;

    % AFE PLL and projections
    omega_pll = dc.omega_s + dc.kp_pll*vq_pll_f + dc.ki_pll*eps_pll;
    s = sin(theta_pll + pi/2);
    c = cos(theta_pll + pi/2);

    ir_pcc =  s*id_afe + c*iq_afe;
    ii_pcc = -c*id_afe + s*iq_afe;
    I_afe_in = ir_pcc + 1j*ii_pcc; % Complex current into AFE

    vd_pcc = s*vr_pcc - c*vi_pcc;
    vq_pcc = c*vr_pcc + s*vi_pcc;

    % AFE controls
    id_ref_afe = dc.kp_dc_afe*(dc.vdc_ups_ref - vdc_ups) + dc.ki_dc_afe*xi_dc_afe;
    iq_ref_afe = 0;
    vd_ref_afe = dc.kp_c_afe*(id_afe - id_ref_afe) + dc.ki_c_afe*gamd_afe + omega_pll*dc.l_afe*iq_afe;
    vq_ref_afe = dc.kp_c_afe*(iq_afe - iq_ref_afe) + dc.ki_c_afe*gamq_afe - omega_pll*dc.l_afe*id_afe;
    md = vd_ref_afe / vdc_ups_safe;
    mq = vq_ref_afe / vdc_ups_safe;

    % DC-link exchange
    i_dc_in  = md*id_afe + mq*iq_afe;
    i_dc_out = mU*iU_cv + mV*iV_cv;

    % State derivatives
    xdot = zeros(21,1);

    % Use dynamic omega_sys for the relative angle
    xdot(1)  = dc.omega_b * (omega_pll - omega_sys);
    xdot(2)  = vq_pll_f;
    xdot(3)  = dc.omega_lp * (vq_pcc - vq_pll_f);

    xdot(4)  = (dc.omega_b/dc.l_afe) * (vd_pcc - md*vdc_ups - dc.r_afe*id_afe + omega_pll*dc.l_afe*iq_afe);
    xdot(5)  = (dc.omega_b/dc.l_afe) * (vq_pcc - mq*vdc_ups - dc.r_afe*iq_afe - omega_pll*dc.l_afe*id_afe);
    xdot(6)  = (dc.vdc_ups_ref - vdc_ups);
    xdot(7)  = (id_afe - id_ref_afe);
    xdot(8)  = (iq_afe - iq_ref_afe);

    xdot(9)  = (dc.omega_b/dc.c_dc) * (i_dc_in - i_dc_out);

    xdot(10) = (dc.omega_b/dc.l_vsi) * (mU*vdc_ups - vU_vsi - dc.r_vsi*iU_cv + dc.omega_vsi*dc.l_vsi*iV_cv);
    xdot(11) = (dc.omega_b/dc.l_vsi) * (mV*vdc_ups - vV_vsi - dc.r_vsi*iV_cv - dc.omega_vsi*dc.l_vsi*iU_cv);
    xdot(12) = (dc.omega_b/dc.c_vsi) * (iU_cv - iU_vsi + dc.omega_vsi*dc.c_vsi*vV_vsi);
    xdot(13) = (dc.omega_b/dc.c_vsi) * (iV_cv - iV_vsi - dc.omega_vsi*dc.c_vsi*vU_vsi);

    xdot(14) = (vU_ref_vsi - vU_vsi);
    xdot(15) = (vV_ref_vsi - vV_vsi);
    xdot(16) = (iU_ref_cv - iU_cv);
    xdot(17) = (iV_ref_cv - iV_cv);

    xdot(18) = (dc.omega_b/dc.c_psu) * (psu_injection_term - i_psu);
    xdot(19) = (dc.vpsu_ref - v_psu);
    xdot(20) = (dc.omega_b/dc.c_eq) * (i_eq - g_load*v_eq);
    xdot(21) = (dc.veq_ref - v_eq);

    % Power drawn at PCC
    p_pcc = vr_pcc*ir_pcc + vi_pcc*ii_pcc;
end
