function [xdot, Iport, aux] = gfl_port_model(t, x, Vport, gfl, pref_cmd, qref_cmd, omega_sys)
%GFL_PORT_MODEL Dynamic port model for a grid-following inverter.

% Use nominal frequency during initialization
if nargin < 7 || isempty(omega_sys)
    omega_sys = gfl.base.ome_sys;
end

b = gfl.base; f = gfl.filt; pll = gfl.pll; out = gfl.outer; in = gfl.inner; dc = gfl.dc;

% Unpack state vector
ir_cv=x(1); ii_cv=x(2); vr_filt=x(3); vi_filt=x(4); ir_g=x(5); ii_g=x(6);
vq_pll=x(7); eps_pll=x(8); the_pll=x(9);
sig_p=x(10); pm=x(11); sig_q=x(12); qm=x(13); gam_d=x(14); gam_q=x(15);

% PLL and dq projections
ome_est = 1 + pll.kp*vq_pll + pll.ki*eps_pll;
s = sin(the_pll + pi/2); c = cos(the_pll + pi/2);

vd_f = s*vr_filt - c*vi_filt; vq_f = c*vr_filt + s*vi_filt;
id_g = s*ir_g - c*ii_g;       iq_g = c*ir_g + s*ii_g;
id_cv= s*ir_cv - c*ii_cv;     iq_cv= c*ir_cv + s*ii_cv;

% Outer loop (P/Q)
p_inst = vr_filt*ir_g + vi_filt*ii_g;
q_inst = vi_filt*ir_g - vr_filt*ii_g;

% Map active power to the d-axis and reactive power to the q-axis
Id_pi =  (out.Kp_p*(pref_cmd - pm) + out.Ki_p*sig_p);
Iq_pi = -(out.Kp_q*(qref_cmd - qm) + out.Ki_q*sig_q);

% Inner loop (current)
Vd_pi = in.kpc*(Id_pi - id_cv) + in.kic*gam_d;
Vq_pi = in.kpc*(Iq_pi - iq_cv) + in.kic*gam_q;

Vd_cv_ref = Vd_pi - ome_est*in.lf*iq_cv + in.kffv*vd_f;
Vq_cv_ref = Vq_pi + ome_est*in.lf*id_cv + in.kffv*vq_f;

md = Vd_cv_ref / dc.Vdc; mq = Vq_cv_ref / dc.Vdc;
vr_cv = ( s*md + c*mq ) * dc.Vdc; vi_cv = (-c*md + s*mq ) * dc.Vdc;

% LCL and grid dynamics in the RI frame
Vr_g = real(Vport); Vi_g = imag(Vport);

xdot = zeros(15,1);
xdot(1:2) = (b.ome_b/f.lf) * [vr_cv - vr_filt - f.rf*ir_cv + omega_sys*f.lf*ii_cv;
                              vi_cv - vi_filt - f.rf*ii_cv - omega_sys*f.lf*ir_cv];
xdot(3:4) = (b.ome_b/f.cf) * [ir_cv - ir_g + omega_sys*f.cf*vi_filt;
                              ii_cv - ii_g - omega_sys*f.cf*vr_filt];
xdot(5:6) = (b.ome_b/f.lg) * [vr_filt - Vr_g - f.rg*ir_g + omega_sys*f.lg*ii_g;
                              vi_filt - Vi_g - f.rg*ii_g - omega_sys*f.lg*ir_g];
xdot(7:9) = [pll.ome_lp*(vq_f - vq_pll); vq_pll; b.ome_b*(ome_est - omega_sys)];
xdot(10:13)= [pref_cmd - pm; out.ome_z*(p_inst - pm); qref_cmd - qm; out.ome_f*(q_inst - qm)];
xdot(14:15)= [Id_pi - id_cv; Iq_pi - iq_cv];

Iport = ir_g + 1j*ii_g;
aux = struct('ome_est',ome_est, 'p_inst',p_inst, 'q_inst',q_inst, 'Id_pi',Id_pi, 'Iq_pi',Iq_pi);
end
