function [xdot, Iport, aux] = gfm_port_model(t, x, Vport, gfm, pref_cmd, qref_cmd, ome_ref, v_ref, omega_sys)
%GFM_PORT_MODEL Dynamic port model for a grid-forming inverter.

% Use nominal frequency when omega_sys is not provided
if nargin < 9 || isempty(omega_sys)
    omega_sys = gfm.base.ome_sys;
end

b = gfm.base; f = gfm.filt; dr = gfm.droop; vi = gfm.virt; in = gfm.inner; dc = gfm.dc;

% Unpack state vector
ir_cv=x(1); ii_cv=x(2); vr_f=x(3); vi_f=x(4); ir_f=x(5); ii_f=x(6);
the=x(7); pm=x(8); qm=x(9); xi_d=x(10); xi_q=x(11); gam_d=x(12); gam_q=x(13);

% dq transformation
sn = sin(the + pi/2); cs = cos(the + pi/2);
T = [sn, -cs; cs, sn];
Vdq_f = T*[vr_f; vi_f]; Idq_f = T*[ir_f; ii_f]; Idq_cv = T*[ir_cv; ii_cv];
vd_f=Vdq_f(1); vq_f=Vdq_f(2); id_f=Idq_f(1); iq_f=Idq_f(2); id_cv=Idq_cv(1); iq_cv=Idq_cv(2);

% Droop and virtual impedance
ome_oc = ome_ref + dr.Kp*(pref_cmd - pm);
v_oc   = v_ref   + dr.Kq*(qref_cmd - qm);

vd_vi_ref = v_oc - vi.rv*id_f + ome_oc*vi.lv*iq_f;
vq_vi_ref = -vi.rv*iq_f - ome_oc*vi.lv*id_f;

% PI loops and modulation
id_cv_ref = in.kpv*(vd_vi_ref - vd_f) + in.kiv*xi_d - f.cf*ome_oc*vq_f;
iq_cv_ref = in.kpv*(vq_vi_ref - vq_f) + in.kiv*xi_q + f.cf*ome_oc*vd_f;

vd_cv_ref = in.kpc*(id_cv_ref - id_cv) + in.kic*gam_d - ome_oc*f.lf*iq_cv;
vq_cv_ref = in.kpc*(iq_cv_ref - iq_cv) + in.kic*gam_q + ome_oc*f.lf*id_cv;

md = vd_cv_ref / dc.Vdc; mq = vq_cv_ref / dc.Vdc;
vr_cv = ( sn*md + cs*mq ) * dc.Vdc; vi_cv = (-cs*md + sn*mq ) * dc.Vdc;

% Dynamics assembly (b.ome_sys is replaced by dynamic omega_sys)
xdot = zeros(13,1);
xdot(1:2) = (b.ome_b/f.lf) * [vr_cv - vr_f - f.rf*ir_cv + omega_sys*f.lf*ii_cv;
                              vi_cv - vi_f - f.rf*ii_cv - omega_sys*f.lf*ir_cv];
xdot(3:4) = (b.ome_b/f.cf) * [ir_cv - ir_f + omega_sys*f.cf*vi_f;
                              ii_cv - ii_f - omega_sys*f.cf*vr_f];
xdot(5:6) = (b.ome_b/f.lg) * [vr_f - real(Vport) - f.rg*ir_f + omega_sys*f.lg*ii_f;
                              vi_f - imag(Vport) - f.rg*ii_f - omega_sys*f.lg*ir_f];
xdot(7)   = b.ome_b*(ome_oc - omega_sys); % Angle dynamics relative to SM
xdot(8:9) = [dr.ome_z*((vr_f*ir_f + vi_f*ii_f) - pm); dr.ome_f*((-vr_f*ii_f + vi_f*ir_f) - qm)];
xdot(10:13)= [vd_vi_ref - vd_f; vq_vi_ref - vq_f; id_cv_ref - id_cv; iq_cv_ref - iq_cv];

Iport = ir_f + 1j*ii_f;
aux = struct('ome_oc',ome_oc, 'v_oc',v_oc, 'p_inst',pm, 'q_inst',qm);
end
