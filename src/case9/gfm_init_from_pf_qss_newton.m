function [x0, refs] = gfm_init_from_pf_qss_newton(Vport_pf, Iport_pf, gfm)
%GFM_INIT_FROM_PF_QSS_NEWTON Initialize a grid-forming inverter from power-flow results.

b = gfm.base; f = gfm.filt; v = gfm.virt;

% Complex QSS calculation
Zg = f.rg + 1j*(b.ome_sys * f.lg);
Vf = Vport_pf + Zg * Iport_pf;

Zv = v.rv + 1j*(b.ome_sys * v.lv);
Vs = Vf + Zv * Iport_pf;

% Power and capacitor injections
S_filt = Vf * conj(Iport_pf);
pref0 = real(S_filt); qref0 = imag(S_filt);

Ic = 1j * b.ome_sys * f.cf * Vf;
Icv = Iport_pf + Ic;
Vcv = Vf + (f.rf + 1j*b.ome_sys*f.lf) * Icv;

refs = struct('pref0', pref0, 'qref0', qref0, 'ome_ref0', 1.0, 'v_ref0', abs(Vs), ...
              'Vf_qss', Vf, 'Vcv_qss', Vcv);

% Initial guess and solve
% x = [ir_cv, ii_cv, vr_f, vi_f, ir_f, ii_f, the, pm, qm, xi_d, xi_q, gam_d, gam_q]
xg = [real(Icv); imag(Icv); real(Vf); imag(Vf); real(Iport_pf); imag(Iport_pf); ...
      angle(Vcv); pref0; qref0; 0; 0; 0; 0];

typX = max(abs(xg), 1e-2); typX(7:9) = 1.0;

opt = optimoptions('fsolve', 'Display', 'none', 'Algorithm', 'levenberg-marquardt', ...
                   'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'TypicalX', typX);

fun = @(x) gfm_port_model(0, x, Vport_pf, gfm, pref0, qref0, 1.0, abs(Vs));
[x0, Fsol, exitflag] = fsolve(fun, xg, opt);

if exitflag <= 0
    [mF, idx] = max(abs(Fsol));
    warning('GFM Newton did not fully converge (exitflag=%d). max|F|=%.3e at eqn #%d', exitflag, mF, idx);
end
end
