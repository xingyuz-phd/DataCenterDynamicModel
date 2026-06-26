function [x0, gfl, refs, info] = gfl_init_from_pf_qss_newton(Vport_pf, Iport_pf, gfl)
%GFL_INIT_FROM_PF_QSS_NEWTON Initialize a grid-following inverter from power-flow results.

base = gfl.base; filt = gfl.filt;

% Step 1: QSS on grid-side inductor
Vr_grid = real(Vport_pf); Vi_grid = imag(Vport_pf);
Ir_g = real(Iport_pf);    Ii_g = imag(Iport_pf);

% Filter voltage from steady-state Lg/Rg
Vr_filt = Vr_grid + filt.rg*Ir_g - base.ome_sys*filt.lg*Ii_g;
Vi_filt = Vi_grid + filt.rg*Ii_g + base.ome_sys*filt.lg*Ir_g;

% Step 2: Setpoints from QSS
Pref0 = Vr_filt*Ir_g + Vi_filt*Ii_g;
Qref0 = Vi_filt*Ir_g - Vr_filt*Ii_g;

gfl.wref0 = 1.0; gfl.pref0 = Pref0; gfl.qref0 = Qref0; gfl.vref0 = hypot(Vr_filt, Vi_filt);
refs = struct('wref0', 1.0, 'vref0', gfl.vref0, 'pref0', Pref0, 'qref0', Qref0);

% Step 3: Initial guess and Newton solve
x_guess = [Ir_g; Ii_g; Vr_filt; Vi_filt; Ir_g; Ii_g; ...
           0; 0; atan2(Vi_filt, Vr_filt); ...
           0; Pref0; 0; Qref0; 0; 0]; % 15x1 state vector

opt = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-12, ...
                   'StepTolerance', 1e-12, 'MaxIterations', 200);

fun = @(x) gfl_port_model(0, x, Vport_pf, gfl, Pref0, Qref0);
[x0, fval, exitflag, output] = fsolve(fun, x_guess, opt);

info = struct('exitflag', exitflag, 'output', output, 'fval_inf', norm(fval, inf));
if exitflag <= 0
    warning('GFL init Newton did not fully converge (exitflag=%d). ||f||_inf=%.3e', exitflag, info.fval_inf);
end
end
