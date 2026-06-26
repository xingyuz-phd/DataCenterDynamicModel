function [x0, sm, refs, info] = sm_init_from_pf_fsolve(Vpf, Ipf, sm, opts)
%SM_INIT_FROM_PF_FSOLVE Initialize the synchronous-machine port model from power-flow results.

if nargin < 4 || isempty(opts)
    opts = optimoptions('fsolve','Display','none', 'FunctionTolerance',1e-12, ...
                        'StepTolerance',1e-12, 'MaxIterations',200);
end
if ~isfield(sm,'ws'); sm.ws = 1.0; end

pref0 = real(Vpf * conj(Ipf));
vt0   = abs(Vpf);

% Find delta via consistency check (inline logic)
vr = real(Vpf); vi = imag(Vpf); ir = real(Ipf); ii = imag(Ipf);
delta_guess = angle(Vpf);

try
    consist_eq = @(d) (sin(d)*vr - cos(d)*vi) - sm.xqp*(cos(d)*ir + sin(d)*ii) - ...
                      (sm.xq - sm.xqp)*(cos(d)*ir + sin(d)*ii);
    delta_guess = fsolve(consist_eq, delta_guess, optimoptions(opts,'MaxIterations',50));
catch
    % keep fallback if it fails
end

% Extract DQ and compute steady states
vd = sin(delta_guess)*vr - cos(delta_guess)*vi;
vq = cos(delta_guess)*vr + sin(delta_guess)*vi;
id = sin(delta_guess)*ir - cos(delta_guess)*ii;
iq = cos(delta_guess)*ir + sin(delta_guess)*ii;

eqp0 = vq + sm.xdp*id;
edp0 = 0.5 * ((vd - sm.xqp*iq) + (sm.xq-sm.xqp)*iq);
efd0 = max(1e-6, eqp0 + (sm.xd - sm.xdp)*id);
vrA0 = (sm.ke + sm.ae*exp(sm.be*efd0)) * efd0;
vref0 = vt0 + vrA0/sm.ka;

x_init = [delta_guess; sm.ws; eqp0; edp0; efd0; 0; vrA0; pref0; pref0];

% Update references and Newton solve
sm.vref = vref0; sm.pref = pref0;
refs = struct('vref',vref0, 'pref',pref0, 'ws',sm.ws);

fun = @(x) sm_port_model(0, x, Vpf, sm, pref0);
[x0, fval, exitflag, output] = fsolve(fun, x_init, opts);

info = struct('exitflag', exitflag, 'output', output, 'residual_inf', norm(fval, inf));
end
