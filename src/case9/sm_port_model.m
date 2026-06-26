function [xdot, Iout, alg] = sm_port_model(t, x, V, sm, pref)
%SM_PORT_MODEL Dynamic port model for the synchronous machine.

% Unpack states
delta=x(1); omega=x(2); eqp=x(3); edp=x(4); efd=x(5); vf=x(6); vrA=x(7); psv=x(8); tm=x(9);

p_ref = pref; if isa(pref,'function_handle'); p_ref = pref(t); end

% dq transformation and stator algebra
vr = real(V); vi = imag(V); vt = hypot(vr,vi);
vd =  sin(delta)*vr - cos(delta)*vi;
vq =  cos(delta)*vr + sin(delta)*vi;

id = (eqp - vq)/sm.xdp;
iq = (vd  - edp)/sm.xqp;

ir =  sin(delta)*id + cos(delta)*iq;
ii = -cos(delta)*id + sin(delta)*iq;
Iout = ir + 1j*ii;

pe = vd*id + vq*iq;
se = sm.ae * exp(sm.be*efd);

% ODE assembly
xdot = zeros(9,1);
xdot(1) = sm.omega_b*(omega - sm.ws);
xdot(2) = (tm - pe - sm.D*(omega - sm.ws)) / (2*sm.H);
xdot(3) = (-eqp - (sm.xd-sm.xdp)*id + efd) / sm.Td0p;
xdot(4) = (-edp + (sm.xq-sm.xqp)*iq) / sm.Tq0p;
xdot(5) = (-(sm.ke+se)*efd + vrA) / sm.Te;
xdot(6) = (-vf + (sm.kf/sm.Te)*vrA - (sm.kf/sm.Te)*(sm.ke+se)*efd) / sm.Tf;
xdot(7) = (-vrA + sm.ka*(sm.vref - vf - vt)) / sm.Ta;
xdot(8) = (-psv + p_ref - (1/sm.r)*(omega - sm.ws)) / sm.Tsv;
xdot(9) = (-tm + psv) / sm.Tch;

alg = struct('vt',vt, 'pe',pe, 'qe',vq*id - vd*iq, 'id',id, 'iq',iq, 'tm',tm);
end
