function [u_self1, v_self1, u_self2, v_self2] = func_self_vel_finite_dipole_walls_regul...
    (xs0, ys1, ys2,  gamma_point, D_tunnel, delta_core)

% function to compute self-induced velocity of a pair of point vortices in a tunnel of certain width,
% separated by a fixed, finite distance (finite dipole model)
% Assumption:
% point # 1 has circulation density of -gamma_elem,
% and point # 2 has circulation density of +gamma_elem
% output: self-induced velocity of both point vortices, with components in both x
% and y directions

%% velocity of point # 1
[u_f11, v_f11] = func_vortex_velocity_walls_regularize(xs0, ys1, D_tunnel, -gamma_point,  xs0, ys1, delta_core); % 1->1

[u_f21, v_f21] = func_vortex_velocity_walls_regularize(xs0, ys2, D_tunnel, gamma_point,  xs0, ys1, delta_core); % 2->1

u_self1 = u_f11 + u_f21;
v_self1 = v_f11 + v_f21;

%% velocity of point # 1
[u_f12, v_f12] = func_vortex_velocity_walls_regularize(xs0, ys1, D_tunnel, -gamma_point,  xs0, ys2, delta_core); % 1->2

[u_f22, v_f22] = func_vortex_velocity_walls_regularize(xs0, ys2, D_tunnel, gamma_point,  xs0, ys2, delta_core); % 2->2

u_self2 = u_f12 + u_f22;
v_self2 = v_f12 + v_f22;
