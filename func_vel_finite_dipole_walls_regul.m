function [u_all, v_all] = func_vel_finite_dipole_walls_regul ...
    (xs0, ys1, ys2, gamma_point, D_tunnel, x, y, delta_core)
% this code computes velocity field of a finite dipole
% 
% input: two vortices locations
%
% Assumption: 
% both vortices have the same x-coordinate
% vortex # 1 has circulation density of -gamma_elem, 
% and vortex # 2 has circulation density of +gamma_elem
%
% output: velocity at the seed locations


    %%  velocity induced by vortex 1
    [u_f1, v_f1] = func_vortex_velocity_walls_regularize(xs0, ys1, D_tunnel, -gamma_point,  x, y, delta_core); % _desingul
    
    %%  velocity induced by vortex 2
    [u_f2, v_f2] = func_vortex_velocity_walls_regularize(xs0, ys2, D_tunnel, gamma_point,  x, y, delta_core);
    
    %% combined
    u_all = u_f1 + u_f2;
    v_all = v_f1 + v_f2;
