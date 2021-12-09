function [u_self1, v_self1, u_self2, v_self2] = func_self_vel_sheet_pair_walls_regul...
    (xe, ye1, ye2, delem, gamma_elem, D_tunnel, delta_core)

% function to compute self-induced velocity of a pair of vortex sheets in a tunnel of certain width,
% separated by a fixed distance
% input of sheet parameters should have already been discretized
% Assumption: 
% sheet # 1 has circulation density of -gamma_elem, 
% and sheet # has circulation density of +gamma_elem
% output: self-induced velocity of both sheets, with components in both x
% and y directions


Nelem = length(delem); % total number of elements

u_self1 = zeros(Nelem, 1);
v_self1 = zeros(Nelem, 1);
u_self2 = zeros(Nelem, 1);
v_self2 = zeros(Nelem, 1);

% iterate over all elements on sheet 1
for Ielem = 1:Nelem % loop over all target elements
    
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    % effect of element i on sheet 1, on the velocity of all points on sheet 1 and 2
    [u_sheet1, v_sheet1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, -Gamma_elem,  xe, ye1, delta_core); % _desingul
    
    [u_sheet2, v_sheet2] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, -Gamma_elem,  xe, ye2, delta_core); % _desingul
    
    u_self1 = u_self1 + u_sheet1(:);
    v_self1 = v_self1 + v_sheet1(:);
    
    u_self2 = u_self2 + u_sheet2(:);
    v_self2 = v_self2 + v_sheet2(:);
    
end % for Ielem

% iterate over all elements on sheet 2
for Ielem = 1:Nelem % loop over all target elements
    
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    % effect of element i on sheet 2, on the velocity of all points on sheet 1 and 2
    [u_sheet1, v_sheet1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye2(Ielem), D_tunnel, Gamma_elem,  xe, ye1, delta_core); % _desingul
    
    [u_sheet2, v_sheet2] = func_vortex_velocity_walls_regularize(xe(Ielem), ye2(Ielem), D_tunnel, Gamma_elem,  xe, ye2, delta_core); % _desingul
    
    u_self1 = u_self1 + u_sheet1(:);
    v_self1 = v_self1 + v_sheet1(:);
    
    u_self2 = u_self2 + u_sheet2(:);
    v_self2 = v_self2 + v_sheet2(:);
    
end % for Ielem
