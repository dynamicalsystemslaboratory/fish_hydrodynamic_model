function [u_all, v_all] = func_velocity_sheet_pair_walls_int ...
    (xe, ye1, ye2, D_tunnel, delem, gamma_elem, x, y, delta_core)
% this code is adapted from func_velocity_sheet_pair_int_v2.m,
% modification: velcoity field takes into consideration presence of walls
% through the use of function func_sheet_velocity_walls_desingul
%
%
% input: discretized element information, including element end point
% coordinates, size, dipole orientation
% input also includes target seed points in the fluid domain
%
% output: velocity at the seed locations

%% assign arrays
Nelem = length(xe); % number of elements

% allocate velocity array at seed locations
u_all = zeros(size(x));
v_all = zeros(size(x));

%% iterate over all elements
for Ielem = 1:Nelem
    
    %% circulation of element Ielem
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    %% compute velocity field
    
    [u_f1, v_f1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  x, y, delta_core); % _desingul
    
    [u_f2, v_f2] = func_vortex_velocity_walls_regularize(xe(Ielem), ye2(Ielem), D_tunnel, -Gamma_elem,  x, y, delta_core);
    
    u = u_f1 + u_f2;
    v = v_f1 + v_f2;
    %% combine all dipoles velocity
    u_all = u_all + u;
    v_all = v_all + v;
    
end
