function [u, v] = func_vortex_velocity_fs_regularize(xf, yf, Gamma, x, y, delta_core)
% velocity field induced by a pair of vortex sheet in free space
% xf, yf, theta: fish position and orientation with respect to x axis
% x, y: target location in space
% Gamma: circulation of dipole (a pair of vortex sheet with infinitesimal length)
% delta_core: core radius of vortex filament; used to desingularize point
% vortex

d2 =  (x - xf).^2 + (y - yf).^2 ; % d^4, 

u = -(0.5*Gamma/pi)*(  y - yf )./( d2 + delta_core^2 ) ; % velocity along x 

v = (0.5*Gamma/pi)*(  x - xf )./( d2 + delta_core^2 ) ; % velocity along y
