function [u, v] = func_vortex_velocity_walls_regularize(xf, yf, H, Gamma, x, y, delta_core)
% velocity field induced by a pair of vortex sheet between 2 parallel walls
% xf, yf, theta: fish position and orientation with respect to x axis
% x, y: target location in space
% H: distance between 2 walls
% Gamma: circulation of dipole (a pair of vortex sheet with infinitesimal length)
% delta_core: core radius of vortex filament; used to desingularize point
% vortex


k = pi/H; % constant k

yim = -yf;

% define common terms
denom1 = cosh(k*(x-xf)) - cos(k*(y-yf)); 
denom2 = cosh(k*(x-xf)) - cos(k*(y-yim)); 

sinhkx = sinh(k*(x-xf));

% compute velocity components

u = -(0.25*Gamma/pi)*k*(   sin(k*(  y - yf ))./( denom1 + 0.5*k^2*delta_core^2 )   -   sin(k*(y-yim))./denom2   ); % velocity along x 

v = (0.25*Gamma/pi)*k*(   sinhkx./( denom1 + 0.5*k^2*delta_core^2 )   -   sinhkx./denom2   ); % velocity along y
