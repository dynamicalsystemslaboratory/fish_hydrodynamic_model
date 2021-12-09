function [u, v] = func_dipole_velocity_fs(xf, yf, theta, r0, v0, x, y)
% velocity field induced by a dipole in free space
% xf, yf, theta: fish position and orientation with respect to x axis
% x, y: target location in space
% r0: characteristic dipole length-scale (on the order of the amplitude of the fish tail beating)
% v0: constant speed of the fish

d4 = ( (x - xf).^2 + (y - yf).^2 ).^2; % d^4, 

d4(d4<1e-10) = 1e-10; % if target point is too close to fish location, substitute with a small number to avoid singularity

u = r0^2 * v0* ( ( (x - xf).^2 - (y - yf).^2 ).*cos(theta)  + 2*(x - xf).*(y - yf).*sin(theta) )./d4; % velocity along x 

v = r0^2 * v0* ( -( (x - xf).^2 - (y - yf).^2 ).*sin(theta)  + 2*(x - xf).*(y - yf).*cos(theta) )./d4;  % velocity along y 