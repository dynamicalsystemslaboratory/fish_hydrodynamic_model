function [u, v] = func_dipole_velocity_wall(xf, yf, theta, r0, v0, h, x, y)
% velocity field induced by a dipole in free space
% xf, yf, theta: fish position and orientation with respect to x axis
% x, y: target location in space
% r0: characteristic dipole length-scale (on the order of the amplitude of the fish tail beating)
% v0: constant speed of the fish
% h: tunnel width

d2 = ( (x - xf).^2 + (y - yf).^2 ); % d^4, 

% d2(d2<1e-10) = 1e-10; % if target point is too close to fish location, substitute with a small number to avoid singularity

A = ( (x - xf) + (y - yf)*1i) / (2*h); 
B = ( (x - xf) + (y + yf)*1i) / (2*h); 

u = 0.25 * (r0^2) * v0* (        (pi^2*exp(-1i*theta)/(2*h^2)) * (  exp(2i*theta)*( (csch(pi*A)).^2 + ( csch(pi*conj(B)) ).^2 ) +  ( csch(pi*conj(A)) ).^2 + (csch(pi*B)).^2      ) ...
                                        + 4*cos(theta)./d2 -  8*(x - xf).*(  (x - xf).*cos(theta) +  (y - yf).*sin(theta) )./(d2.^2)            );  % velocity along x 

v = 0.25 * (r0^2) * v0* (        (1i*pi^2*exp(-1i*theta)/(2*h^2)) * (  exp(2i*theta)*( (csch(pi*A)).^2 - ( csch(pi*conj(B)) ).^2 ) -  ( csch(pi*conj(A)) ).^2 + (csch(pi*B)).^2      ) ...
                                        + 4*sin(theta)./d2 -  8*(y - yf).*(  (x - xf).*cos(theta) +  (y - yf).*sin(theta) )./(d2.^2)            ); % velocity along y 
                                    