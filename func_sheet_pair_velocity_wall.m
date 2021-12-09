function [u, v] = func_sheet_pair_velocity_wall(xf, yf, theta, r0, Gamma, h, x, y)
% velocity field induced by a dipole in free space
% xf, yf, theta: fish position and orientation with respect to x axis
% x, y: target location in space
% r0: characteristic dipole length-scale (on the order of the amplitude of the fish tail beating)
% Gamma: circulation of dipole (a pair of vortex sheet with infinitesimal length)
% h: tunnel width

d2 = ( (x - xf).^2 + (y - yf).^2 ); % d^4, 

% d2(d2<1e-15) = 1e-15; % if target point is too close to fish location, substitute with a small number to avoid singularity

A = ( (x - xf) + (y - yf)*1i) / (2*h); 
B = ( (x - xf) + (y + yf)*1i) / (2*h); 

u = 0.25 * (0.5/pi)*Gamma*r0* (        (pi^2*exp(-1i*theta)/(2*h^2)) * (  exp(2i*theta)*( (csch(pi*A)).^2 + ( csch(pi*conj(B)) ).^2 ) +  ( csch(pi*conj(A)) ).^2 + (csch(pi*B)).^2      ) ...
                                        + 4*cos(theta)./d2 -  8*(x - xf).*(  (x - xf).*cos(theta) +  (y - yf).*sin(theta) )./(d2.^2)            );  % velocity along x 

v = 0.25 * (0.5/pi)*Gamma*r0* (        (1i*pi^2*exp(-1i*theta)/(2*h^2)) * (  exp(2i*theta)*( (csch(pi*A)).^2 - ( csch(pi*conj(B)) ).^2 ) -  ( csch(pi*conj(A)) ).^2 + (csch(pi*B)).^2      ) ...
                                        + 4*sin(theta)./d2 -  8*(y - yf).*(  (x - xf).*cos(theta) +  (y - yf).*sin(theta) )./(d2.^2)            ); % velocity along y 
                      
% if any(real(u)>1,'all')
%     pause(0.1);
% end