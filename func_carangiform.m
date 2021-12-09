function z = func_carangiform(x, x_rgd, t, c1, c2, c3, kL, omega, y_shift, phi_0, phi_1)

% phi0 = 0; % no phase shift


% define y function as piecewise function of x

z = (c1*(x-x_rgd)+c2*((x-x_rgd).^2)).*sin(kL*(x-x_rgd)-omega*t+phi_0) ...
    - c1*(x-x_rgd).*sin(-omega*t+phi_0) + c3*sin(-omega*t+phi_0+phi_1) + y_shift; % tail region

z(x<x_rgd) = c3*sin(-omega*t+phi_0+phi_1) + y_shift; % rigid heading area