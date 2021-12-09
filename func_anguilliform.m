function z = func_anguilliform(x, t, c1, c2, c3, c4, kL, omega, y_shift, phi_0)

% phi0 = 0; % no phase shift


% define y function as piecewise function of x

z = (c1 + c2*x + c3*x.^2 + c4*x.^3 ).*sin(kL*x - omega*t + phi_0) + y_shift; % tail region
