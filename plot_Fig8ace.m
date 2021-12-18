% plot Fig 8 a c e

clear
close all
clc

load('Data_Fig8ace.mat', 'x_all', 'y_all', 'u_all', 'v_all', 'D_all', 'U_all');

%% fish parameters
L_fish = 0.035481891661866; % fish body length
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain

%% all cases
ND = length(U_all);

%% read all data

for ID = 1:ND
    U_swim = U_all(ID); % m/s; inlet speed
    D_tunnel = D_all(ID);
    D_domain = D_tunnel;
    
    %% read simulation data
    x = x_all{ID};
    y = y_all{ID};
    u_fish_mean = u_all{ID};
    v_fish_mean = v_all{ID};
    
    u_fish_mag = sqrt(u_fish_mean.^2 + v_fish_mean.^2); % velocity magnitude
    
    %% plot velocity
    
    figure(20+ID);
    clf;
    hold on;
    pcolor(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mag/U_swim);
    shading flat
    axis equal 
    
    h1=streamslice(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mean/U_swim, v_fish_mean/U_swim, 1.5, 'arrows');
    clr_str = [1,1,1];
    set( h1, 'Color', clr_str )
    set(h1, 'linewidth', 1)
    
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-');
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-');
    xlim([-0.5, 1.5]);
    
    xlabel('$x/L_t$','Interpreter','latex')
    ylabel('$y/L_t$','Interpreter','latex')
    caxis([0, 0.5]); 
    
    colormap jet

end % for ID
