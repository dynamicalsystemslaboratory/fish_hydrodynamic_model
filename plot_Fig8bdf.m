% plot Fig 8 b d f

clear
close all
clc

load('Data_Fig8ace.mat', 'u_all');
u_CFD = u_all; 

load('Data_Fig8bdf.mat', 'x_all', 'y_all', 'u_all', 'v_all', 'ys0', 'xs_all', 'ys1_all', 'ys2_all', 'D_all', 'U_all');

%% tunnel parameters
L_fish = 0.035481891661866; % fish body length

L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain

x_head = L_tunnel/2; % location of the fish head;

ND = length(D_all);

%% read all data

for ID = 1:ND
    U_swim = U_all(ID); % m/s; inlet speed
    D_tunnel = D_all(ID);
    D_domain = D_tunnel;
    
    %% read simulation data
    x = x_all{ID};
    y = y_all{ID};
    u_opt = u_all{ID};
    v_opt = v_all{ID};
    
    xs = xs_all(ID, :);
    ys1 = ys1_all(ID, :);
    ys2 = ys2_all(ID, :);
    %% plot optimal solution
    
    uv_mag = sqrt(u_opt.^2 + v_opt.^2); % velocity magnitude
    
    %% plot optimal fitted velocity field
    Ix_mask = isnan(u_CFD{ID}); 
    uv_mag(Ix_mask) = nan;
    
    figure(70+ID);
    clf;
    hold on;
    pcolor((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, uv_mag/U_swim);
    shading flat
    axis equal
    
    h1=streamslice((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, u_opt/U_swim, v_opt/U_swim, 2, 'arrows');
    clr_str = [1,1,1];
    set( h1, 'Color', clr_str )
    set(h1, 'linewidth', 1)
    
    plot( ([xs(1), xs(end)]-x_head)/L_tunnel, ([ys1(1), ys1(end)]-ys0)/L_tunnel, 'k:')
    plot( ([xs(1), xs(end)]-x_head)/L_tunnel, ([ys2(1), ys2(end)]-ys0)/L_tunnel, 'k-')
    
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-');
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-');
    
    xlim([-0.5, 1.5]);
    xlabel('$x/L_t$','Interpreter','latex');
    ylabel('$y/L_t$','Interpreter','latex')
    colormap('jet');
    caxis([0, 0.5]); % limits for vel magnitude

end % for ID
