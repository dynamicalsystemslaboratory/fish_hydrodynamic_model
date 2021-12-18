% plot Fig 6 bdfh

clear
close all
clc

load('Data_Fig6aceg.mat', 'u_all'); 
u_CFD = u_all; 

load('Data_Fig6bdfh.mat', 'x', 'y', 'u_all', 'v_all', 'ys0', 'xs_all', 'ys1_all', 'ys2_all');

%% tunnel parameters
L_fish = 0.035481891661866; % fish body length
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain
D_domain = D_tunnel;

x_head = L_tunnel/2; % location of the fish head;

%% all swim speed cases
U_all = [0.05, 0.1, 0.25, 0.5]; % m/s; all speeds
NU = length(U_all);

%% read all data
for IU = 1:NU
    U_swim = U_all(IU); % m/s; inlet speed
    
    u_opt = u_all(:, :, IU);
    v_opt = v_all(:, :, IU);
    xs = xs_all(IU, :);
    ys1 = ys1_all(IU, :);
    ys2 = ys2_all(IU, :);
    
    u_simul = u_CFD(:, :, IU);
    %% plot optimal solution
    uv_mag = sqrt(u_opt.^2 + v_opt.^2); % velocity magnitude
    
    %% plot optimal fitted velocity field
    Ix_mask = isnan(u_simul);
    uv_mag(Ix_mask) = nan;
    
    figure(70+IU);
    clf;
    hold on;
    pcolor((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, uv_mag/U_swim);
    shading flat
    h1=streamslice((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, u_opt/U_swim, v_opt/U_swim, 2, 'arrows');
    clr_str = [1,1,1];
    set( h1, 'Color', clr_str )
    set(h1, 'linewidth', 1)
    
    plot( ([xs(1), xs(end)]-x_head)/L_tunnel, ([ys1(1), ys1(end)]-ys0)/L_tunnel, 'k:')
    plot( ([xs(1), xs(end)]-x_head)/L_tunnel, ([ys2(1), ys2(end)]-ys0)/L_tunnel, 'k-')
    
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-' );
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-' );
    
    xlabel('$x/L_t$','Interpreter','latex');
    ylabel('$y/L_t$','Interpreter','latex')
    
    colormap('jet');
    caxis([0, 0.5]); % limits for vel magnitude
    axis equal tight

end % for IU
