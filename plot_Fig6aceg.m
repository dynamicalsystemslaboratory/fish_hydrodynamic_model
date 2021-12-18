% plot Fig 6 b c e g
clear 
close all

load('Data_Fig6aceg.mat', 'x', 'y', 'u_all', 'v_all'); 

%% fish parameters
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain
D_domain = D_tunnel;

%% all cases
U_all = [0.05, 0.1, 0.25, 0.5]; % m/s; all speeds
NU = length(U_all);

%% read all data
for IU = 1:NU
    U_swim = U_all(IU); % m/s; inlet speed
        
    %% read simulation data

    u_fish_mean = u_all(:, :, IU);
    v_fish_mean = v_all(:, :, IU);
    
    u_fish_mag = sqrt(u_fish_mean.^2 + v_fish_mean.^2); % velocity magnitude
    
    %% plot velocity
    
    figure(20+IU);
    clf;
    hold on;
    pcolor(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mag/U_swim);
    shading flat
    %     cx = colorbar;
    axis equal tight
    
    h1=streamslice(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mean/U_swim, v_fish_mean/U_swim, 1.5, 'arrows');
    clr_str = [1,1,1];
    set( h1, 'Color', clr_str )
    set(h1, 'linewidth', 1)
    
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-');
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-');

    xlabel('$x/L_t$','Interpreter','latex')
    ylabel('$y/L_t$','Interpreter','latex')
    caxis([0, 0.5]); % limits for vel magnitude

    colormap jet
    

end % for IU
