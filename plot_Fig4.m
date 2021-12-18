% plot figure 4

clear
close all
clc

FTsz = 20;
LNwd = 2;

%% load CFD data

load('Data_Fig4.mat', 'x', 'y', 'u_all', 'v_all'); 

%% fish parameters
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain
D_domain = D_tunnel;

U_swim = 0.156379646024586; % m/s; fish swimming speed
%% plot instantaneous flow
Icount = 0;

for Iframe = 1:6
    Icount = Icount + 1;
    
    u_fish = u_all(:,:,Icount);
    v_fish = v_all(:,:,Icount); 
    
    u_fish_mag = sqrt(u_fish.^2 + v_fish.^2); % velocity magnitude
    
    %% plot the case
    figure(50+Icount);
    clf
    hold on;
    pcolor(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mag/U_swim);
    shading flat
    axis equal tight
    h1=streamslice(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish/U_swim, v_fish/U_swim, 3, 'arrows');
    clr_str = [1,1,1];
    set( h1, 'Color', clr_str )
    set(h1, 'linewidth', 1)
    
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-');
    plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-');

    xlabel('$x/L_t$','Interpreter','latex')
    ylabel('$y/L_t$','Interpreter','latex')
    axis([-0.5, 1.5, 0, 0.31]);
    set(gca, 'YTick', 0:0.15:0.3);
    
    colormap jet
    
end
