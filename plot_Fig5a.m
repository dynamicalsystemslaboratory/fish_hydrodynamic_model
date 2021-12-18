% plot Fig5a

clear
close all
clc

FTsz = 20;
LNwd = 2;

load('Data_Fig5a.mat', 'x', 'y', 'u_fish_mean', 'v_fish_mean'); 
%% fish parameters
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain
D_domain = D_tunnel;

U_swim = 0.156379646024586; % m/s; fish swimming speed

u_fish_mag = sqrt(u_fish_mean.^2 + v_fish_mean.^2);

%% plot velocity with tunnel flow

figure(20);
clf;
hold on;
pcolor(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mag/U_swim);
shading flat
axis equal tight

h1=streamslice(x/L_tunnel, (y+D_domain/2)/L_tunnel, u_fish_mean/U_swim, v_fish_mean/U_swim, 1.5, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-', 'linewidth', LNwd );
plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-', 'linewidth', LNwd );

xlabel('$x/L_t$','Interpreter','latex')
ylabel('$y/L_t$','Interpreter','latex')

caxis([0, 0.5]); % limits for vel magnitude
axis([-0.5, 1.5, 0, 0.31]);

colormap jet
