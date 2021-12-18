% plot Fig 5 f g

load('Data_Fig5fg.mat', 'x', 'y', 'x_head', 'ys0', 'ys1', 'ys2', 'xs', 'u_opt', 'v_opt', 'u_fish_mean', 'v_fish_mean');

%% vortex sheet parameters
L_fish = 0.035481891661866; % fish body length
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain
D_domain = D_tunnel;

U_swim = 0.156379646024586; % m/s; fish swimming speed

%% plot optimal solution

uv_mag = sqrt(u_opt.^2 + v_opt.^2); % velocity magnitude

diff_uv = (u_opt - u_fish_mean).^2 +  (v_opt - v_fish_mean).^2;

Ix_mask = isnan(u_fish_mean);
uv_mag(Ix_mask) = nan;
diff_uv(Ix_mask) = nan;
%% plot optimal fitted velocity field

figure(59);
clf;
hold on;
pcolor((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, uv_mag/U_swim);
shading flat
h1=streamslice((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, u_opt/U_swim, v_opt/U_swim, 2, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

plot( (xs-x_head)/L_tunnel, (ys1-ys0)/L_tunnel, 'k:', 'linewidth',LNwd)
plot( (xs-x_head)/L_tunnel, (ys2-ys0)/L_tunnel, 'k-', 'linewidth',LNwd)

plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-', 'linewidth', LNwd );
plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-', 'linewidth', LNwd );

xlabel('$x/L_t$','Interpreter','latex')
ylabel('$y/L_t$','Interpreter','latex')

colormap('jet');
caxis([0, 0.5]); % limits for vel magnitude
axis equal


figure(60);
clf;
hold on;
pcolor((x-x_head)/L_tunnel, (y-ys0)/L_tunnel, sqrt(diff_uv)/U_swim);
shading flat

plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [0, 0]/L_tunnel, 'k-', 'linewidth', LNwd );
plot([-L_tunnel/2, 1.5*L_tunnel]/L_tunnel, [D_domain, D_domain]/L_tunnel, 'k-', 'linewidth', LNwd );

xlabel('$x/L_t$','Interpreter','latex')
ylabel('$y/L_t$','Interpreter','latex')
colormap jet
caxis([0, 0.5]); % limits for vel magnitude

axis equal
