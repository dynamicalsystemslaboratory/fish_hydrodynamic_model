% this code is modified from main_fit_sheet_simul_body_v10.m
% flow field is fitted on the side of the body using one dipole

clear
close all
clc

flag_save = 1;

FTsz = 20;
LNwd = 2;
MKsz = 10;

%% directories
save_dir = '../analysis_results_dipole/';
simul_results_dir = '../../fish_swim/analysis_results/';

%% read simulation results
fname = [simul_results_dir, 'simul_vel_opt.mat'];

load(fname, 'x', 'y', 'u_fish_mean', 'v_fish_mean');

x_simul = x-min(x(1,:));
y_simul = y-min(y(:,1));

%% vortex sheet parameters
L_fish = 0.035481891661866; % fish body length
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain

theta = 1*pi; % orientation of sheet pair; angle with x-axis

TBA = (0.002997180904265 + 0.020600558880515 -0.040537635755492 + 0.097600742072607 )*L_fish; % m; corresponding to tail beat amplitude of 4 mm
v0 = 0.156379646024586; % m/s; fish fish swimming speed at lowest flow rate
Gamma = 2*pi*TBA*v0; % m^2/s; treated circulation of one vortex in a dipole; defined in Porfiri PRL 2021
r0 = 2*TBA;

% location of the starting point of vortex sheet pair
x_head = L_tunnel/2; % location of the fish head;
ys0 = D_tunnel/2; %
% xs0 = 0; % 1st dipole position

%% define parameter space
gamma0 = 1e-2; % dimensional circulation

% starting point of the sheet pair
Nxs = 40+1;
% xs_all = linspace(x_head-0.2*L_fish, x_head+1.2*L_fish, Nxs);
xs_all = linspace(0.08, 0.085, Nxs);

% circulation density
Nr0 = 40+1;
r0_all = linspace(1.8e-3, 2.1e-3, Nr0);
% r0_all(1) = 0.5*(r0_all(1)+r0_all(2));  % remove trivial case

%% fluid domain seed points
x = x_simul; % use same seeding points for simulation
y = y_simul;

dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);

%% iterate through parameter space
error_uv = nan(Nxs, Nr0); % assign error function

tic

for Ixs = 1:Nxs
    xs0 = xs_all(Ixs);
    
    fprintf('Computing case xs0: %d out of %d ... \n', Ixs, Nxs );
    
    for Ir0 = 1:Nr0
        
        %% vortex sheet parameters
        r0 = r0_all(Ir0); % (constant) curculation density (per unit length) of one vortex sheet
        
        %% compute velocity field
        [u_fs, v_fs] = func_dipole_velocity_fs(xs0, ys0, theta, r0, v0, x, y);
        [u_wall, v_wall] = func_dipole_velocity_wall(xs0, ys0, theta, r0, v0, D_tunnel, x, y);
        
        u_dipole = u_fs + real(u_wall); % remove imaginary part, which is very small and should be neglected
        v_dipole = v_fs + real(v_wall);
        
        %% compute error between theory and simulation
        diff_uv = (u_dipole - u_fish_mean).^2 +  (v_dipole - v_fish_mean).^2;
        N_uv = sum(~isnan(diff_uv(:))); % number of non-NAN points in the domain
        error_uv(Ixs, Ir0) = sqrt(  sum( diff_uv(:), 'omitnan' ) / ( N_uv )  ); % compute error; see notes for the definition; ignore nan values
        %}
        
        %% plot
        %{
            figure(30);
            clf;
            hold on;
            pcolor(x, y, sqrt(diff_uv));
            %         pcolor(x, y, tor_max_mag);
            shading flat
            colorbar;
            axis equal tight
            xlabel('$x$ (m)','Interpreter','latex')
            ylabel('$y$ (m)','Interpreter','latex')
            set(gcf,'Color','w');
            set(gca,'fontsize',FTsz);
            set(gcf,'position', [189, 546, 1333, 257]);
            colormap jet
            %             caxis([0, 0.3]);
            
            pause(0.0001)
        %}
        
        %{
        uv_mag = sqrt(u_dipole.^2 + v_dipole.^2); % velocity magnitude

        figure(40);
        clf;
        hold on;
        pcolor(x, y, uv_mag);
%         pcolor(x, y, sqrt(diff_uv));
        shading flat
        h1=streamslice(x, y, u_dipole, v_dipole, 4, 'arrows');
        clr_str = [1,1,1];
        set( h1, 'Color', clr_str )
        set(h1, 'linewidth', 1)
        colorbar;
        plot(xs0, ys0, 'kd', 'linewidth', 0.5);
        
        axis equal tight
        xlabel('$x$ (m)','Interpreter','latex')
        ylabel('$y$ (m)','Interpreter','latex')
        set(gcf,'Color','w');
        set(gca,'fontsize',FTsz);
        set(gcf,'position', [207   211   859   444]);
        colormap jet
        %}
        %{
        % flow across width
        y_tar = linspace(0,D_tunnel,1025);
        x_tar = ones(1,length(y_tar))*(xs0 + L_sheet*0.5);
        
        u_int = interp2(x, y, u_dipole, x_tar, y_tar);
        v_int = interp2(x, y, v_dipole, x_tar, y_tar);
        
        figure(67)
        clf;
        hold on
        plot(y_tar, u_int, 'b.-', 'linewidth', LNwd);
        plot(y_tar, v_int, 'r.-', 'linewidth', LNwd);
        xlabel('$y$ (m)','Interpreter','latex')
        ylabel('$u,v$ (m/s)','Interpreter','latex')
        legend('$u$', '$v$', 'Interpreter','latex', 'location', 'best')
        set(gcf,'Color','w');
        set(gca,'fontsize',20);
        box on
        
        pause(0.001)
        
        %}
        
    end % for Ir0
    
end % for Ix0

toc

%% determine optimal solution
[~, Ixmin] = min(error_uv(:));

[r0s, x0s] = meshgrid(r0_all, xs_all);

xs_opt = x0s(Ixmin);
r0_opt = r0s(Ixmin);

Ir0_opt = find(r0_all==r0_opt); % dimension corresponding to the optimal L solution

fprintf('Optimal slice from r0: # %d / %d. \n', Ir0_opt, Nr0);
%% plot simulation
Ix_seedr = 1:4:size(x,1);
Ix_seedc = 1:4:size(x,2);

u_fish_mag = sqrt(u_fish_mean.^2 + v_fish_mean.^2); % velocity magnitude

figure(10);
clf;
hold on;
pcolor(x_simul-x_head, y_simul, u_fish_mag);
shading flat
colorbar('XTick', [0:0.02:0.06,0.07]);
axis equal tight
h1=streamslice(x_simul-x_head, y_simul, u_fish_mean, v_fish_mean, 2, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
caxis([0, 0.07]); % limits for vel magnitude
% set(cx, 'XTick', [0:0.02:0.06, 0.07]);
colormap jet
set(gcf,'position', [189, 546, 1333, 257]);

%
%}
%% visualize error
figure(71)
clf;
hold on;
pcolor(x0s, r0s, error_uv);
shading flat
colorbar;
plot(xs_opt, r0_opt, 'ko', 'markerfacecolor', 'w', 'markersize', MKsz,'linewidth',LNwd);
axis tight
xlabel('$x_0$ (m)','Interpreter','latex')
ylabel('$r_0$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [468   242   637   444]);
colormap jet

%% plot the optimal solution
% compute velocity field
[u_fs, v_fs] = func_dipole_velocity_fs(xs_opt, ys0, theta, r0_opt, v0, x, y);
[u_wall, v_wall] = func_dipole_velocity_wall(xs_opt, ys0, theta, r0_opt, v0, D_tunnel, x, y);

u_opt = u_fs + real(u_wall);
v_opt = v_fs + real(v_wall);

uv_mag = sqrt(u_opt.^2 + v_opt.^2); % velocity magnitude

% compute the error between theory and simulation
diff_uv = (u_opt - u_fish_mean).^2 +  (v_opt - v_fish_mean).^2;

%% plot velocity field
Ix_mask = isnan(u_fish_mag);
uv_mag(Ix_mask) = nan;
diff_uv(Ix_mask) = nan;

figure(59);
clf;
hold on;
pcolor(x-x_head, y, uv_mag);
shading flat
colorbar('XTick', [0:0.02:0.06,0.07]);
h1=streamslice(x-x_head, y, u_opt, v_opt, 2, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

plot(xs_opt-x_head, ys0, 'kd', 'markerfacecolor', 'k', 'linewidth',0.5)

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [94, 327, 1333, 258]);
colormap('jet');
caxis([0, 0.07]);
axis equal tight
% set(clrmp,'Ticks',0:0.02:0.12);
% set(clrmp,'TickLabels', 0:0.02:0.12);

figure(60);
clf;
hold on;
pcolor(x-x_head, y, sqrt(diff_uv));
shading flat
colorbar('XTick', [0:0.02:0.06,0.07]);
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [94, 327, 1333, 258]);
colormap jet
caxis([0, 0.07]); %

%% flow across width
%
y_tar = linspace(0,D_tunnel,1025);
x_tar = ones(1,length(y_tar))*(x_head+0.007);  % xs_opt

u_int = interp2(x, y, u_opt, x_tar, y_tar);
v_int = interp2(x, y, v_opt, x_tar, y_tar);

figure(67)
clf;
hold on
plot(y_tar, u_int, 'b.-', 'linewidth', LNwd);
plot(y_tar, v_int, 'r.-', 'linewidth', LNwd);
plot((ys0 - r0_opt/2)*[1,1], [min(u_int),max(u_int)], 'k--', 'linewidth', LNwd);
plot((ys0 + r0_opt/2)*[1,1], [min(u_int),max(u_int)], 'k--', 'linewidth', LNwd);
xlabel('$y$ (m)','Interpreter','latex')
ylabel('$u,v$ (m/s)','Interpreter','latex')
legend('$u$', '$v$', 'Interpreter','latex', 'location', 'best')
set(gcf,'Color','w');
set(gca,'fontsize',20);
box on
%}

%% cross-width plot, comparison with simulation results
NL_tar = 5;
N_seed = 2^8;

L_tars = x_head + L_fish*linspace(0,1,NL_tar);
L_tars(1) = L_tars(1) + 1e-3;
L_tars(end) = L_tars(end) - 1e-3;

y_tar = linspace(0,D_tunnel,N_seed);

CM=hsv(NL_tar);
figure(75)
clf;
figure(76)
clf;
for IL = 1:NL_tar
    
    x_tar = ones(1,length(y_tar))*L_tars(IL);
    
    u_int = interp2(x, y, u_dipole, x_tar, y_tar);
    v_int = interp2(x, y, v_dipole, x_tar, y_tar);
    
    u_int_simul = interp2(x_simul, y_simul, u_fish_mean, x_tar, y_tar);
    v_int_simul = interp2(x_simul, y_simul, v_fish_mean, x_tar, y_tar);
    
    figure(75)
    hold on
    plot3(x_tar-x_head, y_tar, u_int, '-', 'color', CM(IL, :), 'linewidth', LNwd);
    plot3(x_tar-x_head, y_tar, u_int_simul, '--', 'color', CM(IL, :), 'linewidth', LNwd);
    
    figure(76)
    hold on
    plot3(x_tar-x_head, y_tar, v_int, '-', 'color', CM(IL, :), 'linewidth', LNwd);
    plot3(x_tar-x_head, y_tar, v_int_simul, '--', 'color', CM(IL, :), 'linewidth', LNwd);
    
end

figure(75)
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('$u$ (m/s)','Interpreter','latex')
xlim([0, max(L_tars)-x_head]);
ylim([0, D_tunnel]);
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
view([1,-1,-1])
set(gcf,'position', [422   226   834   572]);
% box on

figure(76)
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
zlabel('$v$ (m/s)','Interpreter','latex')
xlim([0, max(L_tars)-x_head]);
ylim([0, D_tunnel]);
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
view([1,-1,-1])
set(gcf,'position', [422   226   834   572]);
% box on
%% save data
if flag_save > 0
    fname = [save_dir, 'fit_point_dipole_para_opt.mat'];
    
    save(fname, 'xs_opt', 'r0_opt', 'x0s', 'r0s', 'error_uv', 'Ir0_opt', 'u_opt', 'v_opt');
    
    fprintf('Results saved. \n');
end
