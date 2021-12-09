% this code is modified from main_fit_sheet_simul_v11.m
% modified to fit velocity fields at other swimming speeds in the
% parametric study

clear
close all
clc

flag_save = -1;

FTsz = 20;
LNwd = 2;
MKsz = 10;

%% directories
save_dir = '../analysis_results_dipole/';
simul_results_dir = '../../fish_swim/analysis_results/';

%% read simulation results
U = 0.05; % m/s; swimming speed of fish used in simulation

fname = [simul_results_dir, 'simul_vel_opt_paraU_', num2str(U),'.mat'];

load(fname, 'x', 'y', 'u_fish_mean', 'v_fish_mean');

x_simul = x - min(x(1,:));
y_simul = y - min(y(:,1));

%% vortex sheet parameters
L_fish = 0.035481891661866; % fish body length
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain

theta_s = 1*pi; % orientation of sheet pair; angle with x-axis

TBA = (0.002997180904265 + 0.020600558880515 - 0.040537635755492 + 0.097600742072607 )*L_fish; % m; corresponding to tail beat amplitude of 4 mm
U_swim = U; % m/s; fish fish swimming speed at lowest flow rate
% Gamma = 2*pi*TBA*U_swim; % m^2/s; estimated circulation of one vortex in a dipole
% r0 = 2*TBA; %

% location of the starting point of vortex sheet pair
x_head = L_tunnel/2; % location of the fish head;
ys0 = D_tunnel/2; %
% xs0 = 0; % 1st dipole position

%% define parameter space
gamma0 = 1e-2/(L_fish); % initial guess of dimensional circulation density; normalized by body length

% length of sheet pair
NL = 10+1;
L_max = L_fish; % 1.5*L_fish; % max length of sheet pair; make sure that the entire sheet pair is in the fluid domain
% L_sheet_all = linspace(0.00, L_max, NL);
L_sheet_all = linspace(0.028, L_max, NL);
% L_sheet_all(1) = 0.5*(L_sheet_all(1)+L_sheet_all(2)); % remove trivial case

% starting point of the sheet pair
Nxs = 10+1;
% xs_all = linspace(x_head-0.2*L_fish, x_head+1.2*L_fish, Nxs); % coarse
xs_all = linspace(0.073, 0.078, Nxs); % fine

% separation distance
Ndist = 10+1;
% dist_all = linspace(0, TBA, Ndist); % coarse
dist_all = linspace(0.0015, 0.0022, Ndist); % fine
% dist_all(1) = 0.5*(dist_all(1)+dist_all(2));  % remove trivial case

%% fluid domain seed points
x = x_simul; % use same seeding points for simulation
y = y_simul;

dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);

delta_core = min([min(dist_all), dx])/4; %core size; making sure vortex cores on the same sheet and the other sheet would not overlap each other

%% iterate through parameter space
error_uv = nan(NL, Nxs, Ndist); % assign error function
Gamma_true = nan(NL, Nxs, Ndist); % assign error function

tic
for IL = 1:NL
    
    L_sheet = L_sheet_all(IL); % length of one sheet
    
    Nelem = max([round(L_sheet/(1e-3)), 32]); % determine number of elements; set a lower threshold
    
    fprintf('Computing %d of %d values of L... Number of elements: %d. \n', IL, NL, Nelem);
    
    for Ixs = 1:Nxs
        xs0 = xs_all(Ixs);
        
        for Idist = 1:Ndist
            
            %% vortex sheet parameters
            dist_sheet = dist_all(Idist); % (constant) curculation density (per unit length) of one vortex sheet
            
            %% compute velocity field
            %% sheet location and discretization
            % end point of vortex sheet pair
            xs_end = xs0 + L_sheet;
            ys_end = ys0 ;
            
            delem = ones(1, Nelem)*L_sheet/Nelem; % size of each element
            
            % element edge point locations
            xs = linspace(xs0, xs_end, Nelem+1);
            ys = linspace(ys0, ys_end, Nelem+1);
            
            % element mid-point locations
            xe = 0.5*(xs(2:end) + xs(1:end-1));
            ye = 0.5*(ys(2:end) + ys(1:end-1));
            ye1 = ye + dist_sheet/2;
            ye2 = ye - dist_sheet/2;
            %% sheet element quantaties
            gamma_elem = gamma0.*ones(1, Nelem); % curculation density (per unit length) of each element
            
            %% sheet self velocity
            %{
                u_self1 = zeros(Nelem, 1);
                v_self1 = zeros(Nelem, 1);
                u_self2 = zeros(Nelem, 1);
                v_self2 = zeros(Nelem, 1);
                
                % iterate over all elements on sheet 1
                for Ielem = 1:Nelem % loop over all target elements

                    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
                    
                    % effect of element i on sheet 1, on the velocity of all points on sheet 1 and 2
                    [u_sheet1, v_sheet1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, -Gamma_elem,  xe, ye1, delta_core); % _desingul
                    
                    [u_sheet2, v_sheet2] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  xe, ye2, delta_core); % _desingul
                    
                    u_self1 = u_self1 + u_sheet1(:);
                    v_self1 = v_self1 + v_sheet1(:);
                    
                    u_self2 = u_self2 + u_sheet2(:);
                    v_self2 = v_self2 + v_sheet2(:);
                    
                end % for Ielem
                
                % iterate over all elements on sheet 2
                for Ielem = 1:Nelem % loop over all target elements

                    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
                    
                    % effect of element i on sheet 2, on the velocity of all points on sheet 1 and 2
                    [u_sheet1, v_sheet1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye2(Ielem), D_tunnel, -Gamma_elem,  xe, ye1, delta_core); % _desingul
                    
                    [u_sheet2, v_sheet2] = func_vortex_velocity_walls_regularize(xe(Ielem), ye2(Ielem), D_tunnel, Gamma_elem,  xe, ye2, delta_core); % _desingul
                    
                    u_self1 = u_self1 + u_sheet1(:);
                    v_self1 = v_self1 + v_sheet1(:);
                    
                    u_self2 = u_self2 + u_sheet2(:);
                    v_self2 = v_self2 + v_sheet2(:);
                    
                end % for Ielem
            %}
            [u_self1, v_self1, u_self2, v_self2] = func_self_vel_sheet_pair_walls_regul...
                (xe, ye1, ye2, delem, gamma_elem, D_tunnel, delta_core);
            
            u_self_mean = mean(u_self1);
            v_self_mean = mean(v_self1);
            
            %{
                figure(43)
                clf;
                hold on;
                plot(xe, u_self1, 'k-', 'linewidth', LNwd+3);
                plot(xe, v_self1, 'r-', 'linewidth', LNwd+3);
                plot(xe, u_self2, 'm:', 'linewidth', LNwd+3);
                plot(xe, v_self2, 'c:', 'linewidth', LNwd+3);
                xlabel('$x$ (m)','Interpreter','latex')
                ylabel('$u,v$ (m/s)','Interpreter','latex')
                legend('$u_1$', '$v_1$', '$u_2$', '$v_2$', 'Interpreter','latex', 'location', 'best')
                set(gcf,'Color','w');
                set(gca,'fontsize',20);
                box on
                set(gcf,'position', [396   251   752   505]);
            %}
            
            % compute the true value of gamma, based on the linear
            % relationship between velocity and circulation density
            gamma_true = gamma0*(-U_swim)/u_self_mean;
            
            gamma_elem = gamma_true.*ones(1, Nelem); % re-define elemental circulation density
            
            Gamma_true(IL, Ixs, Idist) = gamma_true; % store circulation
            %% validate the density by computing new self-velocity; needed only during debugging process
            %{

                [u_self1, v_self1, u_self2, v_self2] = func_self_vel_sheet_pair_walls_regul...
                    (xe, ye1, ye2, delem, gamma_elem, D_tunnel, delta_core);

                u_self_mean = mean(u_self1);
                v_self_mean = mean(v_self1);
            %}
            %% compute velocity field
            
            [u_sheet, v_sheet] = func_vel_sheet_pair_walls_regul ...
                (xe, ye1, ye2, D_tunnel, delem, gamma_elem, x, y, delta_core);
            %  uv_mag = sqrt(u_sheet.^2 + v_sheet.^2); % velocity magnitude
            
            %% compute error between theory and simulation
            diff_uv = (u_sheet - u_fish_mean).^2 +  (v_sheet - v_fish_mean).^2;
            N_uv = sum(~isnan(diff_uv(:))); % number of non-NAN points in the domain
            error_uv(IL, Ixs, Idist) = sqrt(  sum( diff_uv(:), 'omitnan' ) / ( N_uv )  ); % compute error; see notes for the definition; ignore nan values
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
        figure(40);
        clf;
        hold on;
        pcolor(x, y, uv_mag);
%         pcolor(x, y, sqrt(diff_uv));
        shading flat
        h1=streamslice(x, y, u_sheet, v_sheet, 4, 'arrows');
        clr_str = [1,1,1];
        set( h1, 'Color', clr_str )
        set(h1, 'linewidth', 1)
        colorbar;
        plot(xe, ye+r0/2, 'k.','markersize', MKsz+3, 'linewidth', LNwd+3);
        plot(xe, ye-r0/2, 'kx','markersize', MKsz,'linewidth', LNwd);
        
        axis equal tight
        xlabel('$x$ (m)','Interpreter','latex')
        ylabel('$y$ (m)','Interpreter','latex')
        set(gcf,'Color','w');
        set(gca,'fontsize',FTsz);
        set(gcf,'position', [207   211   859   444]);
        colormap jet
        
        % flow across width
        y_tar = linspace(0,D_tunnel,1025);
        x_tar = ones(1,length(y_tar))*(xs0 + L_sheet*0.5);
        
        u_int = interp2(x, y, u_sheet, x_tar, y_tar);
        v_int = interp2(x, y, v_sheet, x_tar, y_tar);
        
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
            
        end % for Idist
        
    end % for Ix0
    
end % for IL

toc

%% determine optimal solution
[~, Ixmin] = min(error_uv(:));

[x0s, Ls, Ds] = meshgrid(xs_all, L_sheet_all, dist_all);

xs_opt = x0s(Ixmin);
dist_sheet_opt = Ds(Ixmin);
L_opt = Ls(Ixmin);

gamma_opt = Gamma_true(Ixmin); % true circulation corresponding to opt distance

Idist_opt = find(dist_all==dist_sheet_opt); % dimension corresponding to the optimal L solution

fprintf('Optimal slice from distance: # %d / %d. \n', Idist_opt, Ndist);
%% plot simulation
Ix_seedr = 1:4:size(x,1);
Ix_seedc = 1:4:size(x,2);

u_fish_mag = sqrt(u_fish_mean.^2 + v_fish_mean.^2); % velocity magnitude

figure(10);
clf;
hold on;
pcolor(x_simul-x_head, y_simul, u_fish_mag);
shading flat
colorbar; %('XTick', [0:0.02:0.06,0.07]);
axis equal tight
% quiver(x(Ix_seedr, Ix_seedc), y(Ix_seedr, Ix_seedc), u_mean(Ix_seedr, Ix_seedc), v_mean(Ix_seedr, Ix_seedc),'w');
%     set(gcf,'position', [329, 129, 1213, 847]);
h1=streamslice(x_simul-x_head, y_simul, u_fish_mean, v_fish_mean, 2, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
% caxis([0, 0.07]); % limits for vel magnitude
% set(cx, 'XTick', [0:0.02:0.06, 0.07]);
colormap jet
set(gcf,'position', [94, 327, 1333, 258]);

%
%}

%{
figure(13);
clf;
hold on;
pcolor(x, y, tor_max_angle_simul);
shading flat
colorbar;
% h1=streamslice(x, y, u_sheet, v_sheet, 1, 'arrows');
% clr_str = [1,1,1];
% set( h1, 'Color', clr_str )
% set(h1, 'linewidth', 1)
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [189         546        1333         257]);
colormap jet
caxis([-pi/2, pi/2]);

figure(14);
clf;
hold on;
pcolor(x, y, tor_min_angle_simul);
shading flat
colorbar;
% h1=streamslice(x, y, u_sheet, v_sheet, 1, 'arrows');
% clr_str = [1,1,1];
% set( h1, 'Color', clr_str )
% set(h1, 'linewidth', 1)
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [189         546        1333         257]);
colormap jet
caxis([-pi/2, pi/2]);
%}

%% visualize error
figure(71)
clf;
hold on;
pcolor(Ls(:,:,Idist_opt), x0s(:,:,Idist_opt), error_uv(:,:,Idist_opt));
shading flat
colorbar;
plot(L_opt, xs_opt, 'ko', 'markerfacecolor', 'w', 'markersize', MKsz, 'linewidth',LNwd);
axis tight
xlabel('$L$ (m)','Interpreter','latex')
ylabel('$x0$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [468   242   637   444]);
colormap jet

%     pause(0.1)
% end
%% plot the optimal solution
Nelem = max([round(L_opt/(1e-3)), 32]); % determine number of elements; set a lower threshold

% end point of vortex sheet pair
xs_end = xs_opt + L_opt;
ys_end = ys0;

delem = ones(1, Nelem)*L_opt/Nelem; % size of each element

% element end point locations
xs = linspace(xs_opt, xs_end, Nelem+1);
ys = linspace(ys0, ys_end, Nelem+1);
ys1 = ys + dist_sheet_opt/2;
ys2 = ys - dist_sheet_opt/2;

% element mid-point locations
xe = 0.5*(xs(2:end) + xs(1:end-1));
ye = 0.5*(ys(2:end) + ys(1:end-1));
ye1 = ye + dist_sheet_opt/2;
ye2 = ye - dist_sheet_opt/2;

% sheet element quantaties
gamma_elem = gamma_opt.*ones(1, Nelem);

[u_opt, v_opt] = func_vel_sheet_pair_walls_regul ...
    (xe, ye1, ye2, D_tunnel, delem, gamma_elem, x, y, delta_core);

uv_mag = sqrt(u_opt.^2 + v_opt.^2); % velocity magnitude

% compute the error between theory and simulation
diff_uv = (u_opt - u_fish_mean).^2 +  (v_opt - v_fish_mean).^2;
% diff_uv = (tor_max_mag - tor_max_mag_simul).^2 +  (tor_min_mag - tor_min_mag_simul).^2;

%% validate that optimal solution gives matching swimming speed
[u_self1, v_self1, u_self2, v_self2] = func_self_vel_sheet_pair_walls_regul...
    (xe, ye1, ye2, delem, gamma_elem, D_tunnel, delta_core);

fprintf('Check: swimming speed at optimal solution: %18.15f m/s. \n', mean(u_self1) );

%% plot optimal fitted velocity field
Ix_mask = isnan(u_fish_mag);
uv_mag(Ix_mask) = nan;
diff_uv(Ix_mask) = nan;

%{
figure(58);
clf;
hold on;
pcolor(x, y, tor_max_mag);
shading flat
colorbar('XTick', -15:5:15);
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [189, 546, 1333, 257]);
colormap jet
caxis([-15, 15]);

figure(59);
clf;
hold on;
pcolor(x, y, tor_min_mag);
shading flat
colorbar('XTick', -15:5:15);
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [189, 546, 1333, 257]);
colormap jet
caxis([-15, 15]);
%}

figure(59);
clf;
hold on;
pcolor(x-x_head, y, uv_mag);
shading flat
colorbar; %('XTick', [0:0.02:0.06,0.07]);
h1=streamslice(x-x_head, y, u_opt, v_opt, 2, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)

plot(xs-x_head, ys1, 'k:', 'linewidth',LNwd)
plot(xs-x_head, ys2, 'k-', 'linewidth',LNwd)

xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [94, 327, 1333, 258]);
colormap('jet');
% caxis([0, 0.07]);
axis equal tight
% set(clrmp,'Ticks',0:0.02:0.12);
% set(clrmp,'TickLabels', 0:0.02:0.12);

figure(60);
clf;
hold on;
pcolor(x-x_head, y, sqrt(diff_uv));
shading flat
colorbar; %('XTick', [0:0.02:0.06,0.07]);
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [94, 327, 1333, 258]);
colormap jet
% caxis([0, 0.07]); %

%% flow across width
%
y_tar = linspace(0,D_tunnel,1025);
x_tar = ones(1,length(y_tar))*(xs_opt + L_opt*0.5);

u_int = interp2(x, y, u_opt, x_tar, y_tar);
v_int = interp2(x, y, v_opt, x_tar, y_tar);

figure(67)
clf;
hold on
plot(y_tar, u_int, 'b.-', 'linewidth', LNwd);
plot(y_tar, v_int, 'r.-', 'linewidth', LNwd);
plot((ys0-dist_sheet/2)*[1,1], [min(u_int),max(u_int)], 'k--', 'linewidth', LNwd);
plot((ys0+dist_sheet/2)*[1,1], [min(u_int),max(u_int)], 'k--', 'linewidth', LNwd);
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
    
    u_int = interp2(x, y, u_opt, x_tar, y_tar);
    v_int = interp2(x, y, v_opt, x_tar, y_tar);
    
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
% legend('$u$', '$v$', 'Interpreter','latex', 'location', 'best')
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
% legend('$u$', '$v$', 'Interpreter','latex', 'location', 'best')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
view([1,-1,-1])
set(gcf,'position', [422   226   834   572]);
% box on
%% save data
if flag_save > 0
    fname = [save_dir, 'fit_sheet_para_paraU_', num2str(U),'.mat']; 
    
    save(fname, 'L_opt', 'xs_opt', 'dist_sheet_opt', 'gamma_opt', 'delta_core', 'x0s', 'Ls', 'Ds', ...
        'error_uv', 'Idist_opt', 'u_opt', 'v_opt', 'u_self1', 'v_self1', 'u_self2', 'v_self2', ...
        'xe', 'ye1', 'ye2', 'xs', 'ys1', 'ys2');
    
    fprintf('Results saved. \n');
end
