% this code is modified from main_fit_sheet_simul_v9.m
% flow field is fitted on the side of the body
% sheet pair is placed along the fish body

clear
close all
clc

flag_save = -1;

FTsz = 20;
LNwd = 1;
MKsz = 10;

%% directories
save_dir = '../analysis_results_dipole/';
simul_results_dir = '../../fish_swim/analysis_results/';

%% vortex sheet parameters
L_fish = 0.035481891661866; % fish body length
D_tunnel = 0.0457; % m; tunnel diameter
L_tunnel = 0.15; % m; tunnel length

L_domain = 2*L_tunnel; % length of fluid domain

theta_s = 1*pi; % orientation of sheet pair; angle with x-axis

TBA = (0.002997180904265 + 0.020600558880515 - 0.040537635755492 + 0.097600742072607 )*L_fish; % m; corresponding to tail beat amplitude of 4 mm
v0 = 0.156379646024586; % m/s; fish fish swimming speed at lowest flow rate
Gamma = 2*pi*TBA*v0; % m^2/s; treated circulation of one vortex in a dipole; defined in Porfiri PRL 2021
r0 = 2*TBA;

%% fluid domain seed points
dx = 1e-3;
dy = dx;

Nx = round(L_domain/dx);
Ny = round(D_tunnel/dy);

X = linspace(0, L_domain, Nx);
Y = linspace(0, D_tunnel, Ny);

[x, y] = meshgrid(X, Y);

dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);

%% iterate through parameter space
L_sheet = L_domain/2; % length of one sheet

Nelem = 16; % determine number of elements; set a lower threshold

xs0 = L_domain/4;
ys0 = D_tunnel/4; %

gamma = Gamma/L_sheet; % (constant) curculation density (per unit length) of one vortex sheet

%%  core size for regularization
% dist_sheet = 2*TBA;
% dist_sheet = 5e-10; % separation distance between 2 sheets
delta_core = dx/4;

%% compute velocity field
% end point of vortex sheet pair
xs_end = xs0 + L_sheet;
ys_end = ys0;

delem = ones(1, Nelem)*L_sheet/Nelem; % size of each element

% element edge point locations
xs = linspace(xs0, xs_end, Nelem+1);
ys = linspace(ys0, ys_end, Nelem+1);

% element mid-point locations
xe = 0.5*(xs(2:end) + xs(1:end-1));
ye = 0.5*(ys(2:end) + ys(1:end-1));
ye1 = ye; % + dist_sheet/2;
% ye2 = ye - dist_sheet/2;

%% sheet element quantaties
gamma_elem = gamma.*ones(1, Nelem); % curculation density (per unit length) of each element

%% compute velocity field
% allocate velocity array at seed locations
u_all = zeros(size(x));
v_all = zeros(size(x));

% iterate over all elements
for Ielem = 1:Nelem
    
    %% circulation of element Ielem
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    %% compute velocity field
    
    [u_f1, v_f1] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  x, y, delta_core); % _desingul
    
    %     u = u_f1 + u_f2;
    %     v = v_f1 + v_f2;
    
    %% combine all dipoles velocity
    u_all = u_all + u_f1;
    v_all = v_all + v_f1;
    
end

uv_mag = sqrt(u_all.^2 + v_all.^2 );
%% sheet self velocity
u_self = zeros(Nelem, 1);
v_self = zeros(Nelem, 1);

% iterate over all elements
for Ielem = 1:Nelem % loop over all target elements
    
    %% circulation of element Ielem
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    [u_sheet, v_sheet] = func_vortex_velocity_walls_regularize(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  xe, ye1, delta_core); % _desingul
    
    u_self = u_self + u_sheet(:);
    v_self = v_self + v_sheet(:);
    
end % for Ielem

%% plot
%
figure(40);
clf;
hold on;
pcolor(x, y, uv_mag);
%         pcolor(x, y, sqrt(diff_uv));
shading flat
h1=streamslice(x, y, u_all, v_all, 4, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)
colorbar;
% plot(xe, ye+r0/2, 'k.','markersize', MKsz+3, 'linewidth', LNwd+3);
% plot(xe, ye-r0/2, 'kx','markersize', MKsz,'linewidth', LNwd);
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [207   211   859   444]);
colormap jet
caxis([0, 6e-2])


figure(50)
clf;
hold on;
plot(xe, u_self, 'k-', 'linewidth', LNwd+3);
plot(xe, v_self, 'r-', 'linewidth', LNwd+3);
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$u,v$ (m/s)','Interpreter','latex')
legend('$u$ (regularized)', '$v$ (regularized)', 'Interpreter','latex', 'location', 'best')
set(gcf,'Color','w');
set(gca,'fontsize',20);
box on

%% flow across width
%{
y_tar = linspace(0,D_tunnel,1025);
x_tar = ones(1,length(y_tar))*(xs_opt + L_opt*0.5);

u_int = interp2(x, y, u_sheet, x_tar, y_tar);
v_int = interp2(x, y, v_sheet, x_tar, y_tar);

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute velocity without regularization, but with removed singularity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute velocity field
% allocate velocity array at seed locations
u_all_singl = zeros(size(x));
v_all_singl = zeros(size(x));

for Ielem = 1:Nelem
    
    %% circulation of element Ielem
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    %% compute velocity field
    
    [u_f1, v_f1] = func_vortex_velocity_walls(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  x, y); % _desingul
    
    %     u = u_f1 + u_f2;
    %     v = v_f1 + v_f2;
    
    %% combine all dipoles velocity
    u_all_singl = u_all_singl + u_f1;
    v_all_singl = v_all_singl + v_f1;
    
end

uv_mag_singl = sqrt(u_all_singl.^2 + v_all_singl.^2 );

%% sheet self velocity
u_self_singl = zeros(Nelem, 1);
v_self_singl = zeros(Nelem, 1);

% iterate over all elements
for Ielem = 1:Nelem % loop over all target elements
    
    %% circulation of element Ielem
    Gamma_elem = delem(Ielem)*gamma_elem(Ielem); % circulation of one element segment on the vortex sheet pair
    
    % compute velocity induced by other elements Jelem on the same sheet
    for Jelem = 1:Nelem % 
        Gamma_elem_J = delem(Jelem)*gamma_elem(Jelem);  % strength of the other element on the sheet
        
        if Jelem ~= Ielem
            [u_sheet, v_sheet] = func_vortex_velocity_walls(xe(Jelem), ye1(Jelem), D_tunnel, Gamma_elem_J,  xe(Ielem), ye1(Ielem)); % _desingul
        else
            k = pi/D_tunnel;
            yim = -ye1(Ielem);
            %             [u_sheet, v_sheet] = func_vortex_velocity_walls_desingul_self(xe(Ielem), ye1(Ielem), D_tunnel, Gamma_elem,  x, y, delta_core); % _desingul
            u_sheet = (0.25*Gamma_elem/pi)*k*sin(k*(ye1(Ielem)-yim))./ (1 - cos(k*(ye1(Ielem)-yim)));
            v_sheet = 0;
        end
        
        u_self_singl(Ielem) = u_self_singl(Ielem) + u_sheet;
        v_self_singl(Ielem) = v_self_singl(Ielem) + v_sheet;
    end
    
end % for Ielem

%% plot
%
figure(41);
clf;
hold on;
pcolor(x, y, uv_mag_singl);
%         pcolor(x, y, sqrt(diff_uv));
shading flat
h1=streamslice(x, y, u_all, v_all, 4, 'arrows');
clr_str = [1,1,1];
set( h1, 'Color', clr_str )
set(h1, 'linewidth', 1)
colorbar;
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [207   211   859   444]);
colormap jet
caxis([0, 6e-2])

figure(50)
% clf;
hold on;
plot(xe, u_self_singl, 'm:', 'linewidth', LNwd+3);
plot(xe, v_self_singl, 'c:', 'linewidth', LNwd+3);
% xlabel('$x$ (m)','Interpreter','latex')
% ylabel('$u,v$ (m/s)','Interpreter','latex')
legend('$u$ (regularized)', '$v$ (regularized)', '$u$ (de-singularized)', '$v$ (de-singularized)', 'Interpreter','latex', 'location', 'best')
% set(gcf,'Color','w');
% set(gca,'fontsize',20);
box on
set(gcf,'position', [396   251   752   505]);
%% plot difference between de-singular and regularized fields
diff_mag = sqrt((u_all_singl-u_all).^2 + (v_all_singl-v_all).^2 );

figure(42);
clf;
hold on;
pcolor(x, y, diff_mag);
%         pcolor(x, y, sqrt(diff_uv));
shading flat
colorbar;
axis equal tight
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$y$ (m)','Interpreter','latex')
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);
set(gcf,'position', [207   211   859   444]);
colormap jet
caxis([0,1e-4])

figure(52)
clf;
hold on;
plot(xe, abs(u_self_singl-u_self), 'k-', 'linewidth', LNwd+3);
plot(xe, abs(v_self_singl-v_self), 'r-', 'linewidth', LNwd+3);
xlabel('$x$ (m)','Interpreter','latex')
ylabel('$|\Delta u|, |\Delta v|$ (m/s)','Interpreter','latex')
legend('$|\Delta u|$', '$|\Delta v|$', 'Interpreter','latex', 'location', 'best')
set(gcf,'Color','w');
set(gca,'fontsize',20);
box on
set(gcf,'position', [396   251   752   505]);
%% save data

if flag_save > 0
    
%     fname = [save_dir, 'fit_sheet_para_opt.mat'];
%     
%     save(fname, 'xs_opt', 'gamma_opt', 'L_opt');
    
end
