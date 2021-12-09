% modified from v1
% mode shape is changed from 2nd order to 3rd order

clear
close all
clc

flag_save = -1; % >0: save results; <0: do not save results
flag_save_video = -1; % >0: save video; <0: do not save video

FTsz = 15;
LNwd = 2;
MKsz = 5;
%% specify case number
fish_all = [5, 6, 11, 13, 14];

Nfish = length(fish_all);
%% directories
results_dir = '..\analysis_results_dipole\';
track_dir = '..\analysis_results_dipole\body_track\';
save_dir = results_dir;

%% load fitting data
fname = [save_dir, 'fish_shape_fit_anguilli.mat'];

load(fname, 'c1', 'c2', 'c3', 'c4', 'kL', 'omega', 'phase_shift', 'y_shift', 'body_length', 'thick2cord');

%% take the mean within fish individuals (for cases where 2 episodes of swimming were tracked)

c1 = mean(c1, 2, 'omitnan');
c2 = mean(c2, 2, 'omitnan');
c3 = mean(c3, 2, 'omitnan');
c4 = mean(c4, 2, 'omitnan');
kL = mean(kL, 2, 'omitnan');
omega = mean(omega, 2, 'omitnan');
% shift = mean(shift, 2, 'omitnan');
phi_0 = mean(phase_shift, 2, 'omitnan');

thick2cord = mean(thick2cord, 2, 'omitnan');
body_length = mean(body_length, 2, 'omitnan');
%%
c1 = mean(c1);
c2 = mean(c2);
c3 = mean(c3);
c4 = mean(c4);
kL = mean(kL);
omega = mean(omega);
shift = 0; % mean(shift);
phi_0 = 0; %mean(phi_0);

thick2cord = mean(thick2cord); % thickness-to-cord ratio
body_length = mean(body_length);

%% build geometry
% half profile of NACA00xx airfoil
Dy =@(x) 5*thick2cord*(0.2969.*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843 * x.^3 - 0.1015*x.^4) ;

% undulation of fish body
Nx = 128;
Nt = 51;

x = linspace(0,1,Nx);
x = x.^2;
times = linspace(0,1,Nt);

%% plot shape over time
if flag_save_video>0
    vid_fname = [save_dir, 'fitted_anguilli_locom.avi']; %fullfile(workingDir,'shuttle_out.avi')
    outputVideo = VideoWriter(vid_fname);
    outputVideo.FrameRate = 5;
    open(outputVideo);
end

% y_init = theta.* sin(k.*x );
y_tail = nan(Nt,1);

delta_y = Dy(x);
delta_y(end) = 0;

for Itime = 1:Nt
    
    y_0 = func_anguilliform(x, times(Itime), c1, c2, c3, c4, kL, omega, shift, phi_0);
    
    y_tail(Itime) = y_0(end);
    
    y_plus = y_0 + delta_y;
    y_minus = y_0 - delta_y;
    
    figure(20);
    clf;
    hold on;
    plot(x, y_plus, 'b.-');
    plot(x, y_minus, 'r.-');
    axis equal
    ylim([-0.5, 0.5])
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    set(gcf,'Color','w');
    set(gca,'fontsize',FTsz);
    
    %% save video
    
    if flag_save_video>0
        vidfrm = getframe(gcf) ;
        writeVideo(outputVideo,vidfrm);
    end
    
    pause(0.2);
end

if flag_save_video > 0
    close(outputVideo);
end

%% save airfoil shape as initial geometry in comsol simulation
y_0 = zeros(size(x)); % assuming no deformation at t=0

delta_y = Dy(x);
delta_y(end) = 0;

y_plus = y_0 + delta_y;
y_minus = y_0 - delta_y;

figure(21);
clf;
hold on;
plot(x, y_plus, 'b.-');
plot(x, y_minus, 'r.-');
axis equal
ylim([-0.5, 0.5])
xlabel('$x/L$','Interpreter','latex');
ylabel('$y/L$','Interpreter','latex');
set(gcf,'Color','w');
set(gca,'fontsize',FTsz);

if flag_save>0
    savefname = [save_dir, 'y_top_0.txt'];
    dlmwrite(savefname,[x(:), y_plus(:)],'delimiter','\t','precision',5);
    
    savefname = [save_dir, 'y_btm_0.txt'];
    dlmwrite(savefname,[x(:), y_minus(:)],'delimiter','\t','precision',5);
end
