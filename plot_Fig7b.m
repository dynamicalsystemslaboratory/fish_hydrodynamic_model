% plot Fig 7 b

clear
close all
clc

load('Data_Fig7b.mat', 'gamma'); 
%% parameters
L_fish = 0.035481891661866; % fish body length

U_all = [0.05, 0.1, 0.156379646024586, 0.25, 0.5]; % m/s; all speeds

NU = length(U_all); 

%% plot

figure(11)
hold on
plot(U_all, gamma, 'ko-','markerfacecolor', 'k');
xlabel('$U\,(\mathrm{m}/\mathrm{s})$','Interpreter','latex')
ylabel('$\Gamma\,(\mathrm{m^2}/\mathrm{s})$','Interpreter','latex')
