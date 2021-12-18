% plot Fig 7 a

clear
close all
clc

load('Data_Fig7a.mat', 'L_all', 'L_fish', 'dist_all', 'thick_body'); 

%%
U_all = [0.05, 0.1, 0.156379646024586, 0.25, 0.5]'; % m/s; all speeds
NU = length(U_all); 

%% plot

fig = figure(12);
hold on

yyaxis left
plot(U_all, L_all/L_fish, 'bo-','markerfacecolor', 'b');
xlabel('$U\,(\mathrm{m}/\mathrm{s})$','Interpreter','latex')
ylabel('$\mathcal{L}/L$','Interpreter','latex')

yyaxis right
plot(U_all, dist_all/thick_body, 'o--','color', 'r');
ylabel('$r_0/b$','Interpreter','latex')

legend('$\mathcal{L}/L$', '$r_0/b$' ,'Interpreter','latex', 'location', 'northeast')
