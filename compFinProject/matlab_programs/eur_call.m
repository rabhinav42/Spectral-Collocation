clear all;
clc;

%%% model parameters %%%
sig = 0.2;
r = 0.05;
T = 0.5;
E = 10;
m = 5; % no of intervals
N = 8; % no of CGL points - 1
tol = 1e-5;
h0 = 1e-2;
opt = 'c'; % set option to european call
C = 0; % no penalty term
ep = 0;

%%% set up auxilliary variables for solving IVP %%%

[D1, D2] = lagrangeD(N);
Xi = xi(N); % the CGL points
smin = 0;
smax = m*E/floor(m/2);
delx = (smax - smin)/m;
s = [smin:delx:smax];
S = [];
for i=1:m
    s0 = s(i);
    s1 = s(i+1);
    S = [S X(Xi(1:N), s0, s1)];
end
S = [S smax];

%%% solve ivp %%%

[u,t, cpu_t] = solve_ivp(tol, h0, opt, sig, r, E, T, m, N, D1, D2, s, ep, C);

tt = T - t;
cact = call_act(S, tt, E, r, sig, 0, T); % actual solution

% %%% solve ivp with constant step size %%%  uncomment this part for
% %%% obtaining constant step size solution for given h0 as well.
% h0 = 5*1e-5;
% [u, t, cpu_t] = solve_ivp_const(h0, opt, sig, r, E, T, m ,N, D1, D2, s, ep, C); 
% tt = T - t;
% cact = call_act(S, tt, E, r, sig, 0 , T);

%%% print needed values %%%
err = max(abs(u(end,:) - cact(end,:)))
accepted_steps = length(tt)-1
cpu_t

%%% error plots etc %%%
figure(1)
plot(S, abs(u(end,:) - cact(end,:)));
xlabel('S')
ylabel('Error(S)')
saveas(figure(1), "err_eur.jpg");

figure(2)
plot(S, cact(end,:));
hold on;
scatter(S, u(end,:), 'black', 'o');
hold off;
xlabel('S')
ylabel('V(S,T)')
legend('Actual', 'Numerical', 'Location', 'northwest')
saveas(figure(2), "eur_vs.jpg");