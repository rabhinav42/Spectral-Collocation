clear all;
clc;

%%% model parameters %%%
sig = 0.3;
r = 0.1;
T = 1;
E = 1;
m = 7; % no of intervals
N = 8; % no of CGL points - 1
tol = 1e-5;
h0 = 1e-2;
opt = 'p'; % set option to american put
C = r*E; % penalty term parameters
ep = 1e-4;
steps = 1000; % no of steps in binomial tree solution


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
pact = put_act(S, T, E, r, sig, steps);

% %%% solve ivp with constant step size %%%  uncomment this part for
% %%% obtaining constant step size solution for given h0 as well.
% h0 = 1e-4;
% [u, t, cpu_t] = solve_ivp_const(h0, opt, sig, r, E, T, m ,N, D1, D2, s, ep, C); 
% tt = T - t;
% pact = put_act(S, T, E, r, sig, steps);

%%% print needed values %%%
err = max(abs(u(end,:) - pact))
accepted_steps = length(tt)-1
cpu_t

%%% error plots etc %%%
figure(1)
plot(S, abs(u(end,:) - pact));
xlabel('S')
ylabel('Error(S)')
saveas(figure(1), "err_a.jpg");

figure(2)
plot(S, pact);
hold on;
scatter(S, u(end,:), 'black', 'o');
hold off;
xlabel('S')
ylabel('V(S,T)')
legend('Actual', 'Numerical')
saveas(figure(2), "a_vs.jpg");

figure(3)
plot(S, max(E-S,0));
hold on;
plot(S, u(end,:), '.', 'Color', 'black');
hold off;
xlabel('S')
ylabel('V(S,T)')
legend('Payoff', 'Computed V(S,T)')
saveas(figure(3), "vst_a.jpg")