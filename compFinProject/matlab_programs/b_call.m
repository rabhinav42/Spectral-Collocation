clear all;
clc;

%%% model parameters %%%

sig = 0.2;
r = 0.05;
T = 0.5;
E = 100;
B = 120;
m = 6; % no of intervals
N = 8; % no of CGL points - 1
tol = 1e-3;
h0 = 1e-2;
opt = 'b'; % set option to barrier call
C = 0; % no penalty term
ep = 0;

%%% set up auxilliary variables for solving IVP %%%

[D1, D2] = lagrangeD(N);
Xi = xi(N); % the CGL points
smin = 0;
smax = B;
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
bact = uno_act(S, tt, E, r, sig, B, T);

%%% print needed values %%%
err = max(abs(u(end,:) - bact(end,:)))
accepted_steps = length(tt)-1
cpu_t

%%% error plots etc %%%
figure(1)
plot(S, abs(u(end,:) - bact(end,:)));
xlabel('S')
ylabel('Error(S)')
saveas(figure(1), "err_b.jpg");

figure(2)
plot(S, bact(end,:));
hold on;
scatter(S, u(end,:), 'black', 'o');
hold off;
xlabel('S')
ylabel('V(S,T)')
legend('Actual', 'Numerical', 'Location', 'northwest')
saveas(figure(2), "b_vs.jpg");