% Day1: hands-on session 
% Consider the Cauchy problem y'(t) = cos(2y(t)) t \in (0,1], y(0)=0
% The exact solution y = 1/2arcsin((e^{4t}-1)/(e^{4t}+1))
%
% Find the approximate solutions using forward Euler and backward Euler
% for h=1/2, 1/4, 1/8, ..., 1/512.
%
% Compute the error at the t=1 and store them in fe (forward Euler)
% be (backward Euler) and cn(Crank Nicolson).
%
% Estimate the order of convergence by 
% p_i = log(e_i / e_{i-1}))/log(h_i/h_{i-1})
% e_i errors relative to the discretization step h_i.
clear all; close all; clc;
tspan = [0,1]; y0=0; f=@(t,y) cos(2*y);
u=@(t) 0.5*asin( (exp(4*t)-1)./ (exp(4*t)+1));

Nh=2; 
for k=1:10
    [t,ufe] = feuler(f, tspan, y0, Nh);
    fe(k) = abs(ufe(end)- u(t(end)));
    [t,ube] = beuler(f, tspan, y0, Nh);
    be(k) = max(abs(ube- u(t)));
    [t,ucn] = cranknic(f, tspan, y0, Nh);
    cn(k) = max(abs(ucn- u(t)));
    Nh = 2*Nh;
end

pf = log(abs(fe(1:end-1) ./ fe(2:end)))/log(2);
pf(1:2:end)

pb = log(abs(be(1:end-1) ./ be(2:end)))/log(2);
pb(1:2:end)

pcn = log(abs(cn(1:end-1) ./ cn(2:end)))/log(2);
pcn(1:2:end)
%%
figure(1);
xaxis = (2.^(1:10));
plot(xaxis, fe,'r','LineWidth',1.5); hold on;
plot(xaxis, be,'g','LineWidth',1.5); hold on;
plot(xaxis, cn,'b','LineWidth',1.5); grid on;
legend('FE', 'BE', 'CN');
