clear all, close all, clc

% Hilbert transform pair of a (real) even function.
% See https://en.wikipedia.org/wiki/Hilbert_transform
t = linspace(0,50,101);
% ft = @(t) 1./(t.^2+1);
% Ht = @ (t) t./(t.^2+1);
ft = @(t) cos(t);
Ht = @(t) sin(t)
% Paremeters of the numerical Hilbert transform computation
ug = 100;
Nq = 2;
Nu = 1000+1;
Ns = 640;

% Numerical computation of the Hilbert transform
[I1,u] = compute_integral_I1(ft, ug, Nu, t');

I2 = compute_integral_I2(ft, t', ug, Nq, Ns);

% K_num = 2*omega/pi.*(I1+I2);

% 
% figure(5)
% plot(freq,K_num./E_hat_imag)

%%
% Plot Hilbert transform pairs
figure(1)
plot(t,ft(t),'b',t,Ht(t),'r')
%%
figure(2)
plot(t,ft(t),'b')
hold on
plot(t,2/pi*I1,'r+')
hold on
plot(t,Ht(t),'k')



