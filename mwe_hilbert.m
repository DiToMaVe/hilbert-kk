clear all, close all, clc

% Minimimal working example illustrating the numerical computation of the
% Hilbert transform as explained in sections 4.1.1 and 4.1.2 of Capillon,
% Desceliers & Soize (2016).

% Hilbert transform pair of a (real) even function.
% See https://en.wikipedia.org/wiki/Hilbert_transform
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = linspace(0,50,101);
 
ft = @(t) 1./(t.^2+1);
Ht = @ (t) t./(t.^2+1);
 
% ft = @(t) cos(t);
% Ht = @(t) sin(t);

% Paremeters of the numerical Hilbert transform computation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ug = 100;  % ug>1, separation between trapezoidal integration rule, 
            % [I1(omega)], and Gaussian quadrature rule, [I2(omega)]. 
Nq = 20;
Nu = 1000+1;
Ns = 100;

% Numerical computation of the Hilbert transform
[I1,u] = compute_integral_I1(ft, ug, Nu, t');
I2 = compute_integral_I2(ft, t', ug, Nq, Ns);
I = I1+I2;

% Plot Hilbert transform pair: analytical expression and numerical solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Hilbert transform pairs
figure(1)
plot(t,ft(t),'b',t,Ht(t),'r')

% Error between analytical and numerical solutions 
err = Ht(t)-2/pi*I';
max_err = max(abs(err))

% Plot results of numerical computation of Hilbert transform and compare
% with analytical expression
figure(2)
plot(t,2/pi*I1,'b')
hold on 
plot(t,2/pi*I2,'r')
hold on
plot(t,Ht(t),'k')
hold on
plot(t,2/pi*I,'m+')
hold on 
plot(t,err,'g')
legend('2/pi I1','2/pi I2','H(t)','2/pi I')
