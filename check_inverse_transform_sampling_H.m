clear all, close all, clc

% Check inverse transform sampling of H(eta), section 4.1.2 in Capillon,
% Desceliers & Soize (2016).

ug = 2;
mu = 0.5*log((ug-1)/(ug+1));

% Check pdf, cdf and inverse cdf
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta = linspace(-1,1,100);
pdf_H = @(eta) 2*ug./(mu*(eta.^2-2*eta+1-4*ug^2));      % Eq.(58), p.13

% Integrate pdf from eta=-1 to eta=1, see
% https://en.wikipedia.org/wiki/List_of_integrals_of_rational_functions,
% C = -1/(2*mu)*log((1+ug)/|1-ug|) = 1;
cdf_H = @ (eta) -1/mu*atanh((eta-1)/(2*ug))+1;  
inv_cdf_H = @(F) 1+2*ug*tanh(-mu*(F-1));

% % Alternative formulation
% cdf_eta = 1/(2*mu)*log(abs(((eta-1-2*ug)./(eta-1+2*ug))))+1;
% % plot integral
% figure()
% F1 = 1/(4*ug)*log(abs(((eta-1-2*ug)./(eta-1+2*ug))));
% F2 = -2/(4*ug)*atanh((eta-1)/(2*ug));
% plot(eta,F1,'b',eta,F2,'r')

% Plot pdf
figure()
plot(eta,pdf_H(eta),'b')

% Plot cdf
figure()
plot(eta,cdf_H(eta),'b')

% Plot cdf and inverse cdf
figure()
plot(eta,cdf_H(eta),'b',inv_cdf_H(cdf_H(eta)),cdf_H(eta),'r+')

% Check inverse transform sampling of H(eta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns = 500000;

% Sample from uniform distribution (uniform random variable)
urv = rand(1,Ns);
Eta = inv_cdf_H(urv);

% Compare histogram of sampled eta and pdf of eta
figure()
[hy, hx] = hist(Eta,50);            % get histogram values
hy = hy/numel(Eta)/(hx(2)-hx(1));   % normalize histogram
bar(hx, hy,'m')                     % plot histogram
hold on
plot(eta,pdf_H(eta),'LineWidth',2)  % plot pdf



