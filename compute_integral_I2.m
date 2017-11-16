function I2 = compute_integral_I2(integrand, omega, ug, Nq, Ns)

% Computation of random matrix [I2(omega)], Eq. (54), as explained in
% section 4.1.2 of Capillon, Desceliers & Soize (2016). Remark: the factor
% mu is lacking in Eq. (60) (it does appear in Eq. (56)), I believe this is
% an error.

mu = 0.5*log((ug-1)/(ug+1));

% Inverse transform sampling of H(eta)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample from uniform distribution (uniform random variable)
urv = rand(1,Ns);
eta_smpl = 1+2*ug*tanh(-mu*(urv-1));

% Generation of the orthogonal polynomials and computation of the
% coefficients alpha and beta
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute M, the (Nq x Ns) real matrix of realizations of the monomials
eta_smpl_rep = repmat(eta_smpl,Nq,1);
kappa = [1:Nq]';
H_exp = repmat(kappa-1,1,Ns);

M = eta_smpl_rep.^H_exp;

% Compute the (Nq x Nq) real matrix F. It is assumed that ns>Nq in order 
% that matrix F is positive definite
F = 1/(Ns-1)*(M*M');

% Compute the lower triangular (Nq x Nq) real matrix L from the Cholesky
% decomposition LL' of positive-definite matrix F
L = chol(F,'lower'); 

% Compute the (Nq x Ns) real matrix P as the solution of the linear matrix
% equation LP=M
P = L\M;

% Calculate coefficients alpha_k and beta_k for k=0,...,Nq-1, Eq. (80)-(81)

% alpha_k for k=0,...,Nq-1
Alpha = 1/Ns*(P.^2)*eta_smpl';

% beta_k for k=1,...,Nq 
PP = P(1:Nq-1,:).*P(2:Nq,:);
Beta_sqrt = 1/Ns*PP*eta_smpl';

% Construction of the quadrature rule for the computation of I2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the Jacobi matrix
Jac = diag(Alpha)+diag(Beta_sqrt,-1)+diag(Beta_sqrt,1);

% Compute eigenvalues eta_q and eigenvectors v_q of the Jacobi matrix
[V_q,Eta_q] = eig(Jac);
Eta_q = diag(Eta_q);

% The desired abscissa are eta_1,...,eta_Nq and the associated weights
% w_1,...,w_Nq are computed as w_q = (v_q)^{2}_{1} where v_q_{1} is the
% first component of the q-th normalized eigenvector v_q.
W_q = V_q(1,:).^2;

% Change of variable U(eta)=2ug/(1-eta), Eq. (55)
U_cal = @(x) 2*ug./(1-x);

% Compute integral I2
I2 = W_q*integrand(U_cal(Eta_q)*omega');
I2 = mu*I2'; % factor mu not in the paper???

end






















