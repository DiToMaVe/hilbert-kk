function I1 = compute_integral_I1(integrand, ug, Nu, omega)

% Computation of random matrix [I1(omega)], Eq. (52) in Capillon,
% Desceliers & Soize (2016). The limit eta->0+ is taken into account by
% sampling the value of the integrand at the sampling points {ui,
% i=1,...,Nu, ui~=1} with a constant step << 1 of the interval [0, ug], and
% shifting the grid such that the singular point is exactly halfway in
% between two sampling points. The first and last point are then fixed to 0
% and ug respectively. As a result the step is constant except for the
% first and last interval.

du = ug/(Nu-1);     % stepsize (except for the first and last interval)
u = [0:Nu-1]*du;    % sampling points, before shifting the grid
idx_plus = ceil(Nu/ug)+1;
idx_min = idx_plus-1;
u_min = (idx_min-1)*du;

if u_min+0.5*du>1
    u_shift = (u_min+0.5*du)-1;
    u = u-u_shift;          % sampling points after shifting the grid
    u = [0,u(2:end),ug];    % Fixing start and end points
else
    u_shift = 1-(u_min+0.5*du);
    u = u+u_shift;          % sampling points after shifting the grid
    u = [0,u(1:end-1),ug];  % Fixing start and end points
end

% Trapezoidal integration rule
U_rep = repmat(u,length(omega),1); 
I1 = trapz(u,integrand(omega*u)./(1-U_rep.^2),2);
 
% I1_1 = trapz(u(:,1:idx_min),integrand(omega*u(1:idx_min))./(1-U_rep(:,1:idx_min).^2),2);
% I1_2 = trapz(u(:,idx_plus:Nu),integrand(omega*u(idx_plus:Nu))./(1-U_rep(:,idx_plus:Nu).^2),2);

end

