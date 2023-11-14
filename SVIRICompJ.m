function [J] = SVIRICompJ(t, y, omega_y, u, omega_u)

n = length(t);
T = max(t);

J = 1/(T)*trapz(t, sum(repmat(omega_y,1,n).*y) ...
                + sum(repmat(omega_u,1,n).*u.^2));
end