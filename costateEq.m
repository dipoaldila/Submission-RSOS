function dz = costateEq(t,z,tu,u,s,w)
global betax mu p q xi theta epsilon k1 k2 d1 d2 

c = interp1(tu,u',t);
y = interp1(tu,s',t);
b=1;
d=1;
dz = zeros(5, 1);
dz(1) = p*(1-b*c(2)+b*c(2)*xi)*betax.*y(3).*(z(1)-z(2)) + q*(1-b*c(2)+b*c(2)*xi)*betax.*y(3).*(z(1)-z(3)) + mu*z(1);
dz(2) = theta*(z(2)-z(3)) + epsilon*(z(2)-z(4)) + mu*z(2);
dz(3) = -w(3) + p*(1-b*c(2)+b*c(2)*xi)*betax.*y(1).*(z(1)-z(2)) + q*(1-b*c(2)+b*c(2)*xi)*betax.*y(1).*(z(1)-z(3)) + d*c(1).*(z(3)-z(4)) + k1*(z(3)-z(5)) + (mu+d1)*z(3);
dz(4) = -w(4) + k2*(z(4)-z(5)) + (mu+d2)*z(4);
dz(5) = mu*z(5);

