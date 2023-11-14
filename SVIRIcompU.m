function u = SVIRIcompU(y,z,omega_u)

global betax xi  p q ;

n      = size(y,2);
u      = zeros(2,n);

u(1,:) = 1*y(3,:).*(z(3,:)-z(4,:))/(2*omega_u(1));
u(2,:) = 0*-(betax*y(1,:).*y(3,:)*(1-xi).*(z(1,:)-p*z(2,:)-q*z(3,:)))/(2*omega_u(2));
end

%g(1,:) = 2*omega_u(1)*u(1,:) + (y(1,:).*(z(2,:)-z(1,:)) + y(3,:).*(z(4,:)-z(3,:)));
%g(2,:) = 2*omega_u(2)*u(2,:) - (y(2,:).*(z(2,:)-z(1,:)) + y(4,:).*(z(4,:)-z(3,:)));