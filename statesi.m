function dg = statesi(t, g, tu, u)
global delta betax mu p q xi theta epsilon k1 k2 d1 d2 

c = interp1(tu,u',t);

dg = zeros(5, 1);

dg(1) = delta - betax*g(1)*g(3) - mu*g(1);
dg(2) = p*betax*g(1)*g(3) - (theta+epsilon+mu)*g(2);
dg(3) = q*betax*g(1)*g(3) + theta*g(2) - (k1+d1+mu)*g(3);
dg(4) = epsilon*g(2) - (mu+k2+d2)*g(4);
dg(5) = k1*g(3) + k2*g(4) - mu*g(5);
