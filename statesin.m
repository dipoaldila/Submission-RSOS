function dr = statesin(t, r, tu, u)
global delta betax mu p q xi theta epsilon k1 k2 d1 d2 

c = interp1(tu,u',t);
b=0;
d=1;
dr = zeros(5, 1);

dr(1) = delta - (1-b*c(2)+b*c(2)*xi)*betax*r(1)*r(3) - mu*r(1);
dr(2) = (1-b*c(2)+b*c(2)*xi)*p*betax*r(1)*r(3) - (theta+epsilon+mu)*r(2);
dr(3) = (1-b*c(2)+b*c(2)*xi)*q*betax*r(1)*r(3) + theta*r(2) - (d*c(1)+k1+d1+mu)*r(3);
dr(4) = epsilon*r(2) + d*c(1)*r(3) - (mu+k2+d2)*r(4);
dr(5) = k1*r(3) + k2*r(4) - mu*r(5);