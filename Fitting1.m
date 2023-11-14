%function parest3
clear all; clc; close all;

%% data procescosg
%datatotal=xlsread('smoothing1718.xls',1);
%datatotal="datatotal";
load TBIndo.csv; % calling data
%startDate = year('2000'); % starting date on wish
%endDate = year('2020');
%beginDate = year('01-01-2018'); % starting date from JHU
%cut = length(beginDate:startDate) - 1;
%tplot = startDate:endDate;
%len = length(tplot);
%datatotal(1:cut,:) = []; datatotal = datatotal(1:len); % clear data
%% known parameters
% q (parameter yang diketahui)
% p (Yang dicari)
% y (Urutan variabelnya)
q(1) = 100000/72;
q(2)= 1/72; q(3)=0.5; q(4)=0.5; q(5)=1; q(6)=0.01; q(7)=0.01;

Sigma = 1e-2*mean(TBIndo);%Tidak usah dirubah

%% solver
% p(beta,c1,c2)
%options = optimoptions('fmincon','Algorithm',...
           %'interior-point','Display','iter');
options = optimoptions('fmincon','Algorithm',...
           'interior-point','Display','iter','MaxFunEvals',1e+5,'MaxIter',700);
Ain = []; bin = []; Aeq = []; beq = [];
%lb = [0.08,0.04,0,0,1/400]; 
%ub = [0.12,0.06,0.08,0.04,1/300];
%p (chi1, varpi, delta2, Su(0), H(0), A(0), gamma1)
lb = [ 0, 0, 0, 10000, 100, 1000, 1000, 0]; 
ub = [ 10/100000, 1, 1/3, 99000, 10000, 7000, 7000, 700 ];
p0 = lb + (ub-lb).*rand(1,length(lb)); % initial guess
%p0 = [2.99999999971353,0.333333333248506,0.00999999999752192,10063082.4601011,245.597927149500,999999.999820356];
p = fmincon(@(p) Obj(p,q,TBIndo,Sigma),p0,Ain,bin,...
              Aeq,beq,lb,ub,@(p) confun(p,q),options);
%writematrix(M1,'smoothing2021.xls')

%% evaluate
ti = linspace(1,length(TBIndo),length(TBIndo));
tspan = [ti(1),ti(end)];
Y = deval(ode45(@(t,y) model(t,y,p,q), tspan, p(4:8)),ti);

%% plotting
plot(ti,TBIndo,'ko-'); hold on
plot(ti,p(3)*Y(2,:),'bo-');
set(gca,'FontSize',16,'XGrid','on')
legend('data','model output');
xlabel('data point (years)')
ylabel('new incidence')

%datetick('x','mmm dd','keepticks');
%set(gcf,'position',[500 500 1300 500])
grid minor


% input: optimal parameters p, known parameters q
gain = 1e-10; P = length(p); par = p;
for i = 1:P
disp(['iter = ',num2str(i)])
for j = 1:P
  if j==i
    pp = p; pp(j) = par(j)+gain*par(j); pm = par; ...
    pm(j) = par(j)-gain*par(j);
    hessian(i,j) = (Obj(pp,q,TBIndo,Sigma)...
    - 2*Obj(par,q,TBIndo,Sigma) + Obj(pm,q,TBIndo,Sigma))...
    /(gain^2*p(j)^2);
  else
    epsp = gain*double(logical([1:P]==i | [1:P]==j));
    epsi = epsp; epsi(j) = -epsi(j);
    epsj = epsp; epsj(i) = -epsj(i);
    hessian(i,j) = (Obj(par+epsp.*par,q,TBIndo,Sigma)...
    - Obj(par+epsi.*par,q,TBIndo,Sigma)...
    - Obj(par+epsj.*par,q,TBIndo,Sigma)...
    + Obj(par-epsp.*par,q,TBIndo,Sigma))...
    /(4*gain^2*p(i)*p(j));
end
end
end

% confidence interval, alpha = 0.05;
CI = 1.96*sqrt(diag(abs(inv(hessian))));

%end

%% Define the model
function dy = model(t,y,p,q)
% p()
% q ()
% y ()
dy = zeros(5,1);
dy(1) = q(1) - p(1)*y(1)*y(3) - q(2)*y(1);
dy(2) = p(2)*p(1)*y(1)*y(3) - (q(3)+p(3)+q(2))*y(2);
dy(3) = (1-p(2))*p(1)*y(1)*y(3) + q(3)*y(2) - (q(4)+q(6)+q(2))*y(3);
dy(4) = p(3)*y(2) - (q(2)+q(5)+q(7))*y(4);
dy(5) = q(4)*y(3) + q(5)*y(4) - q(2)*y(5);
end


function LL = Obj(p,q,TBIndo,Sigma)
ti = linspace(1,length(TBIndo),length(TBIndo));
tspan = [ti(1),ti(end)];
y0 = p(4:8); 
% solution on a daily basis
Y = deval(ode45(@(t,y) model(t,y,p,q), tspan, y0),ti);
%G = p(2)*Y(1)*Y(3)+((p(3)*p(2))/(1+p(4)*Y(3)))*Y(2)*Y(3);
LL = 0.5*sum((p(3)*Y(2,:)- TBIndo').^2)/Sigma^2;
%LL = 0.5*sum((Y(3,:)- dataHIV1').^2)/Sigma^2;
end

function [cin,ceq] = confun(p,q)
ceq = []; % make = 0
cin = []; % make <= 0
end