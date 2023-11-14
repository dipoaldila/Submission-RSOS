clear all ; close all ; clc ; format long ;

T  = 50;%duration of simulation --> Waktu maksimal simulasi
nT = 51; %Banyak titik diskrit
ti = linspace(0,T,nT)';

global y0 delta betax mu p q xi theta epsilon k1 k2 d1 d2 

%Deklarasi nilai parameter
betax=2.6016*10^(-5);
p =0.9013;
q=1-p;
epsilon=0.2969;

delta = 100000/70; 
mu=1/70;
xi=0.1;
theta=1/2;
k1=1/2;
k2=1;
d1=0.01;
d2=0.01;


u      = zeros(2,nT);%Initial guess 
u(1,:) = 0.1*ones(1,nT);  
u(2,:) = 0.1*ones(1,nT);  

omega_y = zeros(5,1);
omega_y(3) = 10^(-5);
omega_y(4) = 10^(-5);

omega_u = zeros(2,1);
omega_u(1) = 0.1;
omega_u(2) = 0.05;

%y0 = [91513; 1746; 1480; 2671; 385];%initial condition for y\
%Angola
%y0 = [6.835652571588624
 %  0.188672956763733
 %  0.191629976130740
 %  0.032373069198732
 %  2.543531042679899]*10^4;
%Lesotho
%y0 = [4.094050369341764
 %  0.263195648273222
 %  0.318351687474206
 %  0.049188401734191
 %  5.822849027241827]*10^4;
% India
%y0 = [4.330971638659708
 %  0.126637363726976
  % 0.129800150334217
   %0.021319622827213
   %2.211101071810003]*10^4;
%Indonesia
y0 = [3.009873836043487
   0.106792015107841
   0.121646404106515
   0.031191571298874
    2.067846443523603]*10^4;

z0 = zeros(1,5);%final costate 


Ops = odeset('RelTol', 1E-6);

km =29;
Jv = zeros(1,km);

fprintf('\n***');

g = deval(ode45(@(s,x) statesi(s,x,ti,u), [0 T], y0, Ops), ti);%Tanpa Kontrol
 r = deval(ode45(@(s,x) statesin(s,x,ti,u), [0 T], y0, Ops), ti); % Dengan Kontrol konstan berupa tebakan awal
[J0] = SVIRICompJ(ti,g,omega_y,u,omega_u)
 [J] = SVIRICompJ(ti,r,omega_y,u,omega_u)
 yi    = r;
 ui    = u;
Jv(1) = J;

%Proses iterasi
for k=2:km
z = deval(ode45(@(s,x) costateEq(s,x,ti,u,r,omega_y), [T 0], z0, Ops), ti);%Hitung variabel Adjoin
%plot(ti,z(1,:))
un = SVIRIcompU(r,z,omega_u);%Hitung Kontrol 
       u1(1,:) = min(1,max(0,un(1,:)));%Update kontrol
       u1(2,:) = min(0,max(0,un(2,:)));%Update kontrol
u_opt(1,:) = (u1(1,:)+u(1,:))/2;%Percepat konvergensi
u_opt(2,:) = (u1(2,:)+u(2,:))/2;%Percepat konvergensi
y_opt = deval(ode45(@(s,x) statesin(s,x,ti,u_opt), [0 T], y0, Ops), ti);%Update variable state/manusia
[J_opt] = SVIRICompJ(ti, y_opt, omega_y, u_opt, omega_u)
     
y     = y_opt;
u     = u_opt;
J     = J_opt;
figure(1)
subplot(1,2,1)
  plot(ti,u(1,:))
  hold on
  pause(0.0001)
  xlabel('years');
  ylabel('u1');

  subplot(1,2,2)
  plot(ti,u(2,:))
  hold on
  pause(0.0001)
  xlabel('years');
  ylabel('u2');
end

Invaver = trapz(ti,((g(2,:)-y_opt(2,:))+(g(3,:)-y_opt(3,:)) + (g(4,:)-y_opt(4,:))))
Cost = trapz(ti,((u_opt(2,:).*y_opt(3,:))+u_opt(1,:).*(y_opt(2,:))))
Recov = trapz(ti,y_opt(5,:))
IAR = Invaver/Recov
ACER = Cost/Invaver



figure(2)
plot(ti,g(1,:),'bo-',ti,y_opt(1,:),'ro-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('S','Fontsize',16);
legend('without control','with control', 'Fontsize',14);
%title('(a)','Fontsize',14)
grid on

figure(3)
plot(ti,g(2,:),'bo-',ti,y_opt(2,:),'ro-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('E','Fontsize',16);
legend('without control','with control', 'Fontsize',14);
%title('(b)','Fontsize',14)
grid on

figure(4)
plot(ti,g(3,:),'bo-',ti,y_opt(3,:),'ro-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('I_1','Fontsize',16);
legend('without control','with control', 'Fontsize',14);
%title('(c)','Fontsize',14)
grid on

figure(5)
plot(ti,g(4,:),'bo-',ti,y_opt(4,:),'ro-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('I_2','Fontsize',16);
legend('without control','with control', 'Fontsize',14);
%title('(d)','Fontsize',14)
grid on

figure(6)
plot(ti,g(5,:),'bo-',ti,y_opt(5,:),'ro-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('R','Fontsize',16);
legend('without control','with control', 'Fontsize',14);
%title('(e)','Fontsize',14)
grid on


figure(7)
stem(ti,epsilon*y_opt(2,:)+u_opt(2,:).*y_opt(3,:),'filled','ro')
xlim([0 51])
 xlabel('years','Fontsize',16)
 ylabel('New Incidence per 100 000 people','Fontsize',16)
 legend('New detected case per 100 000 individual','Fontsize',14)
 grid on

figure(8)
plot(ti,u_opt(1,:),'mo-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('u_1','Fontsize',16);
legend('u_1(t)','Fontsize',14)
grid on

figure(9)
plot(ti,u_opt(2,:),'co-','LineWidth',1);
xlabel('years','Fontsize',16);
ylabel('u_2','Fontsize',16);
legend('u_2(t)','Fontsize',14)
grid on
