clear all;          close all;                 clc;
g=1; rho = 1;
M=7;            NDX=2^(M+2);     MY=0;       NDY=2^(MY+2);     % *IMPORTANT

N=20;   % NUMBER OF PERIOD 
IPER=512;    % STEP PER PERIOD: if IPER <= 64 -> unstable
ISTP = N*IPER;
NDP=1;   % number of modes

lambda=2*pi*0.05;
kx=2*pi/lambda;         ky=0;                          % *IMPORTANT
k = sqrt(kx^2+ky^2);    as=0.02/k;
h = 0.01;
% omega=(k*g*tanh(k*h))^(1/2);
% T=2*pi/omega;
% dT = T/IPER;
zeta = 0.0;
gamma = 0.01;
mu = k*h;

p = [gamma*tanh(mu), 1i*mu*gamma*zeta, -mu, ...
    -1i*mu^2*gamma*zeta*tanh(mu), mu^2*(1-gamma)*tanh(mu)];
OMEGAroots = roots(p);                            
% OMEGA = OMEGAroots(1)
OMEGA = OMEGAroots(4)                                          % *IMPORTANT
omega = OMEGA/sqrt(h/g);
T=2*pi/real(omega);
dT = T/IPER;

ab = as*cosh(mu)*(1-mu*tanh(mu)/OMEGA^2);
A = -1i*as*(omega^2+g*k)/(2*k*omega);
B = 1i*as*(omega^2-g*k)/(2*k*omega);

x=linspace(0,2*pi-2*pi/NDX,NDX);%0:2*pi/kx/NDX:((2*pi-2*pi/NDX)/kx);
y=linspace(0,2*pi-2*pi/NDY,NDY);%0:2*pi/ky/NDY:((2*pi-2*pi/NDY)/ky);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

eta = as*exp(1i*kx*X+1i*ky*Y);
etab = ab*exp(1i*kx*X+1i*ky*Y);
phi = (A*exp(k*0)+B*exp(-k*0)).*exp(1i*kx*X+1i*ky*Y);
phib = (A*exp(k*(-h+0))+B*exp(-k*(-h+0))).*exp(1i*kx*X+1i*ky*Y);

eta_over_time = @(t) as*exp(1i*(kx*X+ky*Y-omega*t));
etab_over_time = @(t) ab*exp(1i*(kx*X+ky*Y-omega*t));

% figure(10)
% for t = linspace(0,4*T-dT,512)
%     ETA = eta_over_time(t);
%     ETAb = etab_over_time(t);
% %     plot(real(ETA(1,:)));    
%     mesh(X,Y,real(ETA))
%     hold on
%     mesh(X,Y,real(ETAb))
%     hold off
%     colorbar
%     colormap(jet)
%     xlabel('x')
% %     axis([0 nnod 0 nlayer -as-as/6 as+as/6])
%     axis([0 2*pi/kx 0 2*pi/ky -as-as/6 as+as/6])
%     daspect([max(daspect)*[1 1] 1])
% %     view(0,0)
%     caxis([-as as])
%     zlim([-as-as/6 as+as/6])
%     pause(0.001)    
% end

% kb=k*2; lambdab=2*pi/kb; L0=10*lambdab; d=0.31/kb; 

figure(1)
mesh(Y,X,real(eta))
colormap(jet)
title('Surface')
figure(2)
mesh(Y,X,real(etab))
colormap(jet)
title('Bottom')

% WRITING FORT.10
fid=fopen('fort.10','w');
fprintf(fid,'%g %g %g\n',NDP,M,MY);   % # of modes, Np, Nq

% TIME STEP PER PERIOD AND TOTAL TIME STEPS
fprintf(fid,'%g %g\n',IPER,ISTP);

% DIMENSIONS, PERIOD AND DEPTH
fprintf(fid,'%18.15e %18.15e %18.15e %g\n',2*pi,2*pi,T,h);
fprintf(fid,'%18.15e %18.15e %18.15e\n',zeta,gamma,mu);

for j=1:1:NDY
    for i = 1:1:NDX        
        fprintf(fid,'%18.15e\n',real(eta(i,j)),real(phi(i,j)));
    end
end
fclose(fid);

fid=fopen('fort.9','w');
for j=1:1:NDY
    for i = 1:1:NDX
        fprintf(fid,'%18.15e\n',real(etab(i,j)),real(phib(i,j)));
    end
end
fclose(fid);

% step_range =[128, 256, 512, 640, 768, 896, 1024];
% error_s_range = [2.057200569473616e-09, -6.343687688912300e-11, 1.986321986520897e-12, 6.765859145170705e-13, 2.573173156032264e-13, -1.188907190289671e-13, -6.097914418782639e-14];
% error_b_range = [1.406048747348626e-09, 4.335817252559161e-11, -1.333820207061187e-12, -4.142001563303678e-13, -1.758547011713555e-13, 8.072029734444246e-14, 4.111728318934027e-14];
% 
% figure()
% semilogy(T./step_range,abs(error_s_range),'*-','LineWidth',2)
% hold on
% semilogy(T./step_range,abs(error_b_range),'r*-','LineWidth',2)
% hold off
% title('Error Analysis')
% xlabel('time step [s]')
% ylabel('relative error')
% legend('surface','bottom')
% grid on
% 
% figure()
% semilogy(step_range,abs(error_s_range),'*-','LineWidth',2)
% hold on
% semilogy(step_range,abs(error_b_range),'r*-','LineWidth',2)
% hold off
% title('Error Analysis')
% xlabel('T/dt')
% ylabel('relative error')
% legend('surface','bottom')
% grid on

% -----------------------------------

step_range =[128, 256, 512, 640, 768, 896, 1024];
error_s_range = [2.213741989405472e-09, 6.589926345884301e-11, 2.063687326545430e-12, -6.396951833983108e-13, 2.730383510278829e-13, 1.264952973769365e-13, 6.617505102249763e-14];
error_b_range = [1.356530817344740e-09, -4.034820737210304e-11, 1.263509677546156e-12, -4.373035498532000e-13, 1.669277945414615e-13, -7.722846969050421e-14, -4.041229424986642e-14];

figure()
semilogy(T./step_range,abs(error_s_range),'*-','LineWidth',2)
hold on
semilogy(T./step_range,abs(error_b_range),'r*-','LineWidth',2)
hold off
title('Error Analysis. Skew Wave')
xlabel('time step [s]')
ylabel('relative error')
legend('surface','bottom')
grid on

figure()
semilogy(step_range,abs(error_s_range),'*-','LineWidth',2)
hold on
semilogy(step_range,abs(error_b_range),'r*-','LineWidth',2)
hold off
title('Error Analysis. Skew Wave')
xlabel('T/dt')
ylabel('relative error')
legend('surface','bottom')
grid on