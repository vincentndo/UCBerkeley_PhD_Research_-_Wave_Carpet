close ALL;     clear;

% get PARAMETERS
% READING FORT.10
fid=fopen('fort.10','r');
line1 = fgetl(fid);
l1 = textscan(line1,'%f');
NDX = 2^(l1{1}(2)+2);
NDY = 2^(l1{1}(3)+2);

line2 = fgetl(fid);
l2 = textscan(line2,'%f');
nstep = l2{1}(1);
total = l2{1}(2);
nT = total/nstep;
fclose(fid);

% READING FORT.11
fid_eta=fopen('fort.15','r');
% fid_phi=fopen('fort.16','r');
STORE_eta = fscanf(fid_eta,'%f');
% STORE_phi = fscanf(fid_phi,'%f');

eta = zeros(NDY,NDX,total);
etab = zeros(NDY,NDX,total);
% phi = zeros(NDY,NDX,total);

lnum = 0;
for k =1:total
    for j = 1:NDY
        lnum = lnum+1;
%         eta(j,:,k) = STORE_eta((1:NDX)+NDX*(lnum-1));
        eta(j,:,k) = STORE_eta((1:NDX)+NDX*(2*lnum-2));
        etab(j,:,k) = STORE_eta((1:NDX)+NDX*(2*lnum-1));
%         phi(j,:,k) = STORE_phi((1:NDX)+NDX*(lnum-1));
    end
end
fclose(fid_eta);
% fclose(fid_phi);

domain = 1;
x=0:2*pi*domain/NDX:(2*pi*domain-2*pi*domain/NDX);

% f1 = (tanh(20*x-5) + 1)/2;
% f2 = (-tanh(20*x-15) + 1)/2;
% f = [f1(x<0.5) f2(x>=0.5)];
% 
% g1 = (tanh(20*x-15)+1)/2;
% g2 = (-tanh(20*x-2*pi*20+5)+1)/2;
% g = [g1(x<pi+0.25) g2(x>=pi+0.25)];

% f1 = (tanh(5*x-0) + 0);
% f2 = (-tanh(5*x-5) + 1)/2;
% f = [f1(x<0.5) f2(x>=0.5)];
% f = (tanh(40*(x-0.05*pi))+tanh(-20*(x-0.25*pi)))/2;

% g1 = (tanh(5*x-5)+1)/2;
% g2 = (-tanh(5*x-2*pi*5+0)+0)/1;
% g = [g1(x<pi+0.5) g2(x>=pi+0.5)];
% g = (tanh(20*(x-0.25*pi))+tanh(-10*(x-(2-0.05)*pi)))/2;

f = (tanh(40*(x-0.05*pi))+tanh(-20*(x-0.25*pi)))/2;
g = (tanh(20*(x-0.25*pi))+tanh(-10*(x-(2-0.05)*pi)))/2;

figure(4)
plot(x,f,'b',x,g,'r','LineWidth',1.5)
title('Filter f')
figure(5)
plot(x,g,'r','LineWidth',1.5)
title('Filter g')

f_mat = repmat(f,NDY,1);
g_mat = repmat(g,NDY,1);

eta_new = zeros(NDY,NDX,total);
etab_new = zeros(NDY,NDX,total);
% phi_new = zeros(NDY,NDX,total);

for k = 1:total
    eta_new(:,:,k) = eta(:,:,k).*f_mat;
    etab_new(:,:,k) = etab(:,:,k).*f_mat;
%     phi_new(:,:,k) = phi(:,:,k).*f_mat;
end

figure()
mesh(eta_new(:,:,1))

figure()
mesh(eta(:,:,1).*g_mat)

% eta_new_2 = eta.*g_mat;

% WRITING FORT.20
fid=fopen('fort.20','w');

% for i=1:1
%     for j=1:NDX
%         fprintf(fid,'%18.15e %18.15e\n',real(eta_new(i,j)),imag(eta_new(i,j)));
%         fprintf(fid,'%18.15e %18.15e\n',real(phi_new(i,j)),imag(phi_new(i,j)));
%     end
% end

for k=1:total
    for j=1:NDY
        for i = 1:NDX
%             fprintf(fid,'%18.15e %18.15e\n',eta_new(j,i,k),phi_new(j,i,k));
            fprintf(fid,'%18.15e %18.15e\n',eta_new(j,i,k),etab_new(j,i,k));
        end
    end
end

for i=1:NDX
    for j = 1:NDY
        fprintf(fid,'%18.15e %18.15e\n',0,0);
%         fprintf(fid,'%18.15e\n',0);
    end
end

fclose(fid);

% WRITING FORT.23
fid = fopen('fort.23','w');

for i = 1:NDX
    fprintf(fid,'%18.15e %18.15e\n',f(i),g(i));
end
fclose(fid);