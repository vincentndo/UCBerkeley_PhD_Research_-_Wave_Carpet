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

% READING FORT.25

fid=fopen('fort.25','r');
STORE = fscanf(fid, '%f');
ETAmax = max(STORE(1:NDX));

eta = zeros(NDY, NDX, total);
etab = zeros(NDY, NDX, total);
lnum = 0;
for k =1:total
    for i = 1:NDY
        lnum = lnum+1;
        eta(i,:,k) = STORE((1:NDX)+NDX*(2*lnum-2));
        etab(i,:,k) = STORE((1:NDX)+NDX*(2*lnum-1));
    end
end
fclose(fid);

ETAbmin = -max(max(etab(:,:,1)));
for k =2:total
    ETAbmin_new = -max(max(etab(:,:,k)));
    if ETAbmin_new < ETAbmin
        ETAbmin = ETAbmin_new;
    end
end

writerObj = VideoWriter(['3D Wave_Maker_Elevation_' num2str(nT) 'T.avi']);
writerObj.FrameRate = 200;
open(writerObj);

set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

x=0:2*pi/NDX:2*pi-2*pi/NDX;
y=0:2*pi/NDY:2*pi-2*pi/NDY;
[X,Y] = meshgrid(x,y);
for i = 1:total
    mesh(X,Y,eta(:,:,i))
    hold on
    mesh(X,Y,etab(:,:,i) - 0.01)
    hold off
%     surf(X,Y,eta(:,:,i))
    title(['3D Wave Maker Elevation over ' num2str(nT) 'T. Timestep=' num2str(i)])
%     view(0,0)
%     caxis([-ETAmax*1.1 ETAmax*1.1])
    caxis([-ETAmax*1.1 ETAmax*1.1])
    colorbar
    colormap(jet)
%     axis([0 2*pi 0 2*pi -ETAmax*1.05 ETAmax*1.05])
    axis([0 2*pi 0 2*pi (ETAbmin-0.01)*1.05 ETAmax*1.4])
%     axis([-50 NDX+50 -50 NDY+50 -ETAmax*1.1 ETAmax*1.3])
%     axis([-50 NDX+50 -50 NDY+50 (ETAbmin-0.05*0.004)*1.05 ETAmax*1.4])
    frame(i) = getframe(gcf);
    writeVideo(writerObj,frame(i));
end

close(writerObj);

