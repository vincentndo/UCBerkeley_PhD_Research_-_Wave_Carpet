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

% READING FORT.15

fid=fopen('fort.15','r');
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

writerObj = VideoWriter(['3D Wave_Propagation_' num2str(nT) 'T.avi']);
writerObj.FrameRate = 200;
open(writerObj);

set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

x=0:2*pi/NDX:2*pi-2*pi/NDX;
y=0:2*pi/NDY:2*pi-2*pi/NDY;
[Y,X] = meshgrid(x,y);
for i = 1:total
    mesh(eta(:,:,i))
    hold on
    mesh(etab(:,:,i))
    hold off
%     surf(X,Y,eta(:,:,i))
    title(['3D Wave Propagation over ' num2str(nT) 'T'])
    view(0,0)
    caxis([-ETAmax*1.1 ETAmax*1.1])
    colorbar
    colormap(jet)
%     axis([0 2*pi 0 2*pi -ETAmax*1.05 ETAmax*1.05])
    axis([-50 NDX+50 -50 NDY+50 -ETAmax*1.1 ETAmax*1.1])
    frame(i) = getframe(gcf);
    writeVideo(writerObj,frame(i));
end

close(writerObj);

