clc; clear; close all;

%Position = csvread('Position'); Velocity = csvread('Velocity');

load('plotting.mat')

for i= 1:size(x,1)
scatter(x(i,:),zeros(1,size(x,2)),500,abs(v(i,:)),'filled');
colormap jet;
colorbar;
caxis([0 0.5])
xlim([0 5])
% title('Particles colored with velocity, time step = 0.04, Verlet Algorithm')
pause(0.00001);

% Uncomment the following to save

% drawnow;
% frame = getframe(1);
% im = frame2im(frame);
% [imind,cm] = rgb2ind(im,256);
% outfile = 'SevenParticles.gif';
% 
% % On the first loop, create the file. In subsequent loops, append.
% if i==1
%     imwrite(imind,cm,outfile,'gif','DelayTime',0,'loopcount',inf);
% else
%     imwrite(imind,cm,outfile,'gif','DelayTime',0,'writemode','append');
% end

end