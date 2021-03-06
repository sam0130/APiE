clc; clear; close all;

%Position = csvread('Position'); Velocity = csvread('Velocity');

load('plotting.mat')
for i= 1:size(position.x,1)
    % hold on;
    scatter3(position.x(i,:),position.y(i,:), zeros(1,size(position.x,2)),700, ...
        abs(sqrt(velocity.x(i,:).^2 + velocity.y(i,:)).^2),'filled');
    view(2)
    set(gcf, 'Position', get(0, 'Screensize'));
    colormap jet;
    colorbar;
    caxis([0 0.5])
    xlim([-0.5 Np*re+0.5])
    ylim([-0.5 Np*re+0.5])
    % plot([0,0], [-1,1], 'g--','Linewidth',3);
    % plot([2,2], [-1,1], 'g--','Linewidth',3);
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