% function [] = postProc(Sol, solver, bc, timeStep, fileOthdId, ...
% 						wrkDir, problemString, crd, elemType, Flag3D)

tsStart = 1;
tsEnd = 100;

for tempo = tsStart:1:tsEnd
    saveFile = strcat('saveCaseData_',num2str(tempo),'.mat');
    load(saveFile);
    figHandle = figure;
    plot(Sol.p);
    ylim([0 1.5e5]);
    pause(0.1);

     title(sprintf('Cavitation Fraction at time %9.8f',timeStep*solver.dt), 'FontSize', 14)
      
     % Capture the plot as an image 
      frame = getframe(gca); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      
      % Write to the GIF File 
      if timeStep == tsStart 
          imwrite(imind,cm,'pres2.gif','gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,'pres2.gif','gif','WriteMode','append'); 
      end
    
    close;

end