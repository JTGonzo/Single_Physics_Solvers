% function [] = postProcCreateGif(solver, rangeStart, rangeEnd)
function [] = postProcCreateGif2(outFreq, dTime, rangeStart, rangeEnd)

% 	for i = rangeStart:solver.outFreq:rangeEnd
    for i = rangeStart:outFreq:rangeEnd
		
		% x-velocity		
			file = strcat('xVelPlot_',num2str(i),'.fig');
			openfig(file);
            view(0,90), axis equal, axis on;
            hold on;
            xlim([-0.005*sqrt(2) 0.005*sqrt(2)]);
			ylim([-0.005*sqrt(2) 0.005*sqrt(2)]);
			zlim([-4 4]);
            axis tight manual;            
			
% 			title(sprintf('x-Velocity at time %9.8f',i*solver.dt), 'FontSize', 14);
            title(sprintf('x-Velocity at time %9.8f',i*dTime), 'FontSize', 14);
			set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
			colormap('jet');
			colorbar;           
			frame = getframe(gcf); 
			im = frame2im(frame); 
			[imind,cm] = rgb2ind(im,256); 
            
			% Write to the GIF File 
			if i == rangeStart 
				imwrite(imind,cm,'anim_xVelPlot.gif','gif', 'Loopcount',inf); 
			else 
				imwrite(imind,cm,'anim_xVelPlot.gif','gif','WriteMode','append'); 
            end     
    end
	close();	
            
    for i = rangeStart:outFreq:rangeEnd        

		% y-velocity		
			file = strcat('yVelPlot_',num2str(i),'.fig');
			openfig(file);
            view(0,90), axis equal, axis on;
            hold on;            
			xlim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			ylim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			zlim([-4 4]);
            axis tight manual;            
			
% 			title(sprintf('y-Velocity at time %9.8f',i*solver.dt), 'FontSize', 14);
			title(sprintf('y-Velocity at time %9.8f',i*dTime), 'FontSize', 14);
			set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
			colormap('jet');
			colorbar;
			frame = getframe(gcf); 
			im = frame2im(frame); 
			[imind,cm] = rgb2ind(im,256); 
            
			% Write to the GIF File 
			if i == rangeStart 
				imwrite(imind,cm,'anim_yVelPlot.gif','gif', 'Loopcount',inf); 
			else 
				imwrite(imind,cm,'anim_yVelPlot.gif','gif','WriteMode','append'); 
            end
    end
	close();
            
    for i = rangeStart:outFreq:rangeEnd        

		% pressure	
			file = strcat('pressurePlot_',num2str(i),'.fig');
			openfig(file);
            view(0,90), axis equal, axis on; 
            hold on;            
			xlim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			ylim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			zlim([0 4.5e5]);
            axis tight manual;            
			
% 			title(sprintf('Pressure at time %9.8f',i*solver.dt), 'FontSize', 14);
			title(sprintf('Pressure at time %9.8f',i*dTime), 'FontSize', 14);
			set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
			colormap('jet');
			colorbar;
			frame = getframe(gcf); 
			im = frame2im(frame); 
			[imind,cm] = rgb2ind(im,256); 
			
            % Write to the GIF File 
			if i == rangeStart 
				imwrite(imind,cm,'anim_pressure.gif','gif', 'Loopcount',inf); 
			else 
				imwrite(imind,cm,'anim_pressure.gif','gif','WriteMode','append'); 
            end
    end
	close();			
		
	for i = rangeStart:outFreq:rangeEnd
            
		% cavitation fraction	
			file = strcat('cavFracPlot_',num2str(i),'.fig');
			openfig(file);
			view(0,90), axis equal, axis on;  
            hold on;            
			xlim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			ylim([-0.5*sqrt(2) 0.5*sqrt(2)]);
			zlim([0 1]);
            axis tight manual;            
			
% 			title(sprintf('Cavitation Fraction at time %9.8f',i*solver.dt), 'FontSize', 14);
			title(sprintf('Cavitation Fraction at time %9.8f',i*dTime), 'FontSize', 14);
			set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
			colormap('jet');
			colorbar;
			frame = getframe(gcf); 
			im = frame2im(frame); 
			[imind,cm] = rgb2ind(im,256); 
			
            % Write to the GIF File 
			if i == rangeStart 
				imwrite(imind,cm,'anim_cavFracPlot.gif','gif', 'Loopcount',inf); 
			else 
				imwrite(imind,cm,'anim_cavFracPlot.gif','gif','WriteMode','append'); 
			end
			
		  
    end
    close();
end