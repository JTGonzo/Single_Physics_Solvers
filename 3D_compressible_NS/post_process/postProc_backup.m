function [] = postProc(Sol, solver, bc, timeStep, fileOthdId, ...
						wrkDir, problemString, crd, elemType, Flag3D)

    figHandle = figure;
	% ax = gca;
	% ax.XLim = [0 0.0001];	
	% ax.YLim = [0 0.0002];
	% ax.ZLim = [0.0 1.0];	
% % % % %     subplot(1,5,1); hold off;
% % % % %     trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.u(:,1)', ...
% % % % %         'FaceColor', 'interp', 'EdgeColor', 'interp');
% % % % %     view(0,90), axis equal, axis on;
% % % % % %     title('Sloshing tank problem', 'FontSize', 14)
% % % % %     set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
% % % % %     colormap('jet');
% % % % %      colorbar;
% % % % %     
% % % % %     subplot(1,5,2); hold off;
% % % % %     trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.u(:,2)', ...
% % % % %         'FaceColor', 'interp', 'EdgeColor', 'interp');
% % % % %     view(0,90), axis equal, axis on;
% % % % % %     title('Sloshing tank problem', 'FontSize', 14)
% % % % %     set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
% % % % %     colormap('jet');
% % % % %      colorbar;
% % % % %     
% % % % %     subplot(1,5,3); hold off;
% % % % %     trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.p, ...
% % % % %         'FaceColor', 'interp', 'EdgeColor', 'interp');
% % % % %     view(0,90), axis equal, axis on;
% % % % % %     title('Sloshing tank problem', 'FontSize', 14)
% % % % %     set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
% % % % %     colormap('jet');
% % % % %      colorbar; 
% % % % % 	 
% % % % %     subplot(1,5,4); hold off;
    trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.vapFrac', ...
        'FaceColor', 'interp', 'EdgeColor', 'interp');

    xlim([-0.005 0.005]);
    ylim([-0.005 0.005]);
    zlim([0 1]);
    view(0,90), axis equal, axis on;
    %title('Sloshing tank problem', 'FontSize', 14)
    set(gcf, 'PaperUnits','Centimeters','PaperPosition',[0 0 30 20]);
    colormap('jet');
    colorbar;
    
% % % % %     subplot(1,5,5); hold off; 
% % % % %     trisurf(Sol.elem,Sol.node(:,1),Sol.node(:,2),zeros(size(Sol.node,1),1));
% % % % %     view(0,90), axis equal, axis on;
    
     pause(0.1);

     title(sprintf('Cavitation Fraction at time %9.8f',timeStep*solver.dt), 'FontSize', 14)
      % Capture the plot as an image 
      frame = getframe(gca); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if timeStep == 1 
          imwrite(imind,cm,'bubble.gif','gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,'bubble.gif','gif','WriteMode','append'); 
      end
    
    if (mod(timeStep,solver.outFreq)==0)
    saveas(gcf,['./plot_' num2str(timeStep) '.fig']) ;
    end
    
    close;
    
    % Output the *.plt files
    if (mod(timeStep,solver.outFreq)==0)
    crd1 = [Sol.node];
    uu = Sol.u ;
    pp = Sol.p ;
    vapFracvapFrac = Sol.vapFrac ;
    elemelem = Sol.elem ;
    
%     indxFree = find(Sol.type==0) ;
%     crd1(indxFree,:) = [];
%     uu(indxFree,:,:) = [];
%     pp(indxFree,:) = [];
%     vapFracvapFrac(indxFree,:) = [];
%     indx2 = find(crd1(:,1)~=0) ;
%     map = [crd1(:,1) indx2];
%     elemelem = elemelem(:) ;
%     [~,ii] = ismember(elemelem,map(:,1));
%     elemelem = reshape(ii,[],3);
%     crd1(:,1) = [] ;
        
    data = [crd1 uu pp vapFracvapFrac ];
    clear uu pp vapFracvapFrac ;
    filename = sprintf('%s/%s.%d.plt',wrkDir,problemString,timeStep);
    fileId = fopen(filename,'w');

    FileTitle = 'velocity plot';

    fprintf(fileId,' TITLE = \"%s\"\n',FileTitle);

    if (strcmp(FileTitle,'velocity plot') & Flag3D == 1)
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"p\", \"vapFrac\"\n');
    elseif (strcmp(FileTitle,'velocity plot') & Flag3D ~= 1)
        fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"p\", \"<greek>f</greek>\" \n');
    end
    
    if (strcmp(elemType,'3Tri'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',timeStep*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'4Quad'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',timeStep*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'4Tet'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON\n',timeStep*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'6Prism'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',timeStep*solver.dt,size(crd1,1),size(elemelem,1));
    elseif (strcmp(elemType,'8Hex'))
        fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',timeStep*solver.dt,size(crd1,1),size(elemelem,1));
    end
    
    if (Flag3D == 1)
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f \n',data');
    else
        fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f \n',data');
    end

    if (strcmp(elemType,'3Tri'))
        fprintf(fileId,'%d %d %d\n',elemelem');
    elseif (strcmp(elemType,'4Quad'))
        fprintf(fileId,'%d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'4Tet'))
        fprintf(fileId,'%d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'6Prism'))
        fprintf(fileId,'%d %d %d %d %d %d\n',elemelem');
    elseif (strcmp(elemType,'8Hex'))
        fprintf(fileId,'%d %d %d %d %d %d %d %d\n',elemelem');
    end
    clear crd1 elemelem ;
    fclose(fileId);
    end
    
    % Output the *.othd files
    if (mod(timeStep,solver.outOthdFreq)==0)
    data = [Sol.node(bc.left.nodes',2) Sol.vapFrac(bc.left.nodes',1)];
    fprintf(fileOthdId,'%d\n',size(bc.left.nodes',1));
    fprintf(fileOthdId,'%f %12.10f\n',data');
    end

end