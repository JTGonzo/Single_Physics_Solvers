function [] = postProc(Sol, solver, bc, timeStep, fileOthdId, ...
						wrkDir, problemString, crd, elemType, Flag3D)
		
	%% Output .fig files
	if (mod(timeStep,solver.outFreq)==0)
	
		% x-velocity	
			trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.u(:, 1)', ...
				'FaceColor', 'interp', 'EdgeColor', 'interp');
			saveas(gcf,['post_process/post_figures/xVelPlot_' num2str(timeStep) '.fig']) ;
			close();	
	
		% y-velocity	
			trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.u(:, 2)', ...
				'FaceColor', 'interp', 'EdgeColor', 'interp');
			saveas(gcf,['post_process/post_figures/yVelPlot_' num2str(timeStep) '.fig']) ;
			close();

		% pressure	
			trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.p, ...
				'FaceColor', 'interp', 'EdgeColor', 'interp');
			saveas(gcf,['post_process/post_figures/pressurePlot_' num2str(timeStep) '.fig']) ;
			close();			
	
%        % density	
% 			trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Fluid.dens, ...
% 				'FaceColor', 'interp', 'EdgeColor', 'interp');
% 			saveas(gcf,['post_process/post_figures/densityPlot_' num2str(timeStep) '.fig']) ;
% 			close();	
            
		% cavitation fraction	
			trisurf(Sol.elem, Sol.node(:,1), Sol.node(:,2), Sol.vapFrac', ...
				'FaceColor', 'interp', 'EdgeColor', 'interp');
			saveas(gcf,['post_process/post_figures/cavFracPlot_' num2str(timeStep) '.fig']) ;
			close();
			
		% % mesh	
			% trisurf(Sol.elem,Sol.node(:,1),Sol.node(:,2),zeros(size(Sol.node,1),1));
			% saveas(gcf,['./meshPlot_' num2str(timeStep) '.fig']) ;
			% close();			
			
    end	
    
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