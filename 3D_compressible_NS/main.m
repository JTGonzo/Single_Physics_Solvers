%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%               CML Multiphase Fluid FEM Solver                 %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
tic

%% Load simulation data and initialize solver parameters
outDir = 'post_process/data/' ;
addpath('initialize/')
addpath('post_process/')
addpath('solver_cavitation/')
addpath('solver_fluid/')
    
saveRestartID = 1;
restartID = 0;

if(restartID)
	load('saveCaseData_60.mat');
	restartTimeStep = timeStep;
    solver.maxSteps = 100;
else
	newCase();
	load('newCase.mat');
	delete('newCase.mat');
	restartTimeStep = 0;
end

fileOthd = sprintf('%s/%s.othd',outDir,problemString);
fileOthdId = fopen(fileOthd,'w');

%% Time loop starts here
for timeStep = restartTimeStep+1:solver.maxSteps

    fprintf('Time step:%d\n',timeStep);
    fprintf(fileOthdId,'Time %f\n',timeStep*solver.dt);
	
    % Predict the solution
    Sol.u = Sol.uPrev ;
    Sol.uDot = (pmc.gamma - 1)/pmc.gamma * Sol.uDotPrev ;
    Sol.p = Sol.pPrev ;
    Sol.vapFrac = Sol.vapFracPrev ;
    Sol.vapFracDot = (pmc.gamma - 1)/pmc.gamma * Sol.vapFracDotPrev ;
    
    [fluid] = updateDensVisc(Sol, fluid);
	
    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax
	
		% Satisfy boundary conditions
		[Sol] = satisfyBC(Sol, fluid, crd, bc, Flag3D);
		
		% Update Mixture density and viscosity
		[fluid] = updateDensVisc(Sol, fluid);
		
		% Interpolate for alpha values for Gen-alpha
		[Sol] = interpGenAlpha(Sol, pmc);
		
        % Navier Stokes
        [Sol, NSnormIndicator] = navierStokes(solver, fluid, pmc, Sol, cnn, crd, ...
                                               elemType, ndof, nen, nElem, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, timeStep, Flag3D, bc); 
		
        % Cavitation - Schnerr Sauer                                   
		[Sol, cavNormIndicator] = cavitation(Sol, fluid, bc, solver, BCTop, BCBottom,...
								  BCLeft, BCRight, pmc, cnn, crd, elemType, nen,...
								  ndof, nElem, Flag3D);									  
        
        % Check convergence criteria
        if ((NSnormIndicator < solver.nLTol_NS) && (cavNormIndicator < solver.nLTol_cav))
            break;
        end
    end
    fprintf('\n');
    
%% Copy current variables to previous variables
    Sol.uPrev = Sol.u ;
    Sol.uDotPrev = Sol.uDot ;
    Sol.pPrev = Sol.p ;
    Sol.vapFracPrev = Sol.vapFrac ;
    Sol.vapFracDotPrev = Sol.vapFracDot ;

%% Output and vizualize incremental time-step details   
	postProc(Sol, solver, bc, timeStep, fileOthdId, ...
						outDir, problemString, crd, elemType, Flag3D);

    timeArr(timeStep) = timeStep*solver.dt;
    massArr(timeStep) = Sol.massV;
    
    if (mod(timeStep,solver.outFreq)==0)
        saveFile = strcat('case_data_new/saveCaseData_',num2str(timeStep),'.mat');
        save(saveFile); 
    end
toc    
end
fclose(fileOthdId);
fprintf('\n');

%% Save restart data for future problem run
if(saveRestartID)
    restartFile = strcat('case_data_new/restartCaseData_',num2str(timeStep),'.mat');
	save(restartFile);
end