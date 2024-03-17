%% newCase
%  This function reads the geometry *.mat files and sets all the problem
%  /simulation parameters (material, solver settings, output files, etc.)
%
function [] = newCase()

    meshDir = 'mesh/' ; 
    addpath('mesh/')
    addpath('initialize/')
    
	problemString = 'RB_sq_v2_temp' ;
	elemType = '8Hex' ;
	problemType = '3D' ;
    nEmax = 10000 ;
    
	if (strcmp(problemType,'3D'))
		Flag3D = 1;
	else
		Flag3D = 0;
	end
	
%% Get crd, cnn, nbc file data
	crdstr = strcat(meshDir,'coord.mat') ;
	cnnstr = strcat(meshDir,'conn.mat') ;	
	BC1str = strcat(meshDir,'BCLeft.mat') ;
	BC2str = strcat(meshDir,'BCRight.mat') ;
	BC3str = strcat(meshDir,'BCTop.mat') ;
	BC4str = strcat(meshDir,'BCBottom.mat') ;    
	if(~Flag3D)
		BC5str = strcat(meshDir,'BCAll.mat') ;
	else
		BC5str = strcat(meshDir,'BCFront.mat') ;	
		BC6str = strcat(meshDir,'BCBack.mat') ;	
		BC7str = strcat(meshDir,'BCAll.mat') ;	
    end
    
    load(crdstr);
    load(cnnstr);
	load(BC1str);
	load(BC2str);
	load(BC3str);
	load(BC4str);
	load(BC5str);   
	if(Flag3D)
		load(BC6str);	
		load(BC7str);	
    end
    
	ndof = size(coord,1) ;
    nodeId = coord(:,1) ;
    
	crd = coord(:,2:4) ;
	cnn = conn(:,1:end) ;
	nen = size(cnn,2) ;
	nElem = size(cnn,1) ;
	ndof = size(crd,1);
    
%% Define solver specific details and fluid properties
	[solver] = setSolverProps();

	[fluid] = setFluidProps(Flag3D);
    
    [pmc] = genAlphaParams(solver); % Gen-alpha parameters

%% Define the quadrature data for the solver
	[gP, gW, N, Nx, Ny, nQuad] = defineQuadratureData(elemType, Flag3D);

%% Initialize Solution Space Variables and define Boundary Conditions
	[Sol] = initializeVar(fluid, ndof, Flag3D, crd);
	
	[bc, Sol] = defineBC(Sol, BCLeft, BCRight, BCTop, BCBottom, BCFront, BCBack, Flag3D);	

%%  Save the initial solution space and define case initial .mat file
	type = uint8(ones(size(crd,1),1));
	
	if(~Flag3D)
		crd(:,3) = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end

	Sol.node = crd ;
	Sol.elem = cnn ;
	Sol.type = type ;

	if(~Flag3D)
		crd = [Sol.node zeros(size(Sol.node,1),1)] ;
	end
	
	save('newCase.mat');
end