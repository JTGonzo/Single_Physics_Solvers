%% newCase
%  This function defines the initial conditions of the fluid domain
%  and identifies the type of BC on the domain faces
%
function [bc, Sol] = defineBC(Sol, BCLeft, BCRight, BCTop, BCBottom, BCFront, BCBack, Flag3D)

%% Boundary conditions:
	% Dirichlet: 1
	% Neumann: 2
	% Symmetry/slip: 3
	% No boundary (2D): 0
	% variable type for velocity: 1 for ux, 2 for uy, 3 for uz

%% Type of initial condition applied to each face    
    bc.left.type = 2;
	bc.right.type = 2 ;
	bc.top.type = 2 ;
	bc.bottom.type = 2 ;
	bc.front.type = 2 ;
	bc.back.type = 2 ;

% Define left nodes    
	bc.left.nNodes = size(unique(BCLeft),1);
	bc.left.nodes = unique(BCLeft);
	bc.left.nodes(bc.left.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.left.value = [0.0 0.0 0.0] ;
        bc.left.var = [1 2 3] ;
	else
        bc.left.value = [0.0 0.0] ;
        bc.left.var = [1 2] ;
    end

% Define right nodes     
	bc.right.nNodes = size(unique(BCRight),1);
	bc.right.nodes = unique(BCRight);
	bc.right.nodes(bc.right.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.right.value = [0.0 0.0 0.0] ;
        bc.right.var = [1 2 3] ;
	else
        bc.right.value = [0.0 0.0] ;
        bc.right.var = [1 2] ;
    end

% Define top nodes     
	bc.top.nNodes = size(unique(BCTop),1);
	bc.top.nodes = unique(BCTop);
	bc.top.nodes(bc.top.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.top.value = [0.0 0.0 0.0] ;
        bc.top.var = [1 2 3] ;
	else
        bc.top.value = [0.0 0.0] ;
        bc.top.var = [1 2] ;
    end

% Define bottom nodes       
	bc.bottom.nNodes = size(unique(BCBottom),1);
	bc.bottom.nodes = unique(BCBottom);
	bc.bottom.nodes(bc.bottom.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.bottom.value = [0.0 0.0 0.0] ;
        bc.bottom.var = [1 2 3] ;
	else
        bc.bottom.value = [0.0 0.0] ;
        bc.bottom.var = [1 2] ;
    end

% Define front nodes      
	bc.front.nNodes = size(unique(BCFront),1);
	bc.front.nodes = unique(BCFront);
	bc.front.nodes(bc.front.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.front.value = [0.0 0.0 0.0] ;
        bc.front.var = [1 2 3] ;
	else
        bc.front.value = [0.0 0.0] ;
        bc.front.var = [1 2] ;
    end	

% Define back nodes       
	bc.back.nNodes = size(unique(BCBack),1);
	bc.back.nodes = unique(BCBack);
	bc.back.nodes(bc.back.nodes==0) = [];
    
	if (Flag3D == 1)
        bc.back.value = [0.0 0.0 0.0] ;
        bc.back.var = [1 2 3] ;
	else
        bc.back.value = [0.0 0.0] ;
        bc.back.var = [1 2] ;
    end		

%% Define Dirichlet and Neumann boundaries
	dirichlet = [] ;
	neumann = [BCLeft; BCRight; BCTop; BCBottom; BCFront; BCBack];
	neumann = unique(sort(neumann,2),'rows');
	
	Sol.Dirichlet = dirichlet ;
	Sol.Neumann = neumann ;

end