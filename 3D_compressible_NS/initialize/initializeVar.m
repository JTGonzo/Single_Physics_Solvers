%% initializeVar
% This function initializes the mutli-domain/phase solution 
% space 
function [Sol] = initializeVar(fluid, ndof, Flag3D, crd)

%% Initialize Fluid Velocity, Acc and Pressure solution space
	if (Flag3D == 1)
		Sol.u = zeros(ndof,3,1);
		Sol.uDot = zeros(ndof,3,1);
		Sol.u(:,1,:) = fluid.vel0(1) ;
		Sol.u(:,2,:) = fluid.vel0(2) ;
		Sol.u(:,3,:) = fluid.vel0(3) ;
        
		Sol.uAlpha = zeros(ndof,3,1) ;
		Sol.uDotAlpha = zeros(ndof,3,1) ;
        
		Sol.uPrev = Sol.u ;
		Sol.uDotPrev = Sol.uDot ;
	else
		Sol.u = zeros(ndof,2,1);
		Sol.uDot = zeros(ndof,2,1);
		Sol.u(:,1,:) = fluid.vel0(1) ;
		Sol.u(:,2,:) = fluid.vel0(2) ;
        
		Sol.uAlpha = zeros(ndof,2,1) ;
		Sol.uDotAlpha = zeros(ndof,2,1) ;
        
		Sol.uPrev = Sol.u ;
		Sol.uDotPrev = Sol.uDot ;
    end

 %% Initialize Cavitation solver and vapor fraction variable space
	Sol.p = ones(ndof,1);
	Sol.pPrev = Sol.p ;
    
    Sol.vapFrac = ones(ndof,1);
	[Sol] = setVapFracIC(Sol, fluid, crd, Flag3D); % vapFrac Initialization
	Sol.vapFracDot = zeros(ndof,1) ; % edit - new
	
	Sol.vapFracAlpha = zeros(ndof,1) ;
	Sol.vapFracDotAlpha = zeros(ndof,1) ;
	
	Sol.vapFracPrev = Sol.vapFrac ;
	Sol.vapFracDotPrev = Sol.vapFracDot ;
end