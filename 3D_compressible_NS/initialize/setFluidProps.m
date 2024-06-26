%% setFluidProps
%  This function defines the fluid medium material properties, initial
%  and boundary conditions and transition criteria
%
function [fluid] = setFluidProps(Flag3D)

% Material properties of multiphased fluid 
	% fluid.Re = 100 ;
	fluid.dens1_0 = 1000 ;
	fluid.dens2_0 = 0.01389 ;
	fluid.visc1 = 1e-3 ;
	fluid.visc2 = 1e-5;
	fluid.pVap = 2320; %% edit - new - vapor pressure of water - check again
	fluid.surfTens = 0.0 ; %0.07198 ;
	fluid.alphasf = 3*sqrt(2)/4; % 1.0;
	fluid.gamma = 1e-2; % 1.0;
    
% Initial boundary conditions:
	fluid.gravFrc = [0.0, 0.0, 0.0];
	fluid.vel0 = [0.0, 0.0, 0.0];
	fluid.p0 = 1e5 ;
	fluid.pInf = 1e5 ;
	fluid.n0 = 1e+8;
    
% Phase change thresholds  
	fluid.dNuc = 1e-4;
	fluid.gammaV = 1.2;
	fluid.cC = 1e+4;
	fluid.cV = 1e+4;
    
	if(Flag3D == 1)
		fluid.alphaNuc = (pi*fluid.n0*fluid.dNuc^3/6)/(1 + pi*fluid.n0*fluid.dNuc^3/6);
	else % mod for 2D bubble
		fluid.alphaNuc = (pi*fluid.n0*fluid.dNuc^2/4)/(1 + pi*fluid.n0*fluid.dNuc^2/4);
    end	
	fluid.vapFrac0 = fluid.alphaNuc;
end