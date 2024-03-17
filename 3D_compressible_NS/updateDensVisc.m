%% updateDensVisc
%  This function updates the fluid domains material properties (denisty,
%  vsicosity) based on the current setp's solver outputs
%
function [fluid] = updateDensVisc(Sol, fluid)

% 	% Tait's equation for liquid (sea water)
% 	fluid.dens1 = fluid.dens1_0./(1-0.315*log((2996.0 + Sol.p/100000)/2997.0));

% 	% Isentropic update for vapor
% 	fluid.dens2 = fluid.dens2_0*(Sol.p/fluid.pVap).^(1/fluid.gammaV); 
% 	fluid.dens2(fluid.dens2 < fluid.dens2_0) = fluid.dens2_0;
	
	fluid.dens1 = fluid.dens1_0*ones(size(Sol.p));
	fluid.dens2 = fluid.dens2_0*ones(size(Sol.p));	

	densLiquid = fluid.dens1_0;
	densVapor  = fluid.dens2_0;
    
	viscLiquid = fluid.visc1;
	viscVapor  = fluid.visc2;	
	
	fluid.dens = densLiquid.*Sol.vapFrac + densVapor.*(1-Sol.vapFrac) ;
	fluid.visc = viscLiquid.*Sol.vapFrac + viscVapor.*(1-Sol.vapFrac) ;
end