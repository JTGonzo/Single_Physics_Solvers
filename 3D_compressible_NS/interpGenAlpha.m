%% interpGenAlpha
%  This function estimates and intermediate solution for the domain
%  variable per the generalized alpha integration scheme
%
function [Sol] = interpGenAlpha(Sol, pmc)

	Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;
	Sol.uDotAlpha = Sol.uDotPrev + pmc.alphaM.*(Sol.uDot - Sol.uDotPrev) ;
	Sol.vapFracAlpha = Sol.vapFracPrev + pmc.alpha.*(Sol.vapFrac - Sol.vapFracPrev) ;
	Sol.vapFracDotAlpha = Sol.vapFracDotPrev + pmc.alphaM.*(Sol.vapFracDot - Sol.vapFracDotPrev) ;
			
end