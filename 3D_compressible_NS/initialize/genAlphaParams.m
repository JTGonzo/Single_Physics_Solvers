%% genAlphaParams
%  This function initializes the general aloha integration parameters
%
function [pmc] = genAlphaParams(solver)

	pmc.alphaM = 0.5*(3-solver.rhoinfty)/(1+solver.rhoinfty) ;
	pmc.alpha = 1/(1+solver.rhoinfty) ;
	pmc.gamma = 0.5 + pmc.alphaM - pmc.alpha ;

end