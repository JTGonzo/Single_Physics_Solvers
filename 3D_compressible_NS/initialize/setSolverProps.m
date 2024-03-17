%% setSolverProps
%  This function defines th solver, output, threshold
% increment details
%
function [solver] = setSolverProps()

% Nonlinear iteration data:
	solver.nLIterMin = 1 ;
	solver.nLIterMax = 7 ;

% Time parameters    
	solver.dt = 5e-9 ;
	solver.maxSteps = 100 ;

% Naview Stokes parameters      
	solver.rhoinfty = 1.0 ;
	solver.nLTol_NS = 5e-4 ;
    
% Cavitation parameters    
	solver.nLTol_cav = 1e-3 ;
	solver.epsilon = 1.0 ; %0.01 ;
	solver.h = 1.5e-5;  %0.01;= 
	
% Output data generation details:
	solver.outFreq = 1 ;
	solver.outOthdFreq = 1;
end