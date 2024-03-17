%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%        2D Incompressible Navier-Stokes Finite-element Solver            %
%                                                                         % 
%   \rho u_t + \rho u.div u - div(mu*grad u) + grad p = f in \Omega,      %
%                         div u = 0  in \Omega,                           %
%                                                                         %
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,         %
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol, NSnormIndicator] = navierStokes(solver, fluid, pmc, Sol, cnn, crd, ...
                                               elemType, ndof, nen, nElem, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, timeStep, Flag3D, bc)                                          
                                           
% Quadrature rules for elements
[gP, gW, N, Nx, Ny, Nz, nQuad] = defineQuadratureData(elemType, Flag3D);
 
iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); 
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end
        
%% Navier-Stokes equations
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
zzf = zeros(size(cnn));
ux = zeros(size(cnn));
uy = zeros(size(cnn));

if (Flag3D == 1)
	uz = zeros(size(cnn));
end

% Localize the data to each element
for i=1:nen
	xxf(:,i) = crd(cnn(:,i),1);
	yyf(:,i) = crd(cnn(:,i),2);
    ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
	uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
    uxDot(:,i) = Sol.uDotAlpha(cnn(:,i),1,1) ;
	uyDot(:,i) = Sol.uDotAlpha(cnn(:,i),2,1) ;
    
	if (Flag3D == 1)
		zzf(:,i) = crd(cnn(:,i),3);
        uz(:,i) =  Sol.uAlpha(cnn(:,i),3,1) ;
        uzDot(:,i) = Sol.uDotAlpha(cnn(:,i),3,1) ;
    end
       
	pres(:,i) = Sol.p(cnn(:,i),1) ;
	pressure(:,i) = Sol.p(cnn(:,i),1) ;
    
    vapFrac(:,i) = Sol.vapFracAlpha(cnn(:,i),1) ;
	vapFracPrev(:,i) = Sol.vapFracPrev(cnn(:,i),1) ;
	vapFracDot(:,i) = Sol.vapFracDotAlpha(cnn(:,i),1) ;
    
    dens(:,i) = fluid.dens(cnn(:,i),1) ;
	dens1(:,i) = fluid.dens1(cnn(:,i),1) ;
	dens2(:,i) = fluid.dens2(cnn(:,i),1) ;
	
	visc(:,i) = fluid.visc(cnn(:,i),1) ;	
end
        
%% Form element matrix and assemble Galerkin & Petrov-Galerkin terms
[LHS, RHS] = GalerkinTerms(Sol, fluid, bc, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, vapFrac, dens, visc, pressure, dens1, dens2, pmc, solver, iif, jjf, xxf, yyf, zzf,...
                            ux, uy, uz, gW, N, Nx, Ny, Nz, nElem, nQuad, nen, ndof, crd, cnn, timeStep, Flag3D);
                                
% [LHS, RHS] = PetrovGalerkinTerms1(Sol, fluid, vapFrac, dens, visc, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D);
% [LHS, RHS] = PetrovGalerkinTerms2(Sol, fluid, vapFrac, dens, visc, pressure, dens1, dens2, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D);
% [LHS, RHS] = PetrovGalerkinTerms3(Sol, fluid, vapFrac, dens, visc, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D);
% [LHS, RHS] = discontinuityCapTerms(Sol, fluid, vapFrac, dens, visc, pres, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D, cnn, crd);								  
                                      
%% Select the unknown nodal values
zerotypenodes = find(Sol.type == 0);
freeNodes1 = setdiff(1:size(crd,1),[zerotypenodes]);  %unique([bc.left.nodes';bc.right.nodes' ]); 
freeNodes2 = setdiff(1:size(crd,1),[zerotypenodes]);  %unique([bc.left.nodes';bc.right.nodes' ]);
freeNodes3 = setdiff(1:size(crd,1),[zerotypenodes]);
freeNodes4 = setdiff(1:size(crd,1),[zerotypenodes]);

if (Flag3D == 1)
	freeNodes = [freeNodes1';freeNodes2' + size(crd,1); freeNodes3' + 2*size(crd,1); freeNodes4' + 3*size(crd,1)];
    % freeNodes = [freeNodes1';freeNodes2' + size(crd,1); freeNodes3' + 3*size(crd,1)];
else
	freeNodes = [freeNodes1';freeNodes2' + size(crd,1); freeNodes3' + 2*size(crd,1)];
end
        
result = Sol.uAlpha(:,:,1);
result = result(:);
result = [result; Sol.p];

resultDot = Sol.uDotAlpha(:,:,1);
resultDot = resultDot(:);
resultDot = [resultDot; Sol.p];

%% Solve the linear system
% Increment = bicgstab(LHS(freeNodes,freeNodes),RHS(freeNodes),1e-6,2000);
% Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
Increment = LHS\RHS;
% Increment = linsolve(LHS(freeNodes,freeNodes), RHS(freeNodes));
        
%% Update the fluid variables increment
% result(freeNodes) = result(freeNodes) + Increment;
% resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Increment ;
result = result + Increment;
resultDot = resultDot + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Increment ;

if (Flag3D == 1)
	Sol.uAlpha(:,:,1) = reshape(result(1:3*ndof),[],3);
	Sol.uDotAlpha(:,:,1) = reshape(resultDot(1:3*ndof),[],3);
	Sol.p = result((3*ndof+1):(4*ndof));
else
	Sol.uAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
	Sol.uDotAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);
	Sol.p = result((2*ndof+1):(3*ndof));    
end

%Sol.p(Sol.p < (fluid.pVap - 50)) = (fluid.pVap - 50); %tempSuraj

%% Update the step's output solution
Sol.u = Sol.uPrev + (1/pmc.alpha)*( Sol.uAlpha - Sol.uPrev ) ;
Sol.uDot = Sol.uDotPrev + (1/pmc.alphaM)*( Sol.uDotAlpha - Sol.uDotPrev ) ;

%% Check the solution conversion and conditioning
% NSnormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
NSnormIndicator =  norm(Increment)/norm(result) ;

fprintf('NS: %e, ', NSnormIndicator);
clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot
        
end

