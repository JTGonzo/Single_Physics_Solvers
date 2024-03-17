function [Sol, structNormIndicator] = hyperelasticMaterial(solver, solid, pmc, Sol, cnnStructure, coordStructure, ...
                                               elemType, ndofS, nenStructure, nElemStructure, BCStructure, tracSrc)

% Quadrature rules for elements
if strcmp(elemType,'3Tri')
    gP = ...
   [1/3,  1/3
    4/3,  1/3
    1/3,  4/3] ;
    gW = ...
   [2/3,  2/3,  2/3] ;
 
    N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
    N(:,2) = 0.5.*(gP(:,1)) ;
    N(:,3) = 0.5.*(gP(:,2)) ;
    
    Nx(:,1) = -0.5.*ones(3,1) ;
    Nx(:,2) =  0.5.*ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -0.5.*ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  0.5.*ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
end

Nx = Nx' ;
Ny = Ny' ;
nQuad = length(gW) ;
 
iis = zeros(nenStructure^2*nElemStructure,1); 
jjs = zeros(nenStructure^2*nElemStructure,1);
index = 0;
for i = 1:nenStructure
   for j = 1:nenStructure
      iis(index+1:index+nElemStructure) = double(cnnStructure(:,i)); 
      jjs(index+1:index+nElemStructure) = double(cnnStructure(:,j));  
      index = index + nElemStructure;
   end
end
        
% Satisfy boundary conditions
Sol.uS(solid.DirichletU,1,1) = solid.DirichletUval ;
Sol.uS(solid.DirichletV,2,1) = solid.DirichletVval ;
 
% Interpolate for alpha values for Gen-alpha
Sol.uSAlpha = Sol.uSPrev + pmc.alphaS.*(Sol.uS - Sol.uSPrev) ;
Sol.uSDotAlpha = Sol.uSDotPrev + pmc.alphaSM.*(Sol.uSDot - Sol.uSDotPrev) ;
Sol.dispSAlpha = Sol.dispSPrev + pmc.alphaS.*(Sol.dispS - Sol.dispSPrev) ;
        
% Navier-Stokes equations
xxs = zeros(size(cnnStructure));
yys = zeros(size(cnnStructure));
dispSx = zeros(size(cnnStructure));
dispSy = zeros(size(cnnStructure));

for i=1:nenStructure
   xxs(:,i) = coordStructure(cnnStructure(:,i),2);
   yys(:,i) = coordStructure(cnnStructure(:,i),3);
   dispSx(:,i) =  Sol.dispSAlpha(cnnStructure(:,i),1,1) ;
   dispSy(:,i) =  Sol.dispSAlpha(cnnStructure(:,i),2,1) ;
end
        
% Form element matrix and assemble Galerkin terms
[LHS, RHS] = GalerkinTermsStructural(Sol, solid, pmc, solver, iis, jjs, xxs, yys, ...
                                    gW, N, Nx, Ny, nElemStructure, nQuad, nenStructure, ndofS, dispSx,...
                                    dispSy, tracSrc);
                                         
% Solve the linear system
        
% Select the unknown nodal values
freeNodesU = unique(solid.DirichletU);
freeNodesU = setdiff(1:size(coordStructure,1),[freeNodesU]);
freeNodesV = unique(solid.DirichletV);
freeNodesV = setdiff(1:size(coordStructure,1),[freeNodesV]);

freeNodes = [freeNodesU';freeNodesV' + size(coordStructure,1)];
        
result = Sol.uSAlpha(:,:,1);
result = result(:);
resultDot = Sol.uSDotAlpha(:,:,1);
resultDot = resultDot(:);
resultDispS = Sol.dispSAlpha(:,:,1);
resultDispS = resultDispS(:);

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments
result(freeNodes) = result(freeNodes) + Increment;
resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaSM/(pmc.gammaS*pmc.alphaS*solver.dt))*Increment;

Sol.uSAlpha(:,:,1) = reshape(result(1:2*ndofS),[],2);
Sol.uSDotAlpha(:,:,1) = reshape(resultDot(1:2*ndofS),[],2);

% Update the solution
Sol.uS = Sol.uSPrev + (1/pmc.alphaS)*( Sol.uSAlpha - Sol.uSPrev ) ;
Sol.uSDot = Sol.uSDotPrev + (1/pmc.alphaSM)*( Sol.uSDotAlpha - Sol.uSDotPrev ) ;
Sol.dispS = Sol.dispSPrev + (solver.dt*(Sol.uSPrev + solver.dt*((1.0/2.0 - pmc.betaS)*Sol.uSDotPrev + pmc.betaS*Sol.uSDot))) ;
        
structNormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
fprintf('structural: %e, ', structNormIndicator);
clear freeNodesU freeNodesV freeNodes
clear result resultDot
        
end

