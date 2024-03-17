%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                2D Linear Elastic Material Equation Solver               %
%                                                                         % 
%              \nabla . \sigma^m = 0,    in \Omega^s,                     %                      
%              \sigma = grad \eta + grad \eta^T + (div \eta)I,            %
%                                                                         %
%    where \sigma is the stress experienced by the ALE mesh due to        %
%    strain induced by structural movement, \eta is the mesh disp-        %
%    lacement for the fluid nodes and I is the identity tensor.           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol, SVKnormIndicator] = SVKMaterialGlobal(Sol, solver, solid, tracSrc, pmc, cnn, crd, elemType,...
                            nen, ndof, nElem, BCStructure, BCBarFix)

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

velS = Sol.uS;
velSPrev = Sol.uSPrev;
velSDot = Sol.uSDot;
velSDotPrev = Sol.uSDotPrev;
dispS = Sol.dispS;
dispSPrev = Sol.dispSPrev;

% Form the local to global map
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

% Satisfy boundary conditions
bcFix = unique(BCBarFix(:));
velS(bcFix,1,1) = solid.DirichletUval ;
velS(bcFix,2,1) = solid.DirichletVval ;

% Interpolate for alpha values for Gen-alpha
velSAlpha = velSPrev + pmc.alpha.*(velS - velSPrev) ;
velSDotAlpha = velSDotPrev + pmc.alphaM.*(velSDot - velSDotPrev) ;
dispSAlpha = dispSPrev + solver.dt.*0.5.*velSAlpha ;
        
% SVK Nonlinear hyperelastic equation
xxs = zeros(size(cnn));
yys = zeros(size(cnn));
ux = zeros(size(cnn));
uy = zeros(size(cnn));
dispx = zeros(size(cnn));
dispy = zeros(size(cnn));
tracLocalx = zeros(size(cnn));
tracLocaly = zeros(size(cnn));

% Localize the data to each element
for i=1:nen
   xxs(:,i) = crd(cnn(:,i),2);
   yys(:,i) = crd(cnn(:,i),3);
   ux(:,i) =  velSAlpha(cnn(:,i),1,1) ;
   uxDot(:,i) = velSDotAlpha(cnn(:,i),1,1) ;
   uy(:,i) =  velSAlpha(cnn(:,i),2,1) ;
   uyDot(:,i) = velSDotAlpha(cnn(:,i),2,1) ;
   dispx(:,i) = dispSAlpha(cnn(:,i),1);
   dispy(:,i) = dispSAlpha(cnn(:,i),2);
   dispxPrev(:,i) = dispSPrev(cnn(:,i),1);
   dispyPrev(:,i) = dispSPrev(cnn(:,i),2);
end

% Form element matrix and assemble Galerkin terms
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xxs*[Nx(:,p)], xxs*[Ny(:,p)],...
         yys*[Nx(:,p)], yys*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
    locUX  = sum(repmat(N(p,:),nElem,1).*(ux),2);
    locUY  = sum(repmat(N(p,:),nElem,1).*(uy),2);
    locUXDot  = sum(repmat(N(p,:),nElem,1).*(uxDot),2);
    locUYDot  = sum(repmat(N(p,:),nElem,1).*(uyDot),2);
    locdispX  = sum(repmat(N(p,:),nElem,1).*(dispx),2);
    locdispY  = sum(repmat(N(p,:),nElem,1).*(dispy),2);
	locTracX  = sum(repmat(N(p,:),nElem,1).*(tracLocalx),2);
    locTracY  = sum(repmat(N(p,:),nElem,1).*(tracLocaly),2);
    
    locgradXdispX = sum(DNDx.*dispx,2);
    locgradYdispX = sum(DNDy.*dispx,2);
    locgradXdispY = sum(DNDx.*dispy,2);
    locgradYdispY = sum(DNDy.*dispy,2);
    
    locgradXdispXPrev = sum(DNDx.*dispxPrev,2);
    locgradYdispXPrev = sum(DNDy.*dispxPrev,2);
    locgradXdispYPrev = sum(DNDx.*dispyPrev,2);
    locgradYdispYPrev = sum(DNDy.*dispyPrev,2);
    
    
    traceE = 0.5.*( 2.*locgradXdispXPrev + 2.*locgradYdispYPrev +...
                    locgradXdispXPrev.*locgradXdispXPrev + locgradXdispYPrev.*locgradXdispYPrev + ...
                    locgradYdispXPrev.*locgradYdispXPrev + locgradYdispYPrev.*locgradYdispYPrev ) ; 
    traceENew = 0.5.*( 2.*locgradXdispX + 2.*locgradYdispY +...
                    locgradXdispX.*locgradXdispX + locgradXdispY.*locgradXdispY + ...
                    locgradYdispX.*locgradYdispX + locgradYdispY.*locgradYdispY ) ; 
        
    L = solid.lambdaS ;
    Mu = solid.muS ;
    dt = solver.dt ;
    
    F11 = 1 + locgradXdispXPrev ;
    F12 = locgradYdispXPrev ;
    F21 = locgradXdispYPrev ;
    F22 = 1 + locgradYdispYPrev ;
    
    F11New = 1 + locgradXdispX ;
    F12New = locgradYdispX ;
    F21New = locgradXdispY ;
    F22New = 1 + locgradYdispY ;
    
    E11 = 0.5.*(locgradXdispXPrev + locgradXdispXPrev + locgradXdispXPrev.*locgradXdispXPrev + locgradXdispYPrev.*locgradXdispYPrev) ;
    E12 = 0.5.*(locgradYdispXPrev + locgradXdispYPrev + locgradXdispXPrev.*locgradYdispXPrev + locgradXdispYPrev.*locgradYdispYPrev) ;
    E21 = 0.5.*(locgradXdispYPrev + locgradYdispXPrev + locgradYdispXPrev.*locgradXdispXPrev + locgradYdispYPrev.*locgradXdispYPrev) ;
    E22 = 0.5.*(locgradYdispYPrev + locgradYdispYPrev + locgradYdispXPrev.*locgradYdispXPrev + locgradYdispYPrev.*locgradYdispYPrev) ;
    
    
    E11New = 0.5.*(locgradXdispX + locgradXdispX + locgradXdispX.*locgradXdispX + locgradXdispY.*locgradXdispY) ;
    E12New = 0.5.*(locgradYdispX + locgradXdispY + locgradXdispX.*locgradYdispX + locgradXdispY.*locgradYdispY) ;
    E21New = 0.5.*(locgradXdispY + locgradYdispX + locgradYdispX.*locgradXdispX + locgradYdispY.*locgradXdispY) ;
    E22New = 0.5.*(locgradYdispY + locgradYdispY + locgradYdispX.*locgradYdispX + locgradYdispY.*locgradYdispY) ;
    
    sigma11 = L.*traceENew.*(1+locgradXdispX) + 2.*Mu.*((1+locgradXdispX).*E11New + locgradYdispX.*E21New);
    sigma12 = L.*traceENew.*(locgradYdispX) + 2.*Mu.*((1+locgradXdispX).*E12New + locgradYdispX.*E22New);
    sigma21 = L.*traceENew.*(locgradXdispY) + 2.*Mu.*((locgradXdispY).*E11New + (1+locgradYdispY).*E21New);
    sigma22 = L.*traceENew.*(1+locgradYdispY) + 2.*Mu.*((locgradXdispY).*E12New + (1+locgradYdispY).*E22New);
                
    index = 0;
    for i = 1:nen
        for j = 1:nen    
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*solid.dens ;
            Mij = Mij.*volume;
            sA(index+1:index+nElem,p) = Mij;
            
            % Galerkin terms for K11 (Psi1 with U1)
            C1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            C2 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            C3 = gW(p)*(DNDy(:,i).*DNDx(:,j));
            C4 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            
            Aij_1 = ( (L.*(E11.*dt/2 + E22.*dt/2 + F11.^2*dt/2) + ...
                              2*Mu.*(F11.^2*dt/2 + F12.^2*dt/4 + E11*dt/2)).*C1 + ...
                             (L.*(F11.*F12.*dt/2) + ...
                              2*Mu.*(F11.*F12*dt/4 + E21.*dt/2)).*C2 + ...
                             (L.*(F11.*F21*dt/2) + ...
                              2*Mu.*(F11.*F21*dt/2 + F12.*F22*dt/4)).*C3 + ...
                             (L.*(F12.*F21*dt/2) + ...
                              2*Mu.*(F11.*F22*dt/4)).*C4 ).*volume;
            Aij_2 = ( (L.*(F11.*F21*dt/2) + ...
                              2*Mu.*(F21.*F11*dt/2 + F22.*F12*dt/4)).*C1 + ...
                             (L.*(F11.*F22*dt/2) + ...
                              2*Mu.*(F21.*F12*dt/4)).*C2 + ...
                             (L.*(E11.*dt/2 + E22*dt/2 + F21.^2*dt/2) + ...
                              2*Mu.*(F21.^2*dt/2 + F22.^2*dt/4 + E11*dt/2)).*C3 + ...
                             (L.*(F21.*F22*dt/2) + ...
                              2*Mu.*(F21.*F22*dt/4 + E21*dt/2)).*C4 ).*volume;
            Aij_3 =( (L.*(F11.*F12*dt/2) + ...
                              2*Mu.*(F11.*F12*dt/4 + E12*dt/2)).*C1 + ...
                             (L.*(E11*dt/2 + E22*dt/2 + F12.^2*dt/2) + ...
                              2*Mu.*(F11.^2*dt/4 + F12.^2*dt/2 + E22*dt/2)).*C2 + ...
                             (L.*(F11.*F22*dt/2) + ...
                              2*Mu.*(F12.*F21*dt/4)).*C3 + ...
                             (L.*(F12.*F22*dt/2) + ...
                              2*Mu.*(F11.*F21*dt/4 + F12.*F22*dt/2)).*C4 ).*volume;
            Aij_4 = ( (L.*(F12.*F21*dt/2) + ...
                              2*Mu.*(F22.*F11*dt/4)).*C1 + ...
                             (L.*(F12.*F22*dt/2) + ...
                              2*Mu.*(F21.*F11*dt/4 + F22.*F12*dt/2)).*C2 + ...
                             (L.*(F21.*F22*dt/2) + ...
                              2*Mu.*(F21.*F22*dt/4 + E12.*dt/2)).*C3 + ...
                             (L.*(E11*dt/2 + E22*dt/2 + F22.^2*dt/2) + ...
                              2*Mu.*(F21.^2*dt/4 + F22.^2*dt/2 + E22*dt/2)).*C4 ).*volume;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
            sA3(index+1:index+nElem,p) = Aij_3;
            sA4(index+1:index+nElem,p) = Aij_4;
            
            % Galerkin source term (gravity force)
            Aij_5 = gW(p)*(N(p,i)*N(p,j));
            Aij_5 = Aij_5.*solid.dens.*volume ;
            sA5(index+1:index+nElem,p) = Aij_5;			
            
            % Galerkin source term (Stress term)
            Aij_6 = gW(p)*(DNDx(:,i).*sigma11 + DNDy(:,i).*sigma21)*N(p,j);
            Aij_7 = gW(p)*(DNDx(:,i).*sigma12 + DNDy(:,i).*sigma22)*N(p,j);
            Aij_6 = Aij_6.*volume ;
            Aij_7 = Aij_7.*volume ;
            sA6(index+1:index+nElem,p) = Aij_6;
            sA7(index+1:index+nElem,p) = Aij_7;

            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA = sum(sA,2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
sA7 = sum(sA7,2);

% Assemble the matrix
Ms = sparse(iif,jjf,sA,ndof,ndof);  
ZeroF = sparse(ndof,ndof);
        
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A6 = sparse(iif,jjf,sA6,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);

Ms = [Ms ZeroF;...
      ZeroF Ms];
Ks = [A1 A2;
       A3 A4];
Src = [A5.*solid.gravFrc(1) ZeroF;...  
       ZeroF A5.*solid.gravFrc(2)];
Src1 = [A6]*[ones(1*ndof,1)];
Src2 = [A7]*[ones(1*ndof,1)];
   
MStruct = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Ms;

% Left-hand side matrix
LHS = MStruct + Ks;
        
% clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
% clear sA10 sA11 sA12 sA13 sA14 sA15 sA16
% clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16
% clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16

uSAlpha = velSAlpha(:,:,1) ;
uSAlpha = [uSAlpha(:)];
uSDotAlpha = velSDotAlpha(:,:,1) ;
uSDotAlpha = [uSDotAlpha(:)];

gravVec = [ones(2*ndof,1)];

% Right-hand side vector
RHS = -(Ms * uSDotAlpha(:)) - [Src1; Src2];
RHS = RHS + (Src)*gravVec ;
RHS = RHS + tracSrc;

% Solve the linear system
        
% Select the unknown nodal values
freeNodesU = [bcFix] ;
freeNodesU = setdiff(unique(BCStructure(:)),[freeNodesU]);
freeNodesV = [bcFix];
freeNodesV = setdiff(unique(BCStructure(:)),[freeNodesV]);

freeNodes = [freeNodesU;freeNodesV + size(crd,1)];
        
result = velSAlpha(:,:,1);
result = result(:);
resultDot = velSDotAlpha(:,:,1);
resultDot = resultDot(:);

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments
result(freeNodes) = result(freeNodes) + Increment;
resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Increment ;

velSAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
velSDotAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);

% Update the solution
velS = velSPrev + (1/pmc.alpha)*( velSAlpha - velSPrev ) ;
velSDot = velSDotPrev + (1/pmc.alphaM)*( velSDotAlpha - velSDotPrev ) ;
dispS = dispSPrev + solver.dt.*velSPrev + solver.dt.*0.5.*(velS - velSPrev) ;
        
SVKnormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
fprintf('SVK: %e, ', SVKnormIndicator);
clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot

% Map local structure to global data
Sol.uS = velS;
Sol.uSPrev = velSPrev;
Sol.uSDot = velSDot;
Sol.uSDotPrev = velSDotPrev;
Sol.dispS = dispS;
Sol.dispSPrev = dispSPrev;

end