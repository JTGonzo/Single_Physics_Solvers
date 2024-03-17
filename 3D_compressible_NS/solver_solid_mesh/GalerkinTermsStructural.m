% Galerkin terms for Structural Solver

function [LHS, RHS] = GalerkinTermsStructural(Sol, solid, pmc, solver, iis, jjs, xxs, yys, ...
                                    gW, N, Nx, Ny, nElemStructure, nQuad, nenStructure, ndofS, dispSx,...
                                    dispSy, tracSrc)

sA = [] ;
sA1 = zeros(nenStructure^2*nElemStructure,nQuad); 
sA2 = zeros(nenStructure^2*nElemStructure,nQuad);
sA3 = zeros(nenStructure^2*nElemStructure,nQuad);
sA4 = zeros(nenStructure^2*nElemStructure,nQuad);
sA5 = zeros(nenStructure^2*nElemStructure,nQuad);
sA6 = zeros(nenStructure^2*nElemStructure,nQuad);      
sA7 = zeros(nenStructure^2*nElemStructure,nQuad); 
sA8 = zeros(nenStructure^2*nElemStructure,nQuad);
sA9 = zeros(nenStructure^2*nElemStructure,nQuad);
sA10 = zeros(nenStructure^2*nElemStructure,nQuad);
sA11 = zeros(nenStructure^2*nElemStructure,nQuad);
sA12 = zeros(nenStructure^2*nElemStructure,nQuad);       
sA13 = zeros(nenStructure^2*nElemStructure,nQuad);        
sA14 = zeros(nenStructure^2*nElemStructure,nQuad);
sA15 = zeros(nenStructure^2*nElemStructure,nQuad);
sA16 = zeros(nenStructure^2*nElemStructure,nQuad);
sA17 = zeros(nenStructure^2*nElemStructure,nQuad);
sA18 = zeros(nenStructure^2*nElemStructure,nQuad);
sA19 = zeros(nenStructure^2*nElemStructure,nQuad);
sA20 = zeros(nenStructure^2*nElemStructure,nQuad);
sA21 = zeros(nenStructure^2*nElemStructure,nQuad);
sA22 = zeros(nenStructure^2*nElemStructure,nQuad);       
sA23 = zeros(nenStructure^2*nElemStructure,nQuad);        
sA24 = zeros(nenStructure^2*nElemStructure,nQuad);

for p = 1:nQuad  
    J = [xxs*[Nx(:,p)], xxs*[Ny(:,p)],...
         yys*[Nx(:,p)], yys*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nenStructure);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nenStructure);
      
    locgradEta1X = sum(dispSx.*DNDx,2) ;
    locgradEta1Y = sum(dispSx.*DNDy,2) ; 
	locgradEta2X = sum(dispSy.*DNDx,2) ;
    locgradEta2Y = sum(dispSy.*DNDy,2) ; 

	locTraceE = 0.5*(locgradEta1X.*locgradEta1X + locgradEta2X.*locgradEta2X + 2*(locgradEta1X + locgradEta2Y) + locgradEta1Y.*locgradEta1Y + locgradEta2Y.*locgradEta2Y);
            
    index = 0;
    for i = 1:nenStructure
        for j = 1:nenStructure
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*solid.dens ;
            Mij = Mij.*volume;
            sA(index+1:index+nElemStructure,p) = Mij;
			
            % Term 1 \mu^s * \grad \psi : \grad \eta^s
            Aij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_2 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            Aij_3 = gW(p)*(DNDy(:,i).*DNDx(:,j));
			Aij_4 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_1 = Aij_1.*volume*solid.muS;
            Aij_2 = Aij_2.*volume*solid.muS;
            Aij_3 = Aij_3.*volume*solid.muS;
            Aij_4 = Aij_4.*volume*solid.muS;			
            sA1(index+1:index+nElemStructure,p) = Aij_1;
            sA2(index+1:index+nElemStructure,p) = Aij_2;	
            sA3(index+1:index+nElemStructure,p) = Aij_3;
            sA4(index+1:index+nElemStructure,p) = Aij_4;	

            % Term 2 \mu^s * \grad \psi : (\grad \eta^s)^T
            Aij_5 = gW(p)*(DNDx(:,i).*DNDx(:,j));
			Aij_6 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_5 = Aij_5.*volume*solid.muS;
            Aij_6 = Aij_6.*volume*solid.muS;
            sA5(index+1:index+nElemStructure,p) = Aij_5;
            sA6(index+1:index+nElemStructure,p) = Aij_6;	

            % Term 3 \lambda^s * (\grad \psi : trace(E)I)	
			Aij_7 = gW(p)*DNDx(:,i).*(locgradEta1X.*DNDx(:,j) + 2*DNDx(:,j) + locgradEta1Y.*DNDy(:,j));
			Aij_8 = gW(p)*DNDx(:,i).*(locgradEta2X.*DNDx(:,j) + 2*DNDy(:,j) + locgradEta2Y.*DNDy(:,j));
			Aij_9 = gW(p)*DNDy(:,i).*(locgradEta1X.*DNDx(:,j) + 2*DNDx(:,j) + locgradEta1Y.*DNDy(:,j));
			Aij_10 = gW(p)*DNDy(:,i).*(locgradEta2X.*DNDx(:,j) + 2*DNDy(:,j) + locgradEta2Y.*DNDy(:,j));			
            Aij_7 = 0.5*solid.lambdaS*Aij_7.*volume;
            Aij_8 = 0.5*solid.lambdaS*Aij_8.*volume;
            Aij_9 = 0.5*solid.lambdaS*Aij_9.*volume;
            Aij_10 = 0.5*solid.lambdaS*Aij_10.*volume;
            sA7(index+1:index+nElemStructure,p) = Aij_7;
            sA8(index+1:index+nElemStructure,p) = Aij_8;	
            sA9(index+1:index+nElemStructure,p) = Aij_9;
            sA10(index+1:index+nElemStructure,p) = Aij_10;

            % Term 4 \lambda^s * trace(E) * (\grad \psi : \grad \eta^s)	
            Aij_11 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_12 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            Aij_13 = gW(p)*(DNDy(:,i).*DNDx(:,j));
			Aij_14 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_11 = solid.muS*Aij_11.*volume.*locTraceE;
            Aij_12 = solid.muS*Aij_12.*volume.*locTraceE;
            Aij_13 = solid.muS*Aij_13.*volume.*locTraceE;
            Aij_14 = solid.muS*Aij_14.*volume.*locTraceE;			
            sA11(index+1:index+nElemStructure,p) = Aij_11;
            sA12(index+1:index+nElemStructure,p) = Aij_12;	
            sA13(index+1:index+nElemStructure,p) = Aij_13;
            sA14(index+1:index+nElemStructure,p) = Aij_14;	

            % Term 5 \mu^s * (\grad \psi : (\grad \eta^s \grad \eta^s))	
            Aij_15 = gW(p)*(locgradEta1X.*DNDx(:,i).*DNDx(:,j) + locgradEta2X.*DNDy(:,i).*DNDx(:,j));
            Aij_16 = gW(p)*(locgradEta1Y.*DNDx(:,i).*DNDx(:,j) + locgradEta2Y.*DNDy(:,i).*DNDx(:,j));
            Aij_17 = gW(p)*(locgradEta1X.*DNDx(:,i).*DNDy(:,j) + locgradEta2X.*DNDy(:,i).*DNDy(:,j));
			Aij_18 = gW(p)*(locgradEta1Y.*DNDx(:,i).*DNDy(:,j) + locgradEta2Y.*DNDy(:,i).*DNDy(:,j));
            Aij_15 = solid.muS*Aij_15.*volume;
            Aij_16 = solid.muS*Aij_16.*volume;
            Aij_17 = solid.muS*Aij_17.*volume;
            Aij_18 = solid.muS*Aij_18.*volume;			
            sA15(index+1:index+nElemStructure,p) = Aij_15;
            sA16(index+1:index+nElemStructure,p) = Aij_16;	
            sA17(index+1:index+nElemStructure,p) = Aij_17;
            sA18(index+1:index+nElemStructure,p) = Aij_18;		

            % Term 6 \mu^s * (\grad \psi : (\grad \eta^s (\grad \eta^s)^T))	
            Aij_19 = gW(p)*(locgradEta1X.*DNDx(:,i).*DNDx(:,j) + locgradEta1Y.*DNDy(:,i).*DNDx(:,j) + locgradEta2X.*DNDx(:,i).*DNDy(:,j) + locgradEta2Y.*DNDy(:,i).*DNDy(:,j));
            Aij_19 = solid.muS*Aij_19.*volume;
			sA19(index+1:index+nElemStructure,p) = Aij_19;	

            % Term 7 \mu^s * (\grad \psi : (\grad \eta^s \grad \eta^s (\grad \eta^s)^T))	
            Aij_20 = gW(p)*((locgradEta1X.*locgradEta1X + locgradEta1Y.*locgradEta2X).*DNDx(:,i).*DNDx(:,j));
			Aij_21 = gW(p)*((locgradEta1X.*locgradEta1Y + locgradEta1Y.*locgradEta2Y).*DNDx(:,i).*DNDy(:,j));
			Aij_22 = gW(p)*((locgradEta2Y.*locgradEta2Y + locgradEta1Y.*locgradEta2X).*DNDy(:,i).*DNDy(:,j));
			Aij_23 = gW(p)*((locgradEta1X.*locgradEta2X + locgradEta2X.*locgradEta2Y).*DNDy(:,i).*DNDx(:,j));
            Aij_20 = solid.muS*Aij_20.*volume;
			Aij_21 = solid.muS*Aij_21.*volume;
			Aij_22 = solid.muS*Aij_22.*volume;
			Aij_23 = solid.muS*Aij_23.*volume;
			sA20(index+1:index+nElemStructure,p) = Aij_20;
			sA21(index+1:index+nElemStructure,p) = Aij_21;
			sA22(index+1:index+nElemStructure,p) = Aij_22;
			sA23(index+1:index+nElemStructure,p) = Aij_23;			
                    
            % Galerkin source term (gravity force)
            Aij_24 = gW(p)*(N(p,i)*N(p,j));
            Aij_24 = Aij_24.*solid.dens.*volume ;
            sA24(index+1:index+nElemStructure,p) = Aij_24;
                              
            index = index + nElemStructure;
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
sA8 = sum(sA8,2);
sA9 = sum(sA9,2);
sA10 = sum(sA10,2);
sA11 = sum(sA11,2);
sA12 = sum(sA12,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);
sA15 = sum(sA15,2);
sA16 = sum(sA16,2);
sA17 = sum(sA17,2);
sA18 = sum(sA18,2);
sA19 = sum(sA19,2);
sA20 = sum(sA20,2);
sA21 = sum(sA21,2);
sA22 = sum(sA22,2);
sA23 = sum(sA23,2);
sA24 = sum(sA24,2);
 
% Assemble the matrix
Ms = sparse(iis,jjs,sA,ndofS,ndofS);  
ZeroF = sparse(ndofS,ndofS);
        
A1 = sparse(iis,jjs,sA1,ndofS,ndofS);
A2 = sparse(iis,jjs,sA2,ndofS,ndofS);
A3 = sparse(iis,jjs,sA3,ndofS,ndofS);
A4 = sparse(iis,jjs,sA4,ndofS,ndofS);
A5 = sparse(iis,jjs,sA5,ndofS,ndofS);
A6 = sparse(iis,jjs,sA6,ndofS,ndofS);
A7 = sparse(iis,jjs,sA7,ndofS,ndofS);
A8 = sparse(iis,jjs,sA8,ndofS,ndofS);
A9 = sparse(iis,jjs,sA9,ndofS,ndofS);
A10 = sparse(iis,jjs,sA10,ndofS,ndofS);
A11 = sparse(iis,jjs,sA11,ndofS,ndofS);
A12 = sparse(iis,jjs,sA12,ndofS,ndofS);
A13 = sparse(iis,jjs,sA13,ndofS,ndofS);
A14 = sparse(iis,jjs,sA14,ndofS,ndofS);
A15 = sparse(iis,jjs,sA15,ndofS,ndofS);
A16 = sparse(iis,jjs,sA16,ndofS,ndofS);
A17 = sparse(iis,jjs,sA17,ndofS,ndofS);
A18 = sparse(iis,jjs,sA18,ndofS,ndofS);
A19 = sparse(iis,jjs,sA19,ndofS,ndofS);
A20 = sparse(iis,jjs,sA20,ndofS,ndofS);
A21 = sparse(iis,jjs,sA21,ndofS,ndofS);
A22 = sparse(iis,jjs,sA22,ndofS,ndofS);
A23 = sparse(iis,jjs,sA23,ndofS,ndofS);
A24 = sparse(iis,jjs,sA24,ndofS,ndofS);

Ms = [Ms ZeroF ZeroF;...
      ZeroF Ms ZeroF;...
      ZeroF ZeroF ZeroF];
	  
term1 = [A1 A2 ZeroF;...
		 A3 A4 ZeroF;...
		 ZeroF ZeroF ZeroF];
		 
term2 = [A5+A6 ZeroF ZeroF;...
		 ZeroF A5+A6 ZeroF;...
		 ZeroF ZeroF ZeroF];

term3 = [A7 A8 ZeroF;...
		 A9 A10 ZeroF;...
		 ZeroF ZeroF ZeroF];

term4 = [A11 A12 ZeroF;...
		 A13 A14 ZeroF;...
		 ZeroF ZeroF ZeroF];

term5 = [A15 A16 ZeroF;...
		 A17 A18 ZeroF;...
		 ZeroF ZeroF ZeroF];	

term6 = [A19 ZeroF ZeroF;...
		 ZeroF A19 ZeroF;...
		 ZeroF ZeroF ZeroF];	

term7 = [A20+A21+A22+A23 ZeroF ZeroF;...
		 ZeroF A20+A21+A22+A23 ZeroF;...
		 ZeroF ZeroF ZeroF];		 
	  
Src = [A24.*solid.gravFrc(1) ZeroF ZeroF;...  
       ZeroF A24.*solid.gravFrc(2) ZeroF;...;...
       ZeroF ZeroF ZeroF];
   
Mstruct = (pmc.alphaSM/(pmc.gammaS*pmc.alphaS*solver.dt))*Ms;
term1S = (pmc.betaS*solver.dt/pmc.gammaS)*term1;
term2S = (pmc.betaS*solver.dt/pmc.gammaS)*term2;
term3S = (pmc.betaS*solver.dt/pmc.gammaS)*term3;
term4S = (pmc.betaS*solver.dt/pmc.gammaS)*term4;
term5S = (pmc.betaS*solver.dt/pmc.gammaS)*term5;
term6S = (pmc.betaS*solver.dt/pmc.gammaS)*term6;
term7S = (pmc.betaS*solver.dt/pmc.gammaS)*term7;

% Left-hand side matrix
LHS = Mstruct - term1S - term2S - term3S - term4S - term5S - term6S - term7S; 
        
clear Mij Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16 Aij_17 Aij_18 Aij_19 Aij_20 Aij_21 Aij_22 Aij_23 Aij_24
clear sA sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15 sA16 sA17 sA18 sA19 sA20 sA21 sA22 sA23 sA24
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18 A19 A20 A21 A22 A23 A24

uSAlpha = Sol.uSAlpha(:,:,1) ;
uSAlpha = [uSAlpha(:); zeros(ndofS,1)];
uSDotAlpha = Sol.uSDotAlpha(:,:,1) ;
uSDotAlpha = [uSDotAlpha(:); zeros(ndofS,1)];

gravVec = [ones(2*ndofS,1); zeros(ndofS,1)];

% Right-hand side vector
RHS = -(Ms * uSDotAlpha(:)) + (term1 + term2 + term3 + term4 + term5 + term6 + term7)* uSAlpha(:) ;
RHS = RHS + (Src)*gravVec + tracSrc;

end