% Petrov-Galerkin stabilization terms for Navier-Stokes equations
% Part 1: (tauM/rho)(rho*u.gradN).(rho*dU/dt) + 
%         (tauM/rho)(rho*u.gradN).(rho*u.gradU) +
%         (tauM/rho)(gradN).(rho*dU/dt)
%
function [LHS, RHS] = PetrovGalerkinTerms1(Sol, fluid, vapFrac, dens, visc, pmc, solver, iif, jjf, ...
                                          xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
                                          nElem, nQuad, nen, ndof, elemType, ...
                                          LHS, RHS, Flag3D)
                                      
%% Initialize the quadrature solution space                                      
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);
sA4 = zeros(nen^2*nElem,nQuad);
sA5 = zeros(nen^2*nElem,nQuad);
sA6 = zeros(nen^2*nElem,nQuad);      
sA7 = zeros(nen^2*nElem,nQuad); 
sA8 = zeros(nen^2*nElem,nQuad);
sA9 = zeros(nen^2*nElem,nQuad);
sA10 = zeros(nen^2*nElem,nQuad);
sA11 = zeros(nen^2*nElem,nQuad);
sA12 = zeros(nen^2*nElem,nQuad);      
sA13 = zeros(nen^2*nElem,nQuad); 
sA14 = zeros(nen^2*nElem,nQuad);
sA15 = zeros(nen^2*nElem,nQuad);

%% Map to local element space
for p = 1:nQuad  
	if (Flag3D ~= 1)
	   J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
			yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
		if size(J,2)==1
			J = J';
		end
		volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
		xsInv(:,4) = J(:,1).*J(:,4) - J(:,2).*J(:,3) ;
		xsInv(:,4) = 1./xsInv(:,4) ;
		xsInv(:,1) = xsInv(:,4).* J(:,4) ;
		xsInv(:,2) = xsInv(:,4).* (-J(:,2)) ;
		xsInv(:,3) = xsInv(:,4).* (-J(:,3)) ;
		xsInv(:,4) = xsInv(:,4).* J(:,1) ;
		
		if (strcmp(elemType,'3Tri'))
%                     c1 = 1.154700538379252 ;
%                     c2 = 0.577350269189626 ;
%                     
%                     a1 = c1 * xsInv(:,1) + c2 * xsInv(:,3) ;
%                     a2 = c1 * xsInv(:,3) + c2 * xsInv(:,1) ;
%                     gijDown(:,1) = xsInv(:,1).* a1 + xsInv(:,3).* a2 ;
%                     gijDown(:,3) = xsInv(:,2).* a1 + xsInv(:,4).* a2 ;
%                     
%                     a1 = c1 * xsInv(:,2) + c2 * xsInv(:,4) ;
%                     a2 = c1 * xsInv(:,4) + c2 * xsInv(:,2) ;
%                     gijDown(:,2) = xsInv(:,2).* a1 + xsInv(:,4).* a2 ;
			gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
			gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
			gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
            
		elseif (strcmp(elemType,'4Quad'))
			
            gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
			gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
			gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
		
        end
	else
		J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
			 yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
			 zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
		if size(J,2)==1
			J = J';
		end
		volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
				  J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
				  J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
			  
		xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
		xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
		xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
		xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
		xsInv(:,9) = 1./xsInv(:,9) ;
		xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
		xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
		xsInv(:,3) = xsInv(:,9).*xsInv(:,3);

		xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
		xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
		xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));

		xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
		xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
		xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));

		if (strcmp(elemType,'4Tet'))
			c1 = 1.259921049894873E0 ;
			c2 = 6.299605249474365D-01;

			a1 = c1 * xsInv(:,1) + c2 * (xsInv(:,4) + xsInv(:,7)) ;
			a2 = c1 * xsInv(:,4) + c2 * (xsInv(:,1) + xsInv(:,7)) ;
			a3 = c1 * xsInv(:,7) + c2 * (xsInv(:,1) + xsInv(:,4)) ;
			gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;

			a1 = c1 * xsInv(:,2) + c2 * (xsInv(:,5) + xsInv(:,8)) ;
			a2 = c1 * xsInv(:,5) + c2 * (xsInv(:,2) + xsInv(:,8)) ;
			a3 = c1 * xsInv(:,8) + c2 * (xsInv(:,2) + xsInv(:,5)) ;
			gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
			gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;

			a1 = c1 * xsInv(:,3) + c2 * (xsInv(:,6) + xsInv(:,9)) ;
			a2 = c1 * xsInv(:,6) + c2 * (xsInv(:,3) + xsInv(:,9)) ;
			a3 = c1 * xsInv(:,9) + c2 * (xsInv(:,3) + xsInv(:,6)) ;          
			gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* a3;
			gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* a3;
			gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* a3;
            
		elseif (strcmp(elemType,'6Prism'))
			c1 = 1.154700538379252E0 ;
			c2 = 5.773502691896259E-01;

			a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
			a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);

			gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;

			a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
			a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
			gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
			gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);

			a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
			a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
			gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
			gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
			gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
            
		elseif (strcmp(elemType,'8Hex'))
            
			gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,4).*xsInv(:,4) + xsInv(:,7).*xsInv(:,7);
			gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,5).*xsInv(:,5) + xsInv(:,8).*xsInv(:,8);
			gijDown(:,3) = xsInv(:,3).*xsInv(:,3) + xsInv(:,6).*xsInv(:,6) + xsInv(:,9).*xsInv(:,9);
			gijDown(:,4) = xsInv(:,1).*xsInv(:,2) + xsInv(:,4).*xsInv(:,5) + xsInv(:,7).*xsInv(:,8);
			gijDown(:,5) = xsInv(:,2).*xsInv(:,3) + xsInv(:,5).*xsInv(:,6) + xsInv(:,8).*xsInv(:,9);
			gijDown(:,6) = xsInv(:,1).*xsInv(:,3) + xsInv(:,4).*xsInv(:,6) + xsInv(:,7).*xsInv(:,9);
            
		end
	end
         
    volume = abs(volume);
    
	if (Flag3D ~= 1)
		DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
		DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);

		locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
		locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
		
		VG1 = gijDown(:,1) .* locUX + gijDown(:,3) .* locUY ;
		VG2	= gijDown(:,3) .* locUX + gijDown(:,2) .* locUY ;

		VGV	= VG1 .* locUX + VG2 .* locUY ;

		GG = gijDown(:,1).^2 + gijDown(:,2).^2 + 2 * ( gijDown(:,3).^2 );
	else
		DNDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*Nx(:,p)'+ ...
			  (J(:,6).*J(:,7)-J(:,4).*J(:,9))*Ny(:,p)'+ ...
			  (J(:,4).*J(:,8)-J(:,5).*J(:,7))*Nz(:,p)')./repmat(volume,1,nen);
		DNDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*Nx(:,p)'+ ...
			  (J(:,1).*J(:,9)-J(:,7).*J(:,3))*Ny(:,p)'+ ...
			  (J(:,2).*J(:,7)-J(:,1).*J(:,8))*Nz(:,p)')./repmat(volume,1,nen);
		DNDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*Nx(:,p)'+ ...
			  (J(:,3).*J(:,4)-J(:,1).*J(:,6))*Ny(:,p)'+ ...
			  (J(:,1).*J(:,5)-J(:,4).*J(:,2))*Nz(:,p)')./repmat(volume,1,nen);

		locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
		locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);
		locUZ  = sum(repmat(N(p,:),nElem,1).*uz,2);
		
		VG1 = gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
		VG2	= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
		VG3	= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

		VGV	= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

		GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
			2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
	end
	
	localvapFrac  = sum(repmat(N(p,:),nElem,1).*vapFrac,2); % edit - new
	
	dens = sum(repmat(N(p,:),nElem,1).*dens,2);
	visc = sum(repmat(N(p,:),nElem,1).*visc,2);
            
    tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
    tauM = tauM.^(0.5);
    tauM = 1./tauM;
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
            % 1st term
            Aij_1 = gW(p)*(locUX.*DNDx(:,i)*N(p,j));
            Aij_2 = gW(p)*(locUY.*DNDy(:,i)*N(p,j));
            Aij_1 = Aij_1.*volume.*tauM.*dens;
            Aij_2 = Aij_2.*volume.*tauM.*dens;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
                    
            % 2nd term
            Aij_4 = gW(p)*locUX.*locUX.*(DNDx(:,i).*DNDx(:,j));
            Aij_5 = gW(p)*locUX.*locUY.*(DNDx(:,i).*DNDy(:,j));
            Aij_7 = gW(p)*locUY.*locUX.*(DNDy(:,i).*DNDx(:,j));
            Aij_8 = gW(p)*locUY.*locUY.*(DNDy(:,i).*DNDy(:,j)); 
            Aij_4 = Aij_4.*volume.*tauM.*dens;
            Aij_5 = Aij_5.*volume.*tauM.*dens;
            Aij_7 = Aij_7.*volume.*tauM.*dens;
            Aij_8 = Aij_8.*volume.*tauM.*dens;
            sA4(index+1:index+nElem,p) = Aij_4;
            sA5(index+1:index+nElem,p) = Aij_5;
            sA7(index+1:index+nElem,p) = Aij_7;
            sA8(index+1:index+nElem,p) = Aij_8;

            % 3rd term
            Aij_13 = gW(p)*(DNDx(:,i).*N(p,j));
            Aij_14 = gW(p)*(DNDy(:,i).*N(p,j));
            Aij_13 = Aij_13.*volume.*tauM;
            Aij_14 = Aij_14.*volume.*tauM; 
            sA13(index+1:index+nElem,p) = Aij_13;
            sA14(index+1:index+nElem,p) = Aij_14;
			
			if (Flag3D == 1)
				Aij_3 = gW(p)*(locUZ.*DNDz(:,i)*N(p,j));
				Aij_6 = gW(p)*locUX.*locUZ.*(DNDx(:,i).*DNDz(:,j));
				Aij_9 = gW(p)*locUY.*locUZ.*(DNDy(:,i).*DNDz(:,j));
				Aij_10 = gW(p)*locUZ.*locUX.*(DNDz(:,i).*DNDx(:,j));
				Aij_11 = gW(p)*locUZ.*locUY.*(DNDz(:,i).*DNDy(:,j));
				Aij_12 = gW(p)*locUZ.*locUZ.*(DNDz(:,i).*DNDz(:,j));
                Aij_15 = gW(p)*(DNDz(:,i).*N(p,j));               
                Aij_3 = Aij_3.*volume.*tauM.*dens;			
				Aij_6 = Aij_6.*volume.*tauM.*dens;
				Aij_9 = Aij_9.*volume.*tauM.*dens;
				Aij_10 = Aij_10.*volume.*tauM.*dens;
				Aij_11 = Aij_11.*volume.*tauM.*dens;
				Aij_12 = Aij_12.*volume.*tauM.*dens;
                Aij_15 = Aij_15.*volume.*tauM;
                sA3(index+1:index+nElem,p) = Aij_3;
				sA6(index+1:index+nElem,p) = Aij_6;
				sA9(index+1:index+nElem,p) = Aij_9;
				sA10(index+1:index+nElem,p) = Aij_10;
				sA11(index+1:index+nElem,p) = Aij_11;
				sA12(index+1:index+nElem,p) = Aij_12;			
				sA15(index+1:index+nElem,p) = Aij_15;
			end			      
            index = index + nElem;
        end
    end
end

%% Summation of all quadrature data
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);       
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA7 = sum(sA7,2);
sA8 = sum(sA8,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);

if (Flag3D == 1)
	sA3 = sum(sA3,2);
	sA6 = sum(sA6,2);
	sA9 = sum(sA9,2);
	sA10 = sum(sA10,2);
	sA11 = sum(sA11,2);
	sA12 = sum(sA12,2);
	sA15 = sum(sA15,2);
end      

%% Assemble the equation matricies element spaces     
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);
A8 = sparse(iif,jjf,sA8,ndof,ndof);
A13 = sparse(iif,jjf,sA13,ndof,ndof);
A14 = sparse(iif,jjf,sA14,ndof,ndof);

if (Flag3D == 1)
    A3 = sparse(iif,jjf,sA3,ndof,ndof);
    A6 = sparse(iif,jjf,sA6,ndof,ndof);
    A9 = sparse(iif,jjf,sA9,ndof,ndof);
    A10 = sparse(iif,jjf,sA10,ndof,ndof);
    A11 = sparse(iif,jjf,sA11,ndof,ndof);
    A12 = sparse(iif,jjf,sA12,ndof,ndof);
    A15 = sparse(iif,jjf,sA15,ndof,ndof);
end
        
ZeroF = sparse(ndof,ndof);

%% Construct the algebraic equation matrices
if (Flag3D == 1)
	stabK1 = [A1+A2+A3 ZeroF ZeroF ZeroF;...
			  ZeroF A1+A2+A3 ZeroF ZeroF;...
			  ZeroF ZeroF A1+A2+A3 ZeroF;...
			  ZeroF ZeroF ZeroF ZeroF];
          
	stabK2 = [A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF ZeroF ZeroF;...
			  ZeroF A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF ZeroF;...
			  ZeroF ZeroF A4+A5+A6+A7+A8+A9+A10+A11+A12 ZeroF;...
			  ZeroF ZeroF ZeroF ZeroF];
          
	stabG1 = [ZeroF ZeroF ZeroF ZeroF;...
			  ZeroF ZeroF ZeroF ZeroF;...
			  ZeroF ZeroF ZeroF ZeroF;...
			  A13 A14 A15 ZeroF];
else
	stabK1 = [A1+A2 ZeroF ZeroF;...
			  ZeroF A1+A2 ZeroF;...
			  ZeroF ZeroF ZeroF];
          
	stabK2 = [A4+A5+A7+A8 ZeroF ZeroF;...
			  ZeroF A4+A5+A7+A8 ZeroF;...
			  ZeroF ZeroF ZeroF];
          
	stabG1 = [ZeroF ZeroF ZeroF;...
			  ZeroF ZeroF ZeroF;...
			  A13 A14 ZeroF];
end

%% Assemble left and right handed contributions
stabK1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabK1;
stabG1Flow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabG1;

% Left-hand side matrix
LHS = LHS + (stabK1Flow + stabK2 + stabG1Flow);
        
clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15
clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear Aij_13 Aij_14 Aij_15
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 xsInv

uAlpha = Sol.uAlpha(:,:,1) ;
uAlpha = [uAlpha(:); zeros(ndof,1)];
uDotAlpha = Sol.uDotAlpha(:,:,1) ;
uDotAlpha = [uDotAlpha(:); zeros(ndof,1)];
if (size(Sol.p,1) == 1 & size(Sol.p,2)>1)
   Sol.p = [Sol.p]';
end

if (Flag3D == 1)
	p1 = [zeros(ndof*3,1);Sol.p];
	gravVec = [ones(3*ndof,1); zeros(ndof,1)];
else
	p1 = [zeros(ndof*2,1);Sol.p];
	gravVec = [ones(2*ndof,1); zeros(ndof,1)];
end
        
% Right-hand side vector
RHS = RHS -(stabK1 + stabG1) * uDotAlpha(:) - (stabK2) * uAlpha(:) ;                                   
                                      
end