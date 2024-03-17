% Petrov-Galerkin stabilization terms for Allen-Cahn equations
% Part 1: (tauM/rho)(rho*u.gradN).(rho*dU/dt) + 
%         (tauM/rho)(rho*u.gradN).(rho*u.gradU) +
%         (tauM/rho)(gradN).(rho*dU/dt)
%

function [LHS, RHS] = PetrovGalerkinTermsCav(Sol, fluid, vapFrac, pressure, dens1, dens2, dens, pmc, solver, iif, jjf, ...
                                          xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
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
sA16 = zeros(nen^2*nElem,nQuad); 
sA17 = zeros(nen^2*nElem,nQuad); 
sA18 = zeros(nen^2*nElem,nQuad);
sA19 = zeros(nen^2*nElem,nQuad);
sA20 = zeros(nen^2*nElem,nQuad);
sA21 = zeros(nen^2*nElem,nQuad);

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
    
	localvapFrac  = sum(repmat(N(p,:),nElem,1).*vapFrac,2); 
    % localvapFracPrev = sum(repmat(N(p,:),nElem,1).*vapFracPrev,2);   
	localP = sum(repmat(N(p,:),nElem,1).*pressure,2);
	localDens1 = sum(repmat(N(p,:),nElem,1).*dens1,2);
	localDens2 = sum(repmat(N(p,:),nElem,1).*dens2,2);	
	localDens = sum(repmat(N(p,:),nElem,1).*dens,2);
	   
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
    
	timeA = pmc.alpha ;
    
    reac = 0.0;
	% reac = 0.25*( (1/timeA^3)*localvapFrac.^2 - (3/timeA^3 - 4/timeA^2)*(localvapFrac).*localvapFracPrev + ...
					% +(3/timeA^3 - 8/timeA^2 + 6/timeA)*localvapFracPrev.^2 - 2/timeA ) ...
					% - lagPar*0.5*( (1/(3*timeA^2)).*localvapFrac + (1/3)*(-2/timeA^2 + 3/timeA).*localvapFracPrev );
	         
	tauP = solver.h ./ 2 ./ sqrt( locUX.^2 + locUY.^2 );
	% tauP = (2/solver.dt)^2 + (solver.epsilon^4).* GG * 9 + abs(VGV) + reac.^2;
    % tauP = tauP.^(0.5);
    % tauP = 1./tauP;   
    
	index = 0;
	for i = 1:nen
		for j = 1:nen
			% 1st term
			Aij_1 = gW(p)*(locUX.*DNDx(:,i)*N(p,j));
			Aij_2 = gW(p)*(locUY.*DNDy(:,i)*N(p,j));
			Aij_1 = Aij_1.*volume.*tauP ;
			Aij_2 = Aij_2.*volume.*tauP ;
			sA1(index+1:index+nElem,p) = Aij_1;
			sA2(index+1:index+nElem,p) = Aij_2;
			
			% 2nd term
			Aij_4 = gW(p)*locUX.*locUX.*(DNDx(:,i).*DNDx(:,j));
			Aij_5 = gW(p)*locUX.*locUY.*(DNDx(:,i).*DNDy(:,j));
			Aij_7 = gW(p)*locUY.*locUX.*(DNDy(:,i).*DNDx(:,j));
			Aij_8 = gW(p)*locUY.*locUY.*(DNDy(:,i).*DNDy(:,j)); 
			Aij_4 = Aij_4.*volume.*tauP ;
			Aij_5 = Aij_5.*volume.*tauP ;
			Aij_7 = Aij_7.*volume.*tauP ;
			Aij_8 = Aij_8.*volume.*tauP ;
			sA4(index+1:index+nElem,p) = Aij_4;
			sA5(index+1:index+nElem,p) = Aij_5;
			sA7(index+1:index+nElem,p) = Aij_7;
			sA8(index+1:index+nElem,p) = Aij_8;
			
			% 3rd term
			timeA = pmc.alpha ;
			reacJac = 0.0;                      
			Aij_13 = gW(p)*(reacJac).*(locUX.*DNDx(:,i)*N(p,j));
			Aij_14 = gW(p)*(reacJac).*(locUY.*DNDy(:,i)*N(p,j));
			Aij_13 = Aij_13.*volume.*tauP;
			Aij_14 = Aij_14.*volume.*tauP;
			sA13(index+1:index+nElem,p) = Aij_13;
			sA14(index+1:index+nElem,p) = Aij_14;
			
			% Galerkin source term matrices 
			reac = 0.0;                  
			Aij_16 = gW(p)*(reac).*(locUX.*DNDx(:,i)*N(p,j));
			Aij_17 = gW(p)*(reac).*(locUY.*DNDy(:,i)*N(p,j));
			Aij_16 = Aij_16.*volume.*tauP;
			Aij_17 = Aij_17.*volume.*tauP;
			sA16(index+1:index+nElem,p) = Aij_16;
			sA17(index+1:index+nElem,p) = Aij_17;
			
            % Local multi-phase terms
			rhol = localDens1;                                     					
			rhov = localDens2;
			rhom = localDens;
			Pv = fluid.pVap;
			alphaNuc = fluid.alphaNuc;
			n0 = fluid.n0;
			dTime = solver.dt;
			cC = fluid.cC;
			cV = fluid.cV;
			
			% if(Flag3D == 1)
				rB = ((3/(4*pi*n0))*(1 + alphaNuc - localvapFrac)./localvapFrac).^(1/3);
			% else % mod for 2D bubble
				% rB = ((1/(pi*n0))*(1 + alphaNuc - localvapFrac)./localvapFrac).^(1/2);
			% end
			
            % Galerkin continuity source term
			source1 = cC*localvapFrac.*(1-localvapFrac).*(3*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol.*abs(localP-Pv))).*max(localP-Pv, 0);
			source2 = cV*localvapFrac.*(1+alphaNuc-localvapFrac).*(3*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol.*abs(localP-Pv))).*min(localP-Pv, 0);
			source = source1 + source2;	
			% source = 0.0;
			
			Aij_19 = gW(p)*(source).*(locUX.*DNDx(:,i)*N(p,j));
			Aij_20 = gW(p)*(source).*(locUY.*DNDy(:,i)*N(p,j));
			Aij_19 = Aij_19.*volume.*tauP;
			Aij_20 = Aij_20.*volume.*tauP;
			sA19(index+1:index+nElem,p) = Aij_19;
			sA20(index+1:index+nElem,p) = Aij_20;
			
			if (Flag3D == 1)
                Aij_3 = gW(p)*(locUZ.*DNDz(:,i)*N(p,j));
                Aij_6 = gW(p)*locUX.*locUZ.*(DNDx(:,i).*DNDz(:,j));
                Aij_9 = gW(p)*locUY.*locUZ.*(DNDy(:,i).*DNDz(:,j));
                Aij_10 = gW(p)*locUZ.*locUX.*(DNDz(:,i).*DNDx(:,j));
                Aij_11 = gW(p)*locUZ.*locUY.*(DNDz(:,i).*DNDy(:,j));
                Aij_12 = gW(p)*locUZ.*locUZ.*(DNDz(:,i).*DNDz(:,j));
                Aij_15 = gW(p)*(reacJac).*(locUZ.*DNDz(:,i)*N(p,j));
                Aij_18 = gW(p)*(reac).*(locUZ.*DNDz(:,i)*N(p,j));
                Aij_21 = gW(p)*(source).*(locUZ.*DNDz(:,i)*N(p,j));
                Aij_3 = Aij_3.*volume.*tauP ;
                Aij_6 = Aij_6.*volume.*tauP ;
                Aij_9 = Aij_9.*volume.*tauP ;
                Aij_10 = Aij_10.*volume.*tauP ;
                Aij_11 = Aij_11.*volume.*tauP ;
                Aij_12 = Aij_12.*volume.*tauP ;
                Aij_15 = Aij_15.*volume.*tauP;
                Aij_18 = Aij_18.*volume.*tauP;
                Aij_21 = Aij_21.*volume.*tauP;
                sA3(index+1:index+nElem,p) = Aij_3;
                sA6(index+1:index+nElem,p) = Aij_6;
                sA9(index+1:index+nElem,p) = Aij_9;
                sA10(index+1:index+nElem,p) = Aij_10;
                sA11(index+1:index+nElem,p) = Aij_11;
                sA12(index+1:index+nElem,p) = Aij_12;               
                sA15(index+1:index+nElem,p) = Aij_15;
                sA18(index+1:index+nElem,p) = Aij_18;              
                sA21(index+1:index+nElem,p) = Aij_21;
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
sA16 = sum(sA16,2);
sA17 = sum(sA17,2);
sA19 = sum(sA19,2);
sA20 = sum(sA20,2);

if (Flag3D == 1)
	sA3 = sum(sA3,2);
	sA6 = sum(sA6,2);
	sA9 = sum(sA9,2);
	sA10 = sum(sA10,2);
	sA11 = sum(sA11,2);
	sA12 = sum(sA12,2);
	sA15 = sum(sA15,2);
	sA18 = sum(sA18,2);
	sA21 = sum(sA21,2);
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
A16 = sparse(iif,jjf,sA16,ndof,ndof);
A17 = sparse(iif,jjf,sA17,ndof,ndof);
A19 = sparse(iif,jjf,sA19,ndof,ndof);
A20 = sparse(iif,jjf,sA20,ndof,ndof);

if (Flag3D == 1)
	A3 = sparse(iif,jjf,sA3,ndof,ndof);
	A6 = sparse(iif,jjf,sA6,ndof,ndof);
	A9 = sparse(iif,jjf,sA9,ndof,ndof);
	A10 = sparse(iif,jjf,sA10,ndof,ndof);
	A11 = sparse(iif,jjf,sA11,ndof,ndof);
	A12 = sparse(iif,jjf,sA12,ndof,ndof);
	A15 = sparse(iif,jjf,sA15,ndof,ndof);
	A18 = sparse(iif,jjf,sA18,ndof,ndof);
	A21 = sparse(iif,jjf,sA21,ndof,ndof);
end

%% Construct the algebraic equation matrices
ZeroF = sparse(ndof,ndof);

if (Flag3D == 1)
	stabAC1 = [A1+A2+A3];
	stabAC2 = [A4+A5+A6+A7+A8+A9+A10+A11+A12];
	stabAC3 = [A13+A14+A15];
	stabAC4 = [A16+A17+A18];
	stabAC5 = [A19+A20+A21];
else
	stabAC1 = [A1+A2];
	stabAC2 = [A4+A5+A7+A8];
	stabAC3 = [A13+A14];
	stabAC4 = [A16+A17];
	stabAC5 = [A19+A20];
end

%% Assemble left and right handed contributions
% Left-hand side matrix
stabAC1Mat = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*stabAC1;
LHS = LHS + stabAC1Mat + stabAC2 + fluid.gamma * stabAC3;

clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 sA13 sA14 sA15
clear sA16 sA17 sA18 sA19 sA20 sA21
clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear Aij_13 Aij_14 Aij_15 Aij_16 Aij_17 Aij_18 Aij_19 Aij_20 Aij_21
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18 A19 A20 A21
        
% Right-hand side vector
vapFracAlpha = Sol.vapFracAlpha(:,1) ;
vapFracDotAlpha = Sol.vapFracDotAlpha(:,1);

RHS = RHS -(stabAC1) * vapFracDotAlpha(:) - (stabAC2 + fluid.gamma * stabAC4) * vapFracAlpha(:)  ;
RHS = RHS + fluid.gamma *(stabAC5) * ones(ndof,1) ;                                   
                                      
end