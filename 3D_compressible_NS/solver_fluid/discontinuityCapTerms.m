% Petrov-Galerkin stabilization terms for Navier-Stokes equations
% Part 3: (tauM/rho)(gradN).(rho*u.(gradU)) + 
%         (tauM/rho)(gradN).(gradP) +
%         (tauM/rho)(gradN).(-rho.gravVec)
%
function [LHS, RHS] = discontinuityCapTerms(Sol, fluid, vapFrac, dens, visc, pres, pmc, solver, iif, jjf, ...
                                          xxf, yyf, zzf, ux, uy, uz, gW, N, Nx, Ny, Nz, ...
                                          nElem, nQuad, nen, ndof, elemType, ...
                                          LHS, RHS, Flag3D, cnn, crd)

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

%% Initialize the solution variable space
for i=1:nen
	xxf(:,i) = crd(cnn(:,i),1);
	yyf(:,i) = crd(cnn(:,i),2);
	if (Flag3D == 1)
	zzf(:,i) = crd(cnn(:,i),3);
    end
    
    ux(:,i) =  Sol.uAlpha(cnn(:,i),1,1) ;
	uxDot(:,i) = Sol.uDotAlpha(cnn(:,i),1,1) ;
	uy(:,i) =  Sol.uAlpha(cnn(:,i),2,1) ;
	uyDot(:,i) = Sol.uDotAlpha(cnn(:,i),2,1) ;
	if (Flag3D == 1)
	uz(:,i) =  Sol.uAlpha(cnn(:,i),3,1) ;
	uzDot(:,i) = Sol.uDotAlpha(cnn(:,i),3,1) ;
    end
    
    vapFrac(:,i) = Sol.vapFracAlpha(cnn(:,i),1) ;
	pres(:,i) = Sol.p(cnn(:,i),1) ;
end

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
			
%                     c1 =  1.154700538379252 ;
%                     c2 = -5.773502691896259e-01 ;
%                     
%                     a1 = c1 * J(:,1) + c2 * J(:,2) ;
%                     a2 = c1 * J(:,2) + c2 * J(:,1) ;
%                     gijUp(:,1) = J(:,1).* a1 + J(:,2).*a2 ;
%                     
%                     a1 = c1 * J(:,3) + c2 * J(:,4) ;
%                     a2 = c1 * J(:,4) + c2 * J(:,3) ;
%                     gijUp(:,2) = J(:,3).* a1 + J(:,4).* a2 ;
%                     gijUp(:,3) = J(:,1).* a1 + J(:,2).* a2 ; 
			gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) ;
			gijUp(:,2) = J(:,3).*J(:,3) + J(:,4).*J(:,4) ;
			gijUp(:,3) = J(:,1).*J(:,3) + J(:,2).*J(:,4) ;
            
		elseif (strcmp(elemType,'4Quad'))
			gijDown(:,1) = xsInv(:,1).*xsInv(:,1) + xsInv(:,3).*xsInv(:,3) ;
			gijDown(:,2) = xsInv(:,2).*xsInv(:,2) + xsInv(:,4).*xsInv(:,4) ;
			gijDown(:,3) = xsInv(:,2).*xsInv(:,1) + xsInv(:,4).*xsInv(:,3) ;
			gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) ;
			gijUp(:,2) = J(:,3).*J(:,3) + J(:,4).*J(:,4) ;
			gijUp(:,3) = J(:,1).*J(:,3) + J(:,2).*J(:,4) ;
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
			
			gijUp(:,1) = J(:,1).*J(:,1) + J(:,2).*J(:,2) + J(:,3).*J(:,3);
			gijUp(:,2) = J(:,4).*J(:,4) + J(:,5).*J(:,5) + J(:,6).*J(:,6);
			gijUp(:,3) = J(:,7).*J(:,7) + J(:,8).*J(:,8) + J(:,9).*J(:,9);
			gijUp(:,4) = J(:,1).*J(:,4) + J(:,2).*J(:,5) + J(:,3).*J(:,6);
			gijUp(:,5) = J(:,4).*J(:,7) + J(:,5).*J(:,8) + J(:,6).*J(:,9);
			gijUp(:,6) = J(:,1).*J(:,7) + J(:,2).*J(:,8) + J(:,3).*J(:,9);
		end
	end
    
    volume = abs(volume);

	if (Flag3D ~= 1)
		DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
		DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);

		locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
		locUY  = sum(repmat(N(p,:),nElem,1).*uy,2); 
		locUXDot  = sum(repmat(N(p,:),nElem,1).*uxDot,2);
		locUYDot  = sum(repmat(N(p,:),nElem,1).*uyDot,2); 
		
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
		
		locUXDot  = sum(repmat(N(p,:),nElem,1).*uxDot,2);
		locUYDot  = sum(repmat(N(p,:),nElem,1).*uyDot,2);
		locUZDot  = sum(repmat(N(p,:),nElem,1).*uzDot,2);
		
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
	% dens = zeros(size(localvapFrac,1),1) ; % edit - new
	% visc = zeros(size(localvapFrac,1),1) ; % edit - new	
	% dens = fluid.dens1*(localvapFrac) + fluid.dens2*(1-localvapFrac) ; % edit - new
	% visc = fluid.visc1*(localvapFrac) + fluid.visc2*(1-localvapFrac) ; % edit - new	
  
    tauM = (2/solver.dt)^2 + ((visc./dens).^2).* GG * 36 + abs(VGV);    
    tauM = tauM.^(0.5);
    tauM = 1./tauM;
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
			locGradPX  = sum(DNDx.*pres,2);
			locGradPY  = sum(DNDy.*pres,2);
			locGradU(:,1)  = sum(DNDx.*ux,2);
			locGradU(:,2)  = sum(DNDx.*uy,2);
			locGradU(:,4)  = sum(DNDy.*ux,2);
			locGradU(:,5)  = sum(DNDy.*uy,2);
            
			if (Flag3D==1)
				locGradPZ  = sum(DNDz.*pres,2);
				locGradU(:,3)  = sum(DNDx.*uz,2);
				locGradU(:,6)  = sum(DNDy.*uz,2);
				locGradU(:,7)  = sum(DNDz.*ux,2);
				locGradU(:,8)  = sum(DNDz.*uy,2);
				locGradU(:,9)  = sum(DNDz.*uz,2);
				ResU(:,1) = dens.*locUXDot + dens.*(locUX.*locGradU(:,1) + locUY.*locGradU(:,4) + ...
				 locUZ.*locGradU(:,7)) + locGradPX - dens.*fluid.gravFrc(1) ;
				ResU(:,2) = dens.*locUYDot + dens.*(locUX.*locGradU(:,2) + locUY.*locGradU(:,5) + ...
				 locUZ.*locGradU(:,8)) + locGradPY - dens.*fluid.gravFrc(2) ;
				ResU(:,3) = dens.*locUZDot + dens.*(locUX.*locGradU(:,3) + locUY.*locGradU(:,6) + ...
				 locUZ.*locGradU(:,9)) + locGradPZ - dens.*fluid.gravFrc(3) ;
			else
				ResU(:,1) = dens.*locUXDot + dens.*(locUX.*locGradU(:,1) + locUY.*locGradU(:,4)) + ...
				  locGradPX - dens.*fluid.gravFrc(1) ;
				ResU(:,2) = dens.*locUYDot + dens.*(locUX.*locGradU(:,2) + locUY.*locGradU(:,5)) + ...
				  locGradPY - dens.*fluid.gravFrc(2) ; 
			end
	
			if (Flag3D==1)
				gGradV(:,1) = gijUp(:,1).*locGradU(:,1) + gijUp(:,4).*locGradU(:,4) + gijUp(:,6).*locGradU(:,7) ;
				gGradV(:,2) = gijUp(:,4).*locGradU(:,1) + gijUp(:,2).*locGradU(:,4) + gijUp(:,5).*locGradU(:,7) ;
				gGradV(:,3) = gijUp(:,6).*locGradU(:,1) + gijUp(:,5).*locGradU(:,4) + gijUp(:,3).*locGradU(:,7) ;
				gGradV(:,4) = gijUp(:,1).*locGradU(:,2) + gijUp(:,4).*locGradU(:,5) + gijUp(:,6).*locGradU(:,8) ;
				gGradV(:,5) = gijUp(:,4).*locGradU(:,2) + gijUp(:,2).*locGradU(:,5) + gijUp(:,5).*locGradU(:,8) ;
				gGradV(:,6) = gijUp(:,6).*locGradU(:,2) + gijUp(:,5).*locGradU(:,5) + gijUp(:,3).*locGradU(:,8) ;
				gGradV(:,7) = gijUp(:,1).*locGradU(:,3) + gijUp(:,4).*locGradU(:,6) + gijUp(:,6).*locGradU(:,9) ;
				gGradV(:,8) = gijUp(:,4).*locGradU(:,3) + gijUp(:,2).*locGradU(:,6) + gijUp(:,5).*locGradU(:,9) ;
				gGradV(:,9) = gijUp(:,6).*locGradU(:,3) + gijUp(:,5).*locGradU(:,6) + gijUp(:,3).*locGradU(:,9) ;
		
				dcFct = locGradU(:,1).*gGradV(:,1) + locGradU(:,4).*gGradV(:,2) ...
					+ locGradU(:,7).*gGradV(:,3) + locGradU(:,2).*gGradV(:,4) ...
					+ locGradU(:,5).*gGradV(:,5) + locGradU(:,8).*gGradV(:,6) ...
					+ locGradU(:,3).*gGradV(:,7) + locGradU(:,6).*gGradV(:,8) ...
					+ locGradU(:,9).*gGradV(:,9) +eps;
			else
				gGradV(:,1) = gijUp(:,1).*locGradU(:,1) + gijUp(:,3).*locGradU(:,4)  ;
				gGradV(:,2) = gijUp(:,3).*locGradU(:,1) + gijUp(:,2).*locGradU(:,4)  ;
				gGradV(:,4) = gijUp(:,1).*locGradU(:,2) + gijUp(:,3).*locGradU(:,5)  ;
				gGradV(:,5) = gijUp(:,3).*locGradU(:,2) + gijUp(:,2).*locGradU(:,5)  ;
		
				dcFct = locGradU(:,1).*gGradV(:,1) + locGradU(:,4).*gGradV(:,2) ...
					+ locGradU(:,2).*gGradV(:,4) + locGradU(:,5).*gGradV(:,5) + eps ;  
			end
	  
			dcFct(dcFct>eps) = 1./dcFct(dcFct>eps) ;
			dcFct(dcFct<=eps) = 0.0 ;
            
			if (Flag3D==1)
				dcFct = dcFct.*( ResU(:,1).*ResU(:,1) + ResU(:,2).*ResU(:,2) + ResU(:,3).*ResU(:,3) ) ;
			else
				dcFct = dcFct.*( ResU(:,1).*ResU(:,1) + ResU(:,2).*ResU(:,2) ) ;
            end
            
            dcQuad = tauM.*dcFct./(dens.^2) ;
			dcConst = 1./(3 * tauM) ;
			dcFct = min( dcConst, 2.*dcQuad ) ;
				
			if (Flag3D == 1)
				Aij_1 = gW(p)*(DNDx(:,i).*gijUp(:,1).*DNDx(:,j));
				Aij_2 = gW(p)*(DNDx(:,i).*gijUp(:,4).*DNDy(:,j));
				Aij_4 = gW(p)*(DNDy(:,i).*gijUp(:,4).*DNDx(:,j));
				Aij_5 = gW(p)*(DNDy(:,i).*gijUp(:,2).*DNDy(:,j));
                Aij_3 = gW(p)*(DNDx(:,i).*gijUp(:,6).*DNDz(:,j));
				Aij_6 = gW(p)*(DNDy(:,i).*gijUp(:,5).*DNDz(:,j));
				Aij_7 = gW(p)*(DNDz(:,i).*gijUp(:,6).*DNDx(:,j));
				Aij_8 = gW(p)*(DNDz(:,i).*gijUp(:,5).*DNDy(:,j));
				Aij_9 = gW(p)*(DNDz(:,i).*gijUp(:,3).*DNDz(:,j));
				Aij_1 = Aij_1.*volume.*dens.*dcFct;
				Aij_2 = Aij_2.*volume.*dens.*dcFct;
				Aij_4 = Aij_4.*volume.*dens.*dcFct;
				Aij_5 = Aij_5.*volume.*dens.*dcFct;
                Aij_3 = Aij_3.*volume.*dens.*dcFct;
				Aij_6 = Aij_6.*volume.*dens.*dcFct;
				Aij_7 = Aij_7.*volume.*dens.*dcFct;
				Aij_8 = Aij_8.*volume.*dens.*dcFct;
				Aij_9 = Aij_9.*volume.*dens.*dcFct;
				sA1(index+1:index+nElem,p) = Aij_1;
				sA2(index+1:index+nElem,p) = Aij_2;
				sA4(index+1:index+nElem,p) = Aij_4;
				sA5(index+1:index+nElem,p) = Aij_5;
				sA3(index+1:index+nElem,p) = Aij_3;
				sA6(index+1:index+nElem,p) = Aij_6;
				sA7(index+1:index+nElem,p) = Aij_7;
				sA8(index+1:index+nElem,p) = Aij_8;
				sA9(index+1:index+nElem,p) = Aij_9;
			else
				Aij_1 = gW(p)*(DNDx(:,i).*gijUp(:,1).*DNDx(:,j));
				Aij_2 = gW(p)*(DNDx(:,i).*gijUp(:,3).*DNDy(:,j));
				Aij_4 = gW(p)*(DNDy(:,i).*gijUp(:,3).*DNDx(:,j));
				Aij_5 = gW(p)*(DNDy(:,i).*gijUp(:,2).*DNDy(:,j));
				Aij_1 = Aij_1.*volume.*dens.*dcFct;
				Aij_2 = Aij_2.*volume.*dens.*dcFct;
				Aij_4 = Aij_4.*volume.*dens.*dcFct;
				Aij_5 = Aij_5.*volume.*dens.*dcFct;
				sA1(index+1:index+nElem,p) = Aij_1;
				sA2(index+1:index+nElem,p) = Aij_2;
				sA4(index+1:index+nElem,p) = Aij_4;
				sA5(index+1:index+nElem,p) = Aij_5;    
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

if (Flag3D == 1)
	sA3 = sum(sA3,2);
	sA6 = sum(sA6,2);
	sA7 = sum(sA7,2);
	sA8 = sum(sA8,2);
	sA9 = sum(sA9,2);
end
        
%% Assemble the equation matricies element spaces    
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);

if (Flag3D == 1)
	A3 = sparse(iif,jjf,sA3,ndof,ndof);
	A6 = sparse(iif,jjf,sA6,ndof,ndof);
	A7 = sparse(iif,jjf,sA7,ndof,ndof);
	A8 = sparse(iif,jjf,sA8,ndof,ndof);
	A9 = sparse(iif,jjf,sA9,ndof,ndof);
end

ZeroF = sparse(ndof,ndof);

%% Construct the algebraic equation matrices
if (Flag3D == 1)
    stabDC = [A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF ZeroF;...
		  ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF;...
		  ZeroF ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF;...
		  ZeroF ZeroF ZeroF ZeroF];
else
    stabDC = [A1+A2+A4+A5 ZeroF ZeroF;...
		  ZeroF A1+A2+A4+A5 ZeroF;...
		  ZeroF ZeroF ZeroF];
end
 
%% Assemble left and right handed contributions
% Left-hand side matrix
LHS = LHS + (stabDC);
        
clear sA1 sA2 sA3 sA4 sA5 sA6 sA7 sA8 sA9 sA10 sA11 sA12 
clear Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12

uAlpha = Sol.uAlpha(:,:,1) ;
uAlpha = [uAlpha(:); zeros(ndof,1)];
uDotAlpha = Sol.uDotAlpha(:,:,1) ;
uDotAlpha = [uDotAlpha(:); zeros(ndof,1)];

% Right-hand side vector
RHS = RHS - stabDC *uAlpha(:) ;                                                        
end