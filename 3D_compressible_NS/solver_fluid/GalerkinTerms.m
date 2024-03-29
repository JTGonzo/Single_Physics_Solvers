% Galerkin terms for Navier-Stokes equations

function [LHS, RHS] = GalerkinTerms(Sol, fluid, bc, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, vapFrac, dens, visc, pressure, dens1, dens2, pmc, solver, iif, jjf, xxf, yyf, zzf,...
                            ux, uy, uz, gW, N, Nx, Ny, Nz, nElem, nQuad, nen, ndof, crd, cnn, timeStep, Flag3D)

%% Initialize the quadrature solution space 
sA = [] ;
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

%% Map to local element space
for p = 1:nQuad  
    
	if (Flag3D ~= 1)
	   J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
			yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
		if size(J,2)==1
			J = J';
		end
		volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
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
	end
           
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end

	if (Flag3D ~= 1)
		DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
		DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);

		locUX = sum(repmat(N(p,:),nElem,1).*ux,2);
		locUY = sum(repmat(N(p,:),nElem,1).*uy,2);   
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
	end
	
	localvapFrac  = sum(repmat(N(p,:),nElem,1).*vapFrac,2);
	localP = sum(repmat(N(p,:),nElem,1).*pressure,2);
	localDens1 = sum(repmat(N(p,:),nElem,1).*dens1,2);
	localDens2 = sum(repmat(N(p,:),nElem,1).*dens2,2);	
	localDens = sum(repmat(N(p,:),nElem,1).*dens,2);
	localVisc = sum(repmat(N(p,:),nElem,1).*visc,2);	
            
    index = 0;
    for i = 1:nen
        for j = 1:nen
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*localDens ;
            Mij = Mij.*volume;
            sA(index+1:index+nElem,p) = Mij;
                    
            % Galerkin convection term
            Aij_1 = gW(p)*N(p,i)*(locUX.*DNDx(:,j));
            Aij_2 = gW(p)*N(p,i)*(locUY.*DNDy(:,j));
            Aij_1 = Aij_1.*volume.*localDens;
            Aij_2 = Aij_2.*volume.*localDens;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
                    
            % Galerkin viscous diffusion term
            Aij_4 = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Aij_5 = gW(p)*(DNDy(:,i).*DNDy(:,j));
            Aij_7 = gW(p)*(DNDx(:,i).*DNDy(:,j));
            Aij_4 = Aij_4.*volume.*localVisc;
            Aij_5 = Aij_5.*volume.*localVisc;
            Aij_7 = Aij_7.*volume.*localVisc;
            sA4(index+1:index+nElem,p) = Aij_4;
            sA5(index+1:index+nElem,p) = Aij_5;
            sA7(index+1:index+nElem,p) = Aij_7;
                    
            % Galerkin pressure term (gradP)
            Aij_10 = gW(p)*(DNDx(:,i).*N(p,j));
            Aij_11 = gW(p)*(DNDy(:,i).*N(p,j));
            Aij_10 = Aij_10.*volume;
            Aij_11 = Aij_11.*volume; 
            sA10(index+1:index+nElem,p) = Aij_10;
            sA11(index+1:index+nElem,p) = Aij_11;
                    
            % Galerkin source term (gravity force)
            Aij_13 = gW(p)*(N(p,i)*N(p,j));
            Aij_13 = Aij_13.*localDens.*volume ;
            sA13(index+1:index+nElem,p) = Aij_13;
                    
            % Galerkin continuity term
            Aij_14 = gW(p)*(N(p,i)).*(DNDx(:,j)) ;
            Aij_15 = gW(p)*(N(p,i)).*(DNDy(:,j)) ;
            Aij_14 = Aij_14.*volume ;
            Aij_15 = Aij_15.*volume ;
            sA14(index+1:index+nElem,p) = Aij_14;
            sA15(index+1:index+nElem,p) = Aij_15;
			
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
			
			if(Flag3D == 1)
				rB = ((3/(4*pi*n0))*(1 + alphaNuc - localvapFrac)./localvapFrac).^(1/3);
			else % mod for 2D bubble
				rB = ((1/(pi*n0))*(1 + alphaNuc - localvapFrac)./localvapFrac).^(1/2);
            end		            
			
			k1 = cC*(3*rhol.*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol));
			k2 = cV*(3*rhol.*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol));		

            ind1 = find(localP > Pv);
			ind2 = find(localP < Pv);
			ind3 = find(localP == Pv);
			
			% Galerkin continuity source term  															
			source(ind1) = k1(ind1).*localvapFrac(ind1).*(1-localvapFrac(ind1)).*(localP(ind1)-Pv).^(1/2);
			source(ind2) = k2(ind2).*localvapFrac(ind2).*(1+alphaNuc-localvapFrac(ind2)).*(Pv-localP(ind2)).^(1/2);
            source(ind3) = 0;
            source = source';
			source = (1./rhol-1./rhov).*source;            

			Aij_17 = gW(p)*(source).*(N(p,i)*N(p,j));
			Aij_17 = Aij_17.*volume;					
			sA17(index+1:index+nElem,p) = Aij_17;	

			% Galerkin continuity source term Jacobian	
            sourceJac(ind1) = (k1(ind1)/2).*localvapFrac(ind1).*(1-localvapFrac(ind1)).*(localP(ind1) - Pv).^(-1/2);
            sourceJac(ind2) = (k2(ind2)/2).*localvapFrac(ind2).*(1+alphaNuc-localvapFrac(ind2)).*(Pv - localP(ind2)).^(-1/2);
            sourceJac(ind3) = 0;
            sourceJac = sourceJac';
            sourceJac = (1./rhol-1./rhov).*sourceJac;        
            
			Aij_18 = gW(p)*(sourceJac).*(N(p,i)*N(p,j));
			Aij_18 = Aij_18.*volume;					
			sA18(index+1:index+nElem,p) = Aij_18;	
			
			if (Flag3D == 1)
				Aij_3 = gW(p)*(locUZ.*DNDz(:,j)*N(p,i));				
				Aij_6 = gW(p)*(DNDz(:,i).*DNDz(:,j));
				Aij_8 = gW(p)*(DNDy(:,i).*DNDz(:,j));
				Aij_9 = gW(p)*(DNDx(:,i).*DNDz(:,j));
                Aij_12 = gW(p)*(DNDz(:,i).*N(p,j));
                Aij_16 = gW(p)*(N(p,i)).*(DNDz(:,j));
                Aij_3 = Aij_3.*volume.*localDens;
				Aij_6 = Aij_6.*volume.*localVisc;
				Aij_8 = Aij_8.*volume.*localVisc;
				Aij_9 = Aij_9.*volume.*localVisc;
                Aij_12 = Aij_12.*volume;
                Aij_16 = Aij_16.*volume ;
                sA3(index+1:index+nElem,p) = Aij_3;
				sA6(index+1:index+nElem,p) = Aij_6;
				sA8(index+1:index+nElem,p) = Aij_8;
				sA9(index+1:index+nElem,p) = Aij_9;
				sA12(index+1:index+nElem,p) = Aij_12;			
				sA16(index+1:index+nElem,p) = Aij_16;
			end			
            index = index + nElem;         
            clear sourceJac source
        end
    end
end

%% Summation of all quadrature data
sA = sum(sA,2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA7 = sum(sA7,2);
sA10 = sum(sA10,2);
sA11 = sum(sA11,2);
sA13 = sum(sA13,2);
sA14 = sum(sA14,2);
sA15 = sum(sA15,2);
sA17 = sum(sA17,2);
sA18 = sum(sA18,2);

if (Flag3D == 1)
    sA3 = sum(sA3,2);
    sA6 = sum(sA6,2);
    sA8 = sum(sA8,2);
    sA9 = sum(sA9,2);
    sA12 = sum(sA12,2);
    sA16 = sum(sA16,2);
end
 
%% Assemble the equation matricies element spaces
Mf = sparse(iif,jjf,sA,ndof,ndof);  
ZeroF = sparse(ndof,ndof);

A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);
A10 = sparse(iif,jjf,sA10,ndof,ndof);
A11 = sparse(iif,jjf,sA11,ndof,ndof);
A13 = sparse(iif,jjf,sA13,ndof,ndof);
A14 = sparse(iif,jjf,sA14,ndof,ndof);
A15 = sparse(iif,jjf,sA15,ndof,ndof);
A17 = sparse(iif,jjf,sA17,ndof,ndof);
A18 = sparse(iif,jjf,sA18,ndof,ndof);

if (Flag3D == 1)
    A3 = sparse(iif,jjf,sA3,ndof,ndof);
    A6 = sparse(iif,jjf,sA6,ndof,ndof);
    A8 = sparse(iif,jjf,sA8,ndof,ndof);
    A9 = sparse(iif,jjf,sA9,ndof,ndof);
    A12 = sparse(iif,jjf,sA12,ndof,ndof);
    A16 = sparse(iif,jjf,sA16,ndof,ndof);
end

%% Construct the algebraic equation matrices
if (Flag3D == 1)
	Mf = [Mf ZeroF ZeroF ZeroF;...
		  ZeroF Mf ZeroF ZeroF;...
		  ZeroF ZeroF Mf ZeroF;...
		  ZeroF ZeroF ZeroF ZeroF];
      
	Conv = [A1+A2+A3 ZeroF ZeroF ZeroF;...
			ZeroF A1+A2+A3 ZeroF ZeroF;...
			ZeroF ZeroF A1+A2+A3 ZeroF;...
			ZeroF ZeroF ZeroF ZeroF];
        
	Kf = [2*A4+A5+A6-2/3*A4 A7'-2/3*A7 A9'-2/3*A9 ZeroF;...
		  A7-2/3*A7' A4+2*A5+A6-2/3*A5 A8'-2/3*A8 ZeroF;...
		  A9-2/3*A9' A8-2/3*A8' A4+A5+2*A6-2/3*A6 ZeroF;...
		  ZeroF ZeroF ZeroF ZeroF];
      
% 	Kf = [2*A4+A5+A6 A7' A9' ZeroF;...
% 		  A7 A4+2*A5+A6 A8' ZeroF;...
% 		  A9 A8 A4+A5+2*A6 ZeroF;...
% 		  ZeroF ZeroF ZeroF ZeroF];      

	Fp = [ZeroF ZeroF ZeroF A10;...
		  ZeroF ZeroF ZeroF A11;...
		  ZeroF ZeroF ZeroF A12;...
		  ZeroF ZeroF ZeroF ZeroF];
      
	GMassCons = [ZeroF ZeroF ZeroF ZeroF;...
		  ZeroF ZeroF ZeroF ZeroF;...
		  ZeroF ZeroF ZeroF ZeroF;...
		  A14 A15 A16 ZeroF];
      
	Src = [A13.*fluid.gravFrc(1) ZeroF ZeroF ZeroF;...  
		   ZeroF A13.*fluid.gravFrc(2) ZeroF ZeroF;...
		   ZeroF ZeroF A13.*fluid.gravFrc(3) ZeroF;...
		   ZeroF ZeroF ZeroF ZeroF];
       
	SrcC = [ZeroF ZeroF ZeroF ZeroF;...  
		   ZeroF ZeroF ZeroF ZeroF;...
		   ZeroF ZeroF ZeroF ZeroF;...
		   ZeroF ZeroF ZeroF A17];		 
       
	SrcCJac = [ZeroF ZeroF ZeroF ZeroF;...  
		   ZeroF ZeroF ZeroF ZeroF;...
		   ZeroF ZeroF ZeroF ZeroF;...
		   ZeroF ZeroF ZeroF A18];			   		   
else
	Mf = [Mf ZeroF ZeroF;...
		  ZeroF Mf ZeroF;...
		  ZeroF ZeroF ZeroF];
      
	Conv = [A1+A2 ZeroF ZeroF;...
			ZeroF A1+A2 ZeroF;...
			ZeroF ZeroF ZeroF];
        
	Kf = [2*A4+A5 A7' ZeroF;...
		  A7 A4+2*A5  ZeroF;...
		  ZeroF ZeroF ZeroF];
      
	Fp = [ZeroF ZeroF A10;...
		  ZeroF ZeroF A11;...
		  ZeroF ZeroF ZeroF];
      
	GMassCons = [ZeroF ZeroF ZeroF;...
				 ZeroF ZeroF ZeroF;...
				 A14 A15 ZeroF];
             
	Src = [A13.*fluid.gravFrc(1) ZeroF ZeroF;...  
		   ZeroF A13.*fluid.gravFrc(2) ZeroF;...;...
		   ZeroF ZeroF ZeroF];
       
	SrcC = [ZeroF ZeroF ZeroF;...  
		   ZeroF ZeroF ZeroF;...;...
		   ZeroF ZeroF A17];	
       
	SrcCJac = [ZeroF ZeroF ZeroF;...  
				ZeroF ZeroF ZeroF;...;...
				ZeroF ZeroF A18];		   
end		   

%% Assemble left and right handed contributions
% Left-hand side matrix
Mflow = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Mf;
LHS = Mflow + Conv  + Kf - Fp + GMassCons - SrcCJac;
        
% clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
% clear sA10 sA11 sA12 sA13 sA14 sA15 sA16 sA17 sA18
% clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16 Aij_17 Aij_18
% clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A18

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
	tempVec = [zeros(3*ndof,1); ones(ndof,1)];
else
	p1 = [zeros(ndof*2,1);Sol.p];
	gravVec = [ones(2*ndof,1); zeros(ndof,1)];
	tempVec = [zeros(2*ndof,1); ones(ndof,1)];
end

[Length, tracSrc] = IntegratedOutput3D_try2(Sol, fluid, bc, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, crd, cnn, nElem, ndof, nen, Flag3D);

% Right-hand side vector
RHS = -(Mf * uDotAlpha(:)) - (Conv + Kf + GMassCons)* uAlpha(:) ;

% if(timeStep > 1)
  RHS = RHS + (Fp * p1) + (Src)*gravVec + (SrcC)*tempVec - tracSrc;
%     RHS = RHS + (Fp * p1) + (Src)*gravVec + (SrcC)*tempVec;
% else
%       RHS = RHS + (Fp * p1) + (Src)*gravVec - tracSrc;
% end

end