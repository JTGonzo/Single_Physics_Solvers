% Galerkin terms for Allen-Cahn equations

function [LHS, RHS] = GalerkinTermsCav(Sol, fluid, vapFrac, pressure, dens1, dens2, dens, pmc, solver, iif, jjf, xxf, yyf, zzf, ux, uy, uz, ...
                                    gW, N, Nx, Ny, Nz, nElem, nQuad, nen, ndof, Flag3D)

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
	
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
    end
	
    volume = abs(volume);

	if (Flag3D ~= 1)
		DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
		DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);

		locUX  = sum(repmat(N(p,:),nElem,1).*ux,2);
		locUY  = sum(repmat(N(p,:),nElem,1).*uy,2);   
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
        
    index = 0;
    for i = 1:nen
        for j = 1:nen   
                    % Galerkin inertia term
                    Mij = gW(p)*(N(p,i)*N(p,j));
                    Mij = Mij.*volume;
                    sA(index+1:index+nElem,p) = Mij;
                    
                    % Galerkin convection term
                    Aij_1 = gW(p)*(locUX.*DNDx(:,j)*N(p,i));
                    Aij_2 = gW(p)*(locUY.*DNDy(:,j)*N(p,i));
                    Aij_1 = Aij_1.*volume;
                    Aij_2 = Aij_2.*volume;
                    sA1(index+1:index+nElem,p) = Aij_1;
                    sA2(index+1:index+nElem,p) = Aij_2;
                    
                    % Galerkin reaction term 
                    reacJac = 0.0;
                    Aij_7 = gW(p)*(reacJac).*(N(p,i)*N(p,j));
                    Aij_7 = Aij_7.*volume;
                    sA7(index+1:index+nElem,p) = Aij_7;
                    	
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
                    
                    k1 = cC*(3).*sqrt((2/3)./(rhol));
                    k2 = cV*(3).*sqrt((2/3)./(rhol));
                    
                    ind1 = find(localP > Pv);
                    ind2 = find(localP < Pv);
                    ind3 = find(localP == Pv);                      

					% Galerkin source term 
                    reac = 0.0;
                    Aij_8 = gW(p)*(reac).*(N(p,i)*N(p,j));
                    Aij_8 = Aij_8.*volume;					
                    sA8(index+1:index+nElem,p) = Aij_8;
                    
                    % Galerkin continuity source term  
                    source(ind1) = (k1(ind1)./rB(ind1)).*localvapFrac(ind1).*(1-localvapFrac(ind1)).*(localP(ind1)-Pv).^(1/2);
                    source(ind2) = (k2(ind2)./rB(ind2)).*localvapFrac(ind2).*(1+alphaNuc-localvapFrac(ind2)).*(Pv-localP(ind2)).^(1/2);
                    source(ind3) = 0;
					source = source';					

                    Aij_9 = gW(p)*(source).*(N(p,i)*N(p,j));
                    Aij_9 = Aij_9.*volume;					
                    sA9(index+1:index+nElem,p) = Aij_9;
					
					% Galerkin continuity source term Jacobian	
					sourceJac(ind1) = ((3/(4*pi*n0))^(-1/3))*k1(ind1).*((localP(ind1) - Pv).^(1/2)).*((1-2*localvapFrac(ind1)).*(((1+alphaNuc-localvapFrac(ind1))./localvapFrac(ind1)).^(-1/3)) + (-1/3)*localvapFrac(ind1).*(1-localvapFrac(ind1)).*(((1+alphaNuc-localvapFrac(ind1))./localvapFrac(ind1)).^(-4/3)).*((-1-alphaNuc)./(localvapFrac(ind1)).^2));
                    sourceJac(ind2) = ((3/(4*pi*n0))^(-1/3))*k2(ind2).*((Pv - localP(ind2)).^(1/2)).*((1+alphaNuc-2*localvapFrac(ind2)).*(((1+alphaNuc-localvapFrac(ind2))./localvapFrac(ind2)).^(-1/3)) + (-1/3)*localvapFrac(ind2).*(1+alphaNuc-localvapFrac(ind2)).*(((1+alphaNuc-localvapFrac(ind2))./localvapFrac(ind2)).^(-4/3)).*((-1-alphaNuc)./(localvapFrac(ind2)).^2));
                    sourceJac(ind3) = 0;
                    sourceJac = sourceJac';                    					
					
                    Aij_10 = gW(p)*(sourceJac).*(N(p,i)*N(p,j));
                    Aij_10 = Aij_10.*volume;					
                    sA10(index+1:index+nElem,p) = Aij_10;					
					
                    if (Flag3D == 1)
						Aij_3 = gW(p)*(locUZ.*DNDz(:,j)*N(p,i));
                        Aij_6 = gW(p)*(DNDz(:,i).*DNDz(:,j));
						Aij_3 = Aij_3.*volume;
                        Aij_6 = Aij_6.*volume.*solver.epsilon^2;
						sA3(index+1:index+nElem,p) = Aij_3;
					    sA6(index+1:index+nElem,p) = Aij_6;
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
sA8 = sum(sA8,2);
sA9 = sum(sA9,2);
sA10 = sum(sA10,2);

if (Flag3D == 1)
	sA3 = sum(sA3,2);
	sA6 = sum(sA6,2);
end
 
%% Assemble the equation matricies element spaces   
Mf = sparse(iif,jjf,sA,ndof,ndof);

A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);
A8 = sparse(iif,jjf,sA8,ndof,ndof);
A9 = sparse(iif,jjf,sA9,ndof,ndof);
A10 = sparse(iif,jjf,sA10,ndof,ndof);

if (Flag3D == 1)
	A3 = sparse(iif,jjf,sA3,ndof,ndof);
	A6 = sparse(iif,jjf,sA6,ndof,ndof);
end

%% Construct the algebraic equation matrices
% Left-hand side matrix
if (Flag3D == 1)
	Mf     = [Mf];
	Conv   = [A1+A2+A3];
	Kf     = [A4+A5+A6];
	Rf     = [A7];
	src    = [A9];	
	srcJac = [A10];	
else
	Mf     = [Mf];
	Conv   = [A1+A2];
	Kf     = [A4+A5];
	Rf     = [A7];
	src    = [A9];	
	srcJac = [A10];
end

%% Assemble left and right handed contributions
% Left-hand side matrix
LHS = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Mf;
% LHS = LHS + Conv + Kf + Rf - srcJac;
LHS = LHS + Conv + Rf - srcJac;
        
% Right-hand side vector
vapFracAlpha = Sol.vapFracAlpha(:,1) ;
vapFracDotAlpha = Sol.vapFracDotAlpha(:,1);

% RHS = -(Mf * vapFracDotAlpha(:)) - (Conv + Kf) * vapFracAlpha(:) ;
RHS = -(Mf * vapFracDotAlpha(:)) - (Conv) * vapFracAlpha(:) ;
RHS = RHS - (A8) * vapFracAlpha(:) + (src) * ones(ndof,1) ;

clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
clear sA10 sA11 sA12
clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12

end