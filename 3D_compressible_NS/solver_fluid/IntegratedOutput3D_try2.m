%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Calculate the integrated output values of the forces on the structure  %
%                                                                         %
%  This function evaluated the integral of the stress term on the surfa-  %
%  ce of the structure, i.e.,                                             %
%                                                                         %
%                   \int_\Gamma (\sigma . n) d\Gamma,                     %
%                                                                         %
%  where \sigma is the fluid stress, \sigma = -pI + \mu(grad U + grad U^T)%
%  p being the pressure, I is the Identity matrix, \mu is the fluid visc- %
%  osity and grad U is the gradient of velocity. n is the normal to the   %
%  structural surface.                                                    %
%                                                                         %
%  The integral is evaluated by reduced quadrature integration on the two- %
%  dimensional element just beside the surface.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Length, tracSrc] = IntegratedOutput3D_try2(Sol, fluid, bc, BCTop, BCBottom, BCLeft, BCRight, BCFront, BCBack, crd, cnn, nElem, ndof, nen, Flag3D)

%% Reoder element connectivity for reduced integration 
% % Get elements corresponding to boundaries
% cnnBCtop = getElemBC(BCTop, cnn, Flag3D);
% cnnBCbottom = getElemBC(BCBottom, cnn, Flag3D);
% cnnBCleft = getElemBC(BCLeft, cnn, Flag3D);
% cnnBCright = getElemBC(BCRight, cnn, Flag3D);
% if(Flag3D)
% 	cnnBCfront = getElemBC(BCFront, cnn, Flag3D);
% 	cnnBCback = getElemBC(BCBack, cnn, Flag3D);
% end
% 
% % Reorder the element points for reduced integration
% [cnnBCrightNew] = reorderElemBC(BCRight, cnnBCright, Flag3D);
% [cnnBCtopNew] = reorderElemBC(BCTop, cnnBCtop, Flag3D);
% [cnnBCbottomNew] = reorderElemBC(BCBottom, cnnBCbottom, Flag3D);
% [cnnBCleftNew] = reorderElemBC(BCLeft, cnnBCleft, Flag3D);
% % [cnnBCrightNew] = reorderElemBC(BCRight, cnnBCright, Flag3D);
%
% if(Flag3D)
% 	[cnnBCfrontNew] = reorderElemBC(BCFront, cnnBCfront, Flag3D);
% 	[cnnBCbackNew] = reorderElemBC(BCBack, cnnBCback, Flag3D);
% end
% 
% cnnBCnew = [cnnBCtopNew; cnnBCbottomNew; cnnBCleftNew; cnnBCrightNew];
% if(Flag3D)
% 	cnnBCnew = [cnnBCnew; cnnBCfrontNew; cnnBCbackNew];
% end

cnnBCnew = [BCTop; BCBottom; BCLeft; BCRight];

if(Flag3D)
	cnnBCnew = [cnnBCnew; BCFront; BCBack];
end

if(~Flag3D)
	cnnBCnew = [cnnBCnew, ones(size(cnnBCnew, 1), 2)];
else
	cnnBCnew = [cnnBCnew, ones(size(cnnBCnew, 1), 4)];
end

%% Define & Quadrature  
% Shape functions, gauss points and weights for Numerical integration
if(~Flag3D)
	% Gauss points
	gP = [-1/sqrt(3), -1.0
		   1/sqrt(3), -1.0] ;
       
	% Gauss weights
	gW = [1.0, 1.0];
    
	% Shape functions
	N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
	N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
	N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
	N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
	% Derivative of shape functions
	Nx(:,1) = -0.25.*(1-gP(:,2)) ;
	Nx(:,2) =  0.25.*(1-gP(:,2)) ;
	Nx(:,3) =  0.25.*(1+gP(:,2)) ;
	Nx(:,4) = -0.25.*(1+gP(:,2)) ;
	Ny(:,1) = -0.25.*(1-gP(:,1)) ;
	Ny(:,2) = -0.25.*(1+gP(:,1)) ;
	Ny(:,3) =  0.25.*(1+gP(:,1)) ;
	Ny(:,4) =  0.25.*(1-gP(:,1)) ;
	
	Nx = Nx' ;
	Ny = Ny' ;
else
	% Gauss points
	gP = ...
      [-5.7735026918962584E-01, -5.7735026918962584E-01, -1
         5.7735026918962584E-01, -5.7735026918962584E-01, -1
        -5.7735026918962584E-01,  5.7735026918962584E-01, -1
         5.7735026918962584E-01,  5.7735026918962584E-01, -1] ;
     
	% Gauss weights
	gW = [1.0, 1.0, 1.0, 1.0] ;
	
    % Shape functions
	N(:,1) = 0.125.*(1-gP(:,1)).*(1-gP(:,2)).*(1-gP(:,3)) ;
	N(:,2) = 0.125.*(1+gP(:,1)).*(1-gP(:,2)).*(1-gP(:,3)) ;
	N(:,3) = 0.125.*(1+gP(:,1)).*(1+gP(:,2)).*(1-gP(:,3)) ;
	N(:,4) = 0.125.*(1-gP(:,1)).*(1+gP(:,2)).*(1-gP(:,3)) ;
	N(:,5) = 0.125.*(1-gP(:,1)).*(1-gP(:,2)).*(1+gP(:,3)) ;
	N(:,6) = 0.125.*(1+gP(:,1)).*(1-gP(:,2)).*(1+gP(:,3)) ;
	N(:,7) = 0.125.*(1+gP(:,1)).*(1+gP(:,2)).*(1+gP(:,3)) ;
	N(:,8) = 0.125.*(1-gP(:,1)).*(1+gP(:,2)).*(1+gP(:,3)) ;		 
	
    % Derivative of shape functions
	Nx(:,1) = -0.125.*(1-gP(:,2)).*(1-gP(:,3)) ;
	Nx(:,2) =  0.125.*(1-gP(:,2)).*(1-gP(:,3)) ;
	Nx(:,3) =  0.125.*(1+gP(:,2)).*(1-gP(:,3)) ;
	Nx(:,4) = -0.125.*(1+gP(:,2)).*(1-gP(:,3)) ;
	Nx(:,5) = -0.125.*(1-gP(:,2)).*(1+gP(:,3)) ;
	Nx(:,6) =  0.125.*(1-gP(:,2)).*(1+gP(:,3)) ;
	Nx(:,7) =  0.125.*(1+gP(:,2)).*(1+gP(:,3)) ;
	Nx(:,8) = -0.125.*(1+gP(:,2)).*(1+gP(:,3)) ; 
	Ny(:,1) = -0.125.*(1-gP(:,1)).*(1-gP(:,3)) ;
	Ny(:,2) = -0.125.*(1+gP(:,1)).*(1-gP(:,3)) ;
	Ny(:,3) =  0.125.*(1+gP(:,1)).*(1-gP(:,3)) ;
	Ny(:,4) =  0.125.*(1-gP(:,1)).*(1-gP(:,3)) ;
	Ny(:,5) = -0.125.*(1-gP(:,1)).*(1+gP(:,3)) ;
	Ny(:,6) = -0.125.*(1+gP(:,1)).*(1+gP(:,3)) ;
	Ny(:,7) =  0.125.*(1+gP(:,1)).*(1+gP(:,3)) ;
	Ny(:,8) =  0.125.*(1-gP(:,1)).*(1+gP(:,3)) ; 
	Nz(:,1) = -0.125.*(1-gP(:,1)).*(1-gP(:,2)) ;
	Nz(:,2) = -0.125.*(1+gP(:,1)).*(1-gP(:,2)) ;
	Nz(:,3) = -0.125.*(1+gP(:,1)).*(1+gP(:,2)) ;
	Nz(:,4) = -0.125.*(1-gP(:,1)).*(1+gP(:,2)) ;
	Nz(:,5) =  0.125.*(1-gP(:,1)).*(1-gP(:,2)) ;
	Nz(:,6) =  0.125.*(1+gP(:,1)).*(1-gP(:,2)) ;
	Nz(:,7) =  0.125.*(1+gP(:,1)).*(1+gP(:,2)) ;
	Nz(:,8) =  0.125.*(1-gP(:,1)).*(1+gP(:,2)) ;
	
	Nx = Nx' ;
	Ny = Ny' ;
	Nz = Nz' ;
end

% Number of quadrature points
nQuad = size(gW,2);
nElemBC = size(cnnBCnew, 1);

% Form the local to global map
iif = zeros(nen^2*nElemBC,1); 
jjf = zeros(nen^2*nElemBC,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElemBC) = double(cnnBCnew(:,i)); 
      jjf(index+1:index+nElemBC) = double(cnnBCnew(:,j));  
      index = index + nElemBC;
   end
end

xxf = zeros(size(cnnBCnew));
yyf = zeros(size(cnnBCnew));
zzf = zeros(size(cnnBCnew));
ux = zeros(size(cnnBCnew));
uy = zeros(size(cnnBCnew));
if (Flag3D == 1)
	uz = zeros(size(cnnBCnew));
end

%% Localize the data to each element
for i=1:nen
	xxf(:,i) = crd(cnnBCnew(:,i),1);
	yyf(:,i) = crd(cnnBCnew(:,i),2);
	ux(:,i) =  Sol.u(cnnBCnew(:,i),1,1) ;
	uy(:,i) =  Sol.u(cnnBCnew(:,i),2,1) ;
	if(Flag3D)
		zzf(:,i) = crd(cnnBCnew(:,i),3);
		uz(:,i) =  Sol.u(cnnBCnew(:,i),3,1) ;
	end
	% pres(:,i) = Sol.p(cnnBCnew(:,i),1) ;
	pres(:,i) = fluid.pInf ;
	dens(:,i) = fluid.dens(cnnBCnew(:,i),1) ;
	visc(:,i) = fluid.visc(cnnBCnew(:,i),1) ;
end

sA1 = zeros(nen^2*nElemBC,nQuad); 
sA2 = zeros(nen^2*nElemBC,nQuad); 

if(Flag3D)
	sA3 = zeros(nen^2*nElemBC,nQuad); 
end

%% Map to local element space
for p = 1:nQuad  
    % Jacobian for each element
	if (Flag3D ~= 1)
	   J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
			yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
		if size(J,2)==1
			J = J';
		end
		vol =( J(:,1).*J(:,4) -J(:,2).*J(:,3) ); % determinant
		volume = sqrt(J(:,1).^2 + J(:,3).^2); % normalizes the normal vector
		
		% Normal evaluation
		normal(:,1) = -J(:,3)./volume ;
		normal(:,2) = J(:,1)./volume ;		
	else
		J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
			 yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
			 zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
		if size(J,2)==1
			J = J';
		end
		vol =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
				  J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
				  J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
		
		% Normal evaluation
		temp1 = (J(:,5).*J(:,7) - J(:,4).*J(:,8));
		temp2 = (J(:,8).*J(:,1) - J(:,7).*J(:,2));
		temp3 = (J(:,2).*J(:,4) - J(:,1).*J(:,5));
		
		volume = sqrt(temp1.^2 + temp2.^2 + temp3.^2);
		normal(:,1) = temp1./volume ;
		normal(:,2) = temp2./volume ;
		normal(:,3) = temp3./volume ;	
		
		clear temp1 temp2 temp3
	end
    vol = abs(vol);
	
    negJacobian = find(vol<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end	

	if (Flag3D ~= 1)
		DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
		DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen); 
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
	end	
    
    % Pressure and gradients evaluation for velocity
    locP  = sum(repmat(N(p,:),size(cnnBCnew, 1),1).*pres,2);
	localVisc = sum(repmat(N(p,:),nElemBC,1).*visc,2);
    
	locgradUx = sum(ux.*DNDx,2) ;
    locgradUy = sum(ux.*DNDy,2) ; 
    locgradVx = sum(uy.*DNDx,2) ;
    locgradVy = sum(uy.*DNDy,2) ;
	
	sigma11 = -locP + 2*localVisc.*locgradUx;
    sigma12 = localVisc.*(locgradUy + locgradVx);
    sigma21 = localVisc.*(locgradVx + locgradUy);
    sigma22 = -locP + 2*localVisc.*locgradVy;
	
	if(Flag3D)
		locgradUz = sum(ux.*DNDz,2) ;
		locgradVz = sum(uy.*DNDz,2) ;
		locgradWx = sum(uz.*DNDx,2) ;
		locgradWy = sum(uz.*DNDy,2) ;
		locgradWz = sum(uz.*DNDz,2) ;
		
		sigma13 = localVisc.*(locgradUz + locgradWx);
		sigma23 = localVisc.*(locgradVz + locgradWy);
		sigma31 = localVisc.*(locgradWx + locgradUz);
		sigma32 = localVisc.*(locgradWy + locgradVz);
		sigma33 = -locP + 2*localVisc.*locgradWz;
	end
            
    % Length/ Area of the line/ surface integral in 1D/2D
     A0(:,p) = gW(p).*volume ;
    
    index = 0;
    for i = 1:nen
        for j = 1:nen 			
			% Galerkin source term (Stress term)			
			if(~Flag3D)
				Aij_1 = gW(p)*N(p,i)*normal(:, 1).*(sigma11 + sigma21)*N(p,j);
				Aij_2 = gW(p)*N(p,i)*normal(:, 2).*(sigma12 + sigma22)*N(p,j);
				Aij_1 = Aij_1.*volume ;
				Aij_2 = Aij_2.*volume ;
				sA1(index+1:index+nElemBC,p) = Aij_1;
				sA2(index+1:index+nElemBC,p) = Aij_2;
			else
				Aij_1 = gW(p)*N(p,i)*normal(:, 1).*(sigma11 + sigma21 + sigma31)*N(p,j);
				Aij_2 = gW(p)*N(p,i)*normal(:, 2).*(sigma12 + sigma22 + sigma32)*N(p,j);
				Aij_3 = gW(p)*N(p,i)*normal(:, 3).*(sigma13 + sigma23 + sigma33)*N(p,j);
				Aij_1 = Aij_1.*volume ;
				Aij_2 = Aij_2.*volume ;
				Aij_3 = Aij_3.*volume ;
				sA1(index+1:index+nElemBC,p) = Aij_1;
				sA2(index+1:index+nElemBC,p) = Aij_2;
				sA3(index+1:index+nElemBC,p) = Aij_3;				
            end            
            index = index + nElemBC;			
        end  
    end
end

%% Summation of all quadrature data
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
A0 = sum(A0,2);
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);

%% Assemble solution vectors
tracSrc1 = [A1]*[ones(1*ndof,1)];
tracSrc2 = [A2]*[ones(1*ndof,1)];
tracSrc = [tracSrc1; tracSrc2];

if(Flag3D)
	sA3 = sum(sA3,2);
	A3 = sparse(iif,jjf,sA3,ndof,ndof);
	tracSrc3 = [A3]*[ones(1*ndof,1)];
	tracSrc = [tracSrc; tracSrc3];
end

tracSrc = [tracSrc; zeros(1*ndof,1)];

Length = [sum(A0,1)];

end
