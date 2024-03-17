%% 
function [Sol, cavNormIndicator] = cavitation(Sol, fluid, bc, solver, BCTop, BCBottom,...
								  BCLeft, BCRight, pmc, cnn, crd, elemType, nen,...
								  ndof, nElem, Flag3D)                                          
                                           
% Interpolate for alpha values for velocity
Sol.uAlpha = Sol.uPrev + pmc.alpha.*(Sol.u - Sol.uPrev) ;

% Quadrature rules for elements
[gP, gW, N, Nx, Ny, Nz, nQuad] = defineQuadratureData(elemType, Flag3D);

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
% % % % bcLeft = unique(BCLeft(:));
% % % % bcRight = unique(BCRight(:));
% % % % bcTop = unique(BCTop(:));
% % % % bcBottom = unique(BCBottom(:));
% % % % BCBound = [bcLeft; bcRight; bcTop; bcBottom];
% % % % BCBound = unique(BCBound(:));
% % % % disp(unique(BCCylinder(:)),1) = 0.0 ;
% % % % disp(unique(BCCylinder(:)),2) = 0.0 ;
% % % % disp(BCBound,1) = zeros(size(BCBound,1),1);
% % % % disp(BCBound,2) = zeros(size(BCBound,1),1);
% % % % bcStruct = unique(BCStructure(:));
% % % % disp(bcStruct,1) = Sol.dispS(bcStruct,1);
% % % % disp(bcStruct,2) = Sol.dispS(bcStruct,2);
        
%% Allen-Cahn equations
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
    
	if (Flag3D == 1)
		zzf(:,i) = crd(cnn(:,i),3);	
		uz(:,i) =  Sol.uAlpha(cnn(:,i),3,1) ;
    end
   
	pressure(:,i) = Sol.p(cnn(:,i),1) ;
    
    vapFrac(:,i) = Sol.vapFracAlpha(cnn(:,i),1) ;
	vapFracPrev(:,i) = Sol.vapFracPrev(cnn(:,i),1) ;
	vapFracDot(:,i) = Sol.vapFracDotAlpha(cnn(:,i),1) ;
	
    dens(:,i) = fluid.dens(cnn(:,i),1) ;
    dens1(:,i) = fluid.dens1(cnn(:,i),1) ;
	dens2(:,i) = fluid.dens2(cnn(:,i),1) ;	
end
        
%% Form element matrix and assemble Galerkin &  Petrov-Galerkin terms
 [LHS, RHS] = GalerkinTermsCav(Sol, fluid, vapFrac, pressure, dens1, dens2, dens, pmc, solver, iif, jjf, xxf, yyf, zzf, ux, uy, uz, ...
                                    gW, N, Nx, Ny, Nz, nElem, nQuad, nen, ndof, Flag3D);
                                
%  [LHS, RHS] = PetrovGalerkinTermsCav(Sol, fluid, vapFrac, pressure, dens1, dens2, dens, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D);
%
%  [LHS, RHS] = ppvTermsCav(Sol, fluid, vapFrac, vapFracDot, pressure, dens1, dens2, dens, pmc, solver, iif, jjf, ...
%                                           xxf, yyf, ux, uy, gW, N, Nx, Ny, ...
%                                           nElem, nQuad, nen, ndof, elemType, ...
%                                           LHS, RHS, Flag3D);							                                        
      
%% Select the unknown nodal values
zerotypenodes = find(Sol.type == 0);
freeNodes = setdiff(1:size(crd,1),zerotypenodes); % setdiff(1:size(crd,1),[bc.top.nodes';bc.bottom.nodes';zerotypenodes]);
freeNodes = [freeNodes'];

result = Sol.vapFracAlpha(:,1);
result = result(:);

resultDot = Sol.vapFracDotAlpha(:,1);
resultDot = resultDot(:);

%% Solve the linear system
% increment = bicgstab(LHS(freeNodes,freeNodes),RHS(freeNodes),1e-6,2000);
increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);

%% Update the cavitation variables increment
result(freeNodes) = result(freeNodes) + increment;
resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*increment ;

Sol.vapFracAlpha(:,:,1) = reshape(result(1:ndof),[],1);
Sol.vapFracDotAlpha(:,:,1) = reshape(resultDot(1:ndof),[],1);

%% Update the step's output solution
Sol.vapFrac = Sol.vapFracPrev + (1/pmc.alpha)*( Sol.vapFracAlpha - Sol.vapFracPrev ) ;
Sol.vapFracDot = Sol.vapFracDotPrev + (1/pmc.alphaM)*( Sol.vapFracDotAlpha - Sol.vapFracDotPrev ) ;

vapFracMax = 1.0;
vapFracMin = fluid.alphaNuc;

% Clip to vapFracMax and vapFracMin
Sol.vapFrac(Sol.vapFrac>vapFracMax) = vapFracMax;
Sol.vapFrac(Sol.vapFrac<vapFracMin) = vapFracMin;

%% Based on other solutions compute fluid mass distribution
mV = [];
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
	localvapFrac  = sum(repmat(N(p,:),nElem,1).*vapFrac,2);
    localDens = sum(repmat(N(p,:),nElem,1).*dens,2);
	mVij = gW(p)*localvapFrac.*volume.*localDens ;

	mV(1:nElem,p) = mVij ;
end

mV = sum(sum(mV,2));
MassvapFrac = mV ;
Sol.massV = mV ;

clear mV mVij

%% Check the solution conversion and conditioning
cavNormIndicator = norm(increment)/norm(result(freeNodes)) ;

fprintf('Cavitation: %e, Mass: %e', cavNormIndicator, MassvapFrac);
clear result resultDot
end

