
function multiPhaseFem_partitioned_AFEM

%% Time loop starts here
for timeStep = 1:solver.maxSteps
    
    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax
 
        % Galerkin terms for Navier-Stokes equations

        for p = 1:nQuad  
            
            index = 0;
            for i = 1:nen
                for j = 1:nen                  
					
                    % Galerkin source term (continuity source)
					rhov = fluid.dens2;
					rhol = fluid.dens1;
					rhom = rhol*localvapFracPrev + rhov*(1-localvapFracPrev);
					Pv = fluid.pVap;
                    alphaNuc = fluid.alphaNuc;
					dTime = solver.dt;
					rB = ((3/(4*pi*n0))*(1 + alphaNuc - localvapFracPrev)./localvapFracPrev).^(1/3);
                    
					source1 = cC*localvapFracPrev.*(1-localvapFracPrev).*(3*rhol*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol*abs(localP-Pv))).*max(localP-Pv, 0);
					source2 = cV*localvapFracPrev.*(1+alphaNuc-localvapFracPrev).*(3*rhol*rhov./(rhom.*rB)).*sqrt((2/3)./(rhol*abs(localP-Pv))).*min(localP-Pv, 0);

                    %source = (1/rhol-1/rhov)*(source1 + source2);
					
					source = 0.0;
                   
                    Aij_17 = gW(p)*(source).*(N(p,i)*N(p,j));
                    Aij_17 = Aij_17.*volume;					
                    sA17(index+1:index+nElem,p) = Aij_17;					

                    index = index + nElem;
                end
            end
            
            
            
              % Surface tension
        
           if (Flag3D ~= 1)
                gradXvapFrac = sum(DNDx.*vapFrac,2);
                gradYvapFrac = sum(DNDy.*vapFrac,2);
                gradPNS = gradXvapFrac.^2+gradYvapFrac.^2;
            else
                gradXvapFrac = sum(DNDx.*vapFrac,2);
                gradYvapFrac = sum(DNDy.*vapFrac,2);
                gradZvapFrac = sum(DNDz.*vapFrac,2);
                gradPNS = gradXvapFrac.^2+gradYvapFrac.^2+gradZvapFrac.^2;
            end
            
            if (Flag3D ~=1)
                st_x(:,:,p) = gW(p)*DNDx.*repmat(gradPNS,1,nen).*repmat(volume,1,nen);
                st_y(:,:,p) = gW(p)*DNDy.*repmat(gradPNS,1,nen).*repmat(volume,1,nen);
            else
                st_x(:,:,p) = gW(p)*DNDx.*repmat(gradPNS,1,nen).*repmat(volume,1,nen);
                st_y(:,:,p) = gW(p)*DNDy.*repmat(gradPNS,1,nen).*repmat(volume,1,nen);
                st_z(:,:,p) = gW(p)*DNDz.*repmat(gradPNS,1,nen).*repmat(volume,1,nen);
            end
            
            if (Flag3D ~=1)
                std_x(:,:,p) = gW(p)*(DNDx.*repmat(gradXvapFrac.*gradXvapFrac,1,nen)+DNDy.*repmat(gradYvapFrac.*gradXvapFrac,1,nen)).*repmat(volume,1,nen);
                std_y(:,:,p) = gW(p)*(DNDx.*repmat(gradXvapFrac.*gradYvapFrac,1,nen)+DNDy.*repmat(gradYvapFrac.*gradYvapFrac,1,nen)).*repmat(volume,1,nen);
            else
                std_x(:,:,p) = gW(p)*(DNDx.*repmat(gradXvapFrac.*gradXvapFrac,1,nen)+DNDy.*repmat(gradYvapFrac.*gradXvapFrac,1,nen)+DNDz.*repmat(gradZvapFrac.*gradXvapFrac,1,nen)).*repmat(volume,1,nen);
                std_y(:,:,p) = gW(p)*(DNDx.*repmat(gradXvapFrac.*gradYvapFrac,1,nen)+DNDy.*repmat(gradYvapFrac.*gradYvapFrac,1,nen)+DNDz.*repmat(gradZvapFrac.*gradYvapFrac,1,nen)).*repmat(volume,1,nen);
                std_z(:,:,p) = gW(p)*(DNDx.*repmat(gradXvapFrac.*gradZvapFrac,1,nen)+DNDy.*repmat(gradYvapFrac.*gradZvapFrac,1,nen)+DNDz.*repmat(gradZvapFrac.*gradZvapFrac,1,nen)).*repmat(volume,1,nen);
            end
            
            if (Flag3D ~=1)
                sTx  = sum(st_x,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTy  = sum(st_y,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTdx = sum(std_x,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTdy = sum(std_y,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                STx  = sparse(cnn,1,sTx,ndof,1);
                STy  = sparse(cnn,1,sTy,ndof,1);
                STDx = sparse(cnn,1,sTdx,ndof,1);
                STDy = sparse(cnn,1,sTdy,ndof,1);
                STP  = [STx;STy;zeros(ndof,1)];
                STD  = [STDx;STDy;zeros(ndof,1)];
            else
                sTx  = sum(st_x,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTy  = sum(st_y,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTz  = sum(st_z,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                
                sTdx = sum(std_x,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTdy = sum(std_y,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                sTdz = sum(std_z,3)*fluid.surfTens*solver.epsilon*fluid.alphasf;
                
                STx  = sparse(cnn,1,sTx,ndof,1);
                STy  = sparse(cnn,1,sTy,ndof,1);
                STz  = sparse(cnn,1,sTz,ndof,1);
                STDx = sparse(cnn,1,sTdx,ndof,1);
                STDy = sparse(cnn,1,sTdy,ndof,1);
                STDz = sparse(cnn,1,sTdz,ndof,1);
                STP  = [STx;STy;STz;zeros(ndof,1)];
                STD  = [STDx;STDy;STDz;zeros(ndof,1)];
            end
            
        end

		SrcC = [A17 ZeroF ZeroF;...  
               ZeroF A17 ZeroF;...;...
               ZeroF ZeroF ZeroF];	   
			   
        end
   
        RHS = RHS + (Fp * p1) + (Src)*gravVec + SrcC*gravVec - STP + STD;
        
        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 1: (tauM/rho)(rho*u.gradN).(rho*dU/dt) + 
        %         (tauM/rho)(rho*u.gradN).(rho*u.gradU) +
        %         (tauM/rho)(gradN).(rho*dU/dt)
        %

        for p = 1:nQuad  
            
            index = 0;
            for i = 1:nen
                for j = 1:nen

      


                    index = index + nElem;
                end
            end
        end



        
        % Assemble the matrix        
        

        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 2: (grad.N).tauC.(rho*(grad.U)) + 
        %         (tauM/rho)(rho*u.gradN).(gradP) +
        %         (tauM/rho)(-rho*gravFrc)

        for p = 1:nQuad  

            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term
                    
                    % 2nd term
                    
                    % 3rd term   

                    index = index + nElem;
                end
            end
        end  
   
        % Assemble the matrix              
        
        % Petrov-Galerkin stabilization terms for Navier-Stokes equations
        % Part 3: (tauM/rho)(gradN).(rho*u.(gradU)) + 
        %         (tauM/rho)(gradN).(gradP) +
        %         (tauM/rho)(gradN).(-rho.gravVec)

        for p = 1:nQuad  
            
            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % 1st term

                    % 2nd term
                    
                    % 3rd term

                    index = index + nElem;
                end
            end
        end
     
        % Assemble the matrix        
      
        % Discontinuity capturing stabilization terms for Navier-Stokes equations
        % Part 4: 
        %

        
        for p = 1:nQuad  

            index = 0;
            for i = 1:nen
                for j = 1:nen
                    


                    index = index + nElem;
                end
            end
        end

        
        %% Allen-Cahn equation
        % Interpolate for alpha values for velocity
        


        
        % % Lagrange multiplier term and mass
        % sB1 = [];
        % sB2 = [];
         sB3 = [];
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
             localvapFracPrev = sum(repmat(N(p,:),nElem,1).*vapFracPrev,2);
             Bij3 = gW(p)*localvapFracPrev.*volume ;

             sB3(1:nElem,p) = Bij3 ;
         end

         sB3 = sum(sum(sB3,2));

         MassvapFrac = sB3 ;
	clear sB3 Bij3
        
        % Galerkin terms for Allen-Cahn equation

        for p = 1:nQuad  			

            index = 0;
            for i = 1:nen
                for j = 1:nen
                    % Galerkin inertia term
                    
                    % Galerkin convection term
                    
                    % Galerkin diffusion term
                    Aij_4 = gW(p)*(DNDx(:,i).*DNDx(:,j));
                    Aij_5 = gW(p)*(DNDy(:,i).*DNDy(:,j));
                    % Aij_4 = Aij_4.*volume;
                    % Aij_5 = Aij_5.*volume;
                    Aij_4 = 0.0;
                    Aij_5 = 0.0;					
                    sA4(index+1:index+nElem,p) = Aij_4;
                    sA5(index+1:index+nElem,p) = Aij_5;
                    
                    % Galerkin reaction term 
                    
                    % Galerkin source term matrices 

                    index = index + nElem;
                end
            end
        end
        
        % Assemble the matrix
        

        

        

        
        
        
       
        % Linear system


        
        
        % Update the solution
		
        % Check convergence criteria

        
        clear linSol1 linSol2
        clear Sol.uNew Sol.uDotNew Sol.uNewPrev Sol.uDotNewPrev Sol.vapFracAlpha vapFrac
        clear uxDot uyDot uzDot pres locGradU ResU gGradV
        
%         fsT = sum(etaRvapFrac,2) ;
%         xx = sum(xx,2) ;
% 
        
%         Sol.node(:,3) = [];

        theta = 0.5 ;
        theta_c = 0.05 ;
%         etaR = computeEtaR(Sol,fsT,solver.epsilon^2) ;
% 
        activeNode = find(Sol.type~=0); 
        nonActive = setdiff(size(crd,1),activeNode);
%         fprintf('EtaR: %e, ',sqrt(sum(etaR))) ;
        fprintf('nElem: %d, ',size(Sol.elem,1)) ;
        fprintf('ndofs: %d\n',size(activeNode,1)) ;
%         crd = [Sol.node zeros(size(Sol.node,1),1)] ;
%         Sol.node = [Sol.node zeros(size(Sol.node,1),1)] ;
        
        flagCoarsen = 0 ;
        flagRefine = 0 ;
        flagBreak = 0 ;
        if ( normIndicator1 <= solver.nLTol1 && normIndicator2 <= solver.nLTol1 ) %(nLiter == solver.nLIterMax))
            flagRefine = 0 ;
            flagCoarsen = 0 ;
            flagBreak = 1 ;
        end
     
        if (flagBreak == 1)
            break;
        end
        
    end
    fprintf('\n');
    % Copy current to previous
    
        for i=1:nen
        vapFraccheckphase(:,i)=Sol.vapFrac(cnn(:,i),1);
        uvali(:,i)=Sol.u(cnn(:,i),2);
    end
    vapFraclocal=vapFraccheckphase*N(:,nen);
    ind4=find(vapFraclocal< 0);
    for p=1:nQuad
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
            
            bubbleEv(:,p) = gW(p).*volume(ind4,1);
            uint(:,p) = gW(p)*uvali(ind4,:)*N(:,p).*volume(ind4,1);
            cymass(:,p)=gW(p)*yyf(ind4,:)*N(:,p).*volume(ind4,1);
    end
    bubblevolume(timeStep)=sum(sum(bubbleEv,2));
    uaverage(timeStep) = sum(sum(uint,2))/bubblevolume(timeStep);
    comy(timeStep)=sum(sum(cymass,2))/bubblevolume(timeStep);
%     massrec(timeStep) = MassvapFrac;
%     save('uaverage.mat','uaverage','comy','bubblevolume','massrec');
    clear bubbleEv uint cymass ind4 vapFracchecphase uvali
%     figure(4);
    % GravapFracc representation 
    toc
end

fclose(fileOthdId);
toc
fprintf('\n');
end
