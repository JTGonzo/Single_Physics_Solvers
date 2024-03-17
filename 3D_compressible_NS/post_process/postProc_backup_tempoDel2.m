% function [] = postProc(Sol, solver, bc, timeStep, fileOthdId, ...
% 						wrkDir, problemString, crd, elemType, Flag3D)
tempoCount = 1;
for tempo = 1:1:10
    
    saveFile = strcat('saveCaseData_',num2str(tempo),'.mat');
    load(saveFile);
      
      [gP, gW, N, Nx, Ny, Nz, nQuad] = defineQuadratureData(elemType, Flag3D);
 
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

%% Navier-Stokes equations
        xxf = zeros(size(cnn));
        yyf = zeros(size(cnn));
        zzf = zeros(size(cnn));

        for i=1:nen
            xxf(:,i) = crd(cnn(:,i),1);
            yyf(:,i) = crd(cnn(:,i),2);
            if (Flag3D == 1)
                zzf(:,i) = crd(cnn(:,i),3);
            end
            vapFrac(:,i) = Sol.vapFracAlpha(cnn(:,i),1) ;	
        end
%       sA = zeros(nen^2*nElem,nQuad);
% % % % % % % % % %     for p = 1:nQuad  
% % % % % % % % % %         if (Flag3D ~= 1)
% % % % % % % % % %            J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
% % % % % % % % % %                 yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
% % % % % % % % % %             if size(J,2)==1
% % % % % % % % % %                 J = J';
% % % % % % % % % %             end
% % % % % % % % % %             volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
% % % % % % % % % %         else
% % % % % % % % % %             J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)], xxf*[Nz(:,p)],...
% % % % % % % % % %                  yyf*[Nx(:,p)], yyf*[Ny(:,p)], yyf*[Nz(:,p)], ...
% % % % % % % % % %                  zzf*[Nx(:,p)], zzf*[Ny(:,p)], zzf*[Nz(:,p)]];
% % % % % % % % % %             if size(J,2)==1
% % % % % % % % % %                 J = J';
% % % % % % % % % %             end
% % % % % % % % % %             volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
% % % % % % % % % %                       J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
% % % % % % % % % %                       J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
% % % % % % % % % %         end
% % % % % % % % % % 
% % % % % % % % % %         volume = abs(volume);
% % % % % % % % % % 
% % % % % % % % % %         negJacobian = find(volume<0);
% % % % % % % % % %         if ~isempty(negJacobian)
% % % % % % % % % %            disp('Mesh deformed, Negative Jacobian');
% % % % % % % % % %            exit
% % % % % % % % % %         end
% % % % % % % % % % 
% % % % % % % % % %         localvapFrac  = sum(repmat(N(p,:),nElem,1).*vapFrac,2);	
% % % % % % % % % % 
% % % % % % % % % %         index = 0;
% % % % % % % % % %         for i = 1:nen
% % % % % % % % % %             for j = 1:nen
% % % % % % % % % %                 % Galerkin inertia term
% % % % % % % % % % %                 if(localvapFrac < 1.0)
% % % % % % % % % %                     Vij = gW(p)*localvapFrac.*volume ;
% % % % % % % % % %                     totalVol = gW(p)*volume;
% % % % % % % % % % %                 end    
% % % % % % % % % %                 index = index + nElem;
% % % % % % % % % % 
% % % % % % % % % %             end
% % % % % % % % % %         end
% % % % % % % % % %     end
    
        
    for i=1:nen
        vapFraccheckphase(:,i)=Sol.vapFrac(cnn(:,i),1);
%         uvali(:,i)=Sol.u(cnn(:,i),2);
    end
    
    vapFraclocal=vapFraccheckphase*N(:,nen);
    ind4=find(vapFraclocal< 0.95);
    ind5=find(vapFraclocal>= 0.95);
    
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
            
            bubbleEv(ind4,p) = gW(p)*volume(ind4,1);
            bubbleEv(ind5,p) = 0.0;
%             bubbleEv = gW(p)*(1-vapFraclocal).*volume;
                
%             uint(:,p) = gW(p)*uvali(ind4,:)*N(:,p).*volume(ind4,1);
%             cymass(:,p)=gW(p)*yyf(ind4,:)*N(:,p).*volume(ind4,1);
    end
    
    bubbleV=sum(sum(bubbleEv,2));
    bubblevolume(tempoCount) = (3/(4*pi)*bubbleV)^(1/3);
%     uaverage(timeStep) = sum(sum(uint,2))/bubblevolume(timeStep);
%     comy(timeStep)=sum(sum(cymass,2))/bubblevolume(timeStep);
%     massrec(timeStep) = MassPhi;
%     save('uaverage.mat','uaverage','comy','bubblevolume','massrec');

    % Summation of all quadrature data
%     vBub = sum(sum(Vij,2)); 
%     totalVol = sum(sum(totalVol,2));
    
%     figHandle = figure;
%     plot(timeStep*solver.dt, (3/(4*pi)*bubblevolume)^(1/3), '*');

      clear bubbleEv
%     ylim([950 1000]);
%     pause(0.1);

%      title(sprintf('Cavitation Fraction at time %9.8f',timeStep*solver.dt), 'FontSize', 14)
%       % Capture the plot as an image 
%       frame = getframe(gca); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if timeStep == 20 
%           imwrite(imind,cm,'dens1.gif','gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,'dens1.gif','gif','WriteMode','append'); 
%       end
    hold on;
    tempoCount = tempoCount + 1;
%     close;

end