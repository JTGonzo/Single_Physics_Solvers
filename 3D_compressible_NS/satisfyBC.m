%% satisfyBC
%  This function imposes the fluid domain boundary conditions
%  for the multiphase solution variables
%
function [Sol] = satisfyBC(Sol, fluid, crd, bc, Flag3D)

%% Identify the fboundary and free nodes
        zerotypenodes = find(Sol.type == 0);
		
        bc.left.nodes = find(Sol.node(:,1)< -0.001+1e-8) ;
        bc.left.nodes = setdiff(bc.left.nodes, zerotypenodes)' ;
        
        bc.right.nodes = find(Sol.node(:,1)> 0.001-1e-8) ;
        bc.right.nodes = setdiff(bc.right.nodes, zerotypenodes)' ;
        
        bc.top.nodes = find(Sol.node(:,2)> 0.001-1e-8) ;
        bc.top.nodes = setdiff(bc.top.nodes, zerotypenodes)' ;
        
        bc.bottom.nodes = find(Sol.node(:,2)< -0.001+1e-8) ;
        bc.bottom.nodes = setdiff(bc.bottom.nodes, zerotypenodes)' ;
        
		if(Flag3D)
			bc.front.nodes = find(Sol.node(:,3)> 0.001-1e-8) ;
			bc.front.nodes = setdiff(bc.front.nodes, zerotypenodes)' ;
			bc.back.nodes = find(Sol.node(:,3)< -0.001+1e-8) ;
			bc.back.nodes = setdiff(bc.back.nodes, zerotypenodes)' ;
		end
		
        % bc.left.nodes = find((sqrt((Sol.node(:,1).^2) + (Sol.node(:,2).^2)) > 0.5*sqrt(2) - 1e-5) & (Sol.node(:, 2) <= 0.5) & (Sol.node(:, 2) >= -0.5) & (Sol.node(:, 1) < 0)) ;
        % bc.left.nodes = setdiff(bc.left.nodes, zerotypenodes)' ;
        % bc.right.nodes = find((sqrt((Sol.node(:,1).^2) + (Sol.node(:,2).^2)) > 0.5*sqrt(2) - 1e-5) & (Sol.node(:, 2) <= 0.5) & (Sol.node(:, 2) >= -0.5) & (Sol.node(:, 1) > 0)) ;
        % bc.right.nodes = setdiff(bc.right.nodes, zerotypenodes)' ;
        % bc.top.nodes = find((sqrt((Sol.node(:,1).^2) + (Sol.node(:,2).^2)) > 0.5*sqrt(2) - 1e-5) & (Sol.node(:, 1) <= 0.5) & (Sol.node(:, 1) >= -0.5) & (Sol.node(:, 2) > 0)) ;
        % bc.top.nodes = setdiff(bc.top.nodes, zerotypenodes)' ;
        % bc.bottom.nodes = find((sqrt((Sol.node(:,1).^2) + (Sol.node(:,2).^2)) > 0.5*sqrt(2) - 1e-5) & (Sol.node(:, 1) <= 0.5) & (Sol.node(:, 1) >= -0.5) & (Sol.node(:, 2) < 0)) ;
        % bc.bottom.nodes = setdiff(bc.bottom.nodes, zerotypenodes)' ;	
				
%% Satisfy Dirichlet boundary condition		
% at left boundary nodes
        if (bc.left.type == 1)
            dirichletNodes = bc.left.nodes' ;
            if (Flag3D == 1)
				Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.left.var(3),1) = bc.left.value(3).*ones(size(bc.left.nodes',1),1) ;                
            else
				Sol.u(dirichletNodes,bc.left.var(1),1) = bc.left.value(1).*ones(size(bc.left.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.left.var(2),1) = bc.left.value(2).*ones(size(bc.left.nodes',1),1) ;
            end
        elseif (bc.left.type == 3)
            dirichletNodes = bc.left.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.left.nodes',1),1) ;
        end       
        dirichletNodes = bc.left.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.left.nodes',1),1) ;

% at right boundary nodes        
        if (bc.right.type == 1)
            dirichletNodes = bc.right.nodes' ;
            if (Flag3D == 1)
				Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.right.var(3),1) = bc.right.value(3).*ones(size(bc.right.nodes',1),1) ;
                Sol.vapFrac(dirichletNodes,1) = ones(size(bc.right.nodes',1),1) ;
            else
				Sol.u(dirichletNodes,bc.right.var(1),1) = bc.right.value(1).*ones(size(bc.right.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.right.var(2),1) = bc.right.value(2).*ones(size(bc.right.nodes',1),1) ;
            end
        elseif (bc.right.type == 3)
            dirichletNodes = bc.right.nodes' ;
            Sol.u(dirichletNodes,1,1) = zeros(size(bc.right.nodes',1),1) ;
        end 
        dirichletNodes = bc.right.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.right.nodes',1),1) ;

% at top boundary nodes        
        if (bc.top.type == 1)
            dirichletNodes = bc.top.nodes' ;
            if (Flag3D == 1)
				Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.top.var(3),1) = bc.top.value(3).*ones(size(bc.top.nodes',1),1) ;
                Sol.vapFrac(dirichletNodes,1) = ones(size(bc.top.nodes',1),1) ;
            else
				Sol.u(dirichletNodes,bc.top.var(1),1) = bc.top.value(1).*ones(size(bc.top.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.top.var(2),1) = bc.top.value(2).*ones(size(bc.top.nodes',1),1) ;
            end
        elseif (bc.top.type == 3)
            dirichletNodes = bc.top.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.top.nodes',1),1) ;
        end
        dirichletNodes = bc.top.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.top.nodes',1),1) ;        

% at bottom boundary nodes          
        if (bc.bottom.type == 1)
            dirichletNodes = bc.bottom.nodes' ;
            if (Flag3D == 1)
				Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.bottom.var(3),1) = bc.bottom.value(3).*ones(size(bc.bottom.nodes',1),1) ;
                Sol.vapFrac(dirichletNodes,1) = ones(size(bc.bottom.nodes',1),1) ; %fluid.vapFrac0*ones(size(bc.bottom.nodes',1),1) ;
            else
				Sol.u(dirichletNodes,bc.bottom.var(1),1) = bc.bottom.value(1).*ones(size(bc.bottom.nodes',1),1) ;
				Sol.u(dirichletNodes,bc.bottom.var(2),1) = bc.bottom.value(2).*ones(size(bc.bottom.nodes',1),1) ;
				Sol.vapFrac(dirichletNodes,1) = fluid.vapFrac0*ones(size(bc.bottom.nodes',1),1) ;
            end
        elseif (bc.bottom.type == 3)
            dirichletNodes = bc.bottom.nodes' ;
            Sol.u(dirichletNodes,2,1) = zeros(size(bc.bottom.nodes',1),1) ;
        end       
        dirichletNodes = bc.bottom.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.bottom.nodes',1),1) ;        

% at front boundary nodes         
        if (bc.front.type == 1)
            dirichletNodes = bc.front.nodes' ;
			Sol.u(dirichletNodes,bc.front.var(1),1) = bc.front.value(1).*ones(size(bc.front.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.front.var(2),1) = bc.front.value(2).*ones(size(bc.front.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.front.var(3),1) = bc.front.value(3).*ones(size(bc.front.nodes',1),1) ;
            Sol.vapFrac(dirichletNodes,1) = ones(size(bc.front.nodes',1),1) ; %fluid.vapFrac0*ones(size(bc.front.nodes',1),1) ;
        elseif (bc.front.type == 3)
            dirichletNodes = bc.front.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.front.nodes',1),1) ;
        end        
        dirichletNodes = bc.front.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.front.nodes',1),1) ;        

% at back boundary nodes            
        if (bc.back.type == 1)
            dirichletNodes = bc.back.nodes' ;
			Sol.u(dirichletNodes,bc.back.var(1),1) = bc.back.value(1).*ones(size(bc.back.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.back.var(2),1) = bc.back.value(2).*ones(size(bc.back.nodes',1),1) ;
            Sol.u(dirichletNodes,bc.back.var(3),1) = bc.back.value(3).*ones(size(bc.back.nodes',1),1) ;
            Sol.vapFrac(dirichletNodes,1) = ones(size(bc.back.nodes',1),1) ; %fluid.vapFrac0*ones(size(bc.back.nodes',1),1) ;		
        elseif (bc.back.type == 3)
            dirichletNodes = bc.back.nodes' ;
            Sol.u(dirichletNodes,3,1) = zeros(size(bc.back.nodes',1),1) ;
        end       
        dirichletNodes = bc.back.nodes' ;
        Sol.vapFrac(dirichletNodes,1) = (1-fluid.alphaNuc)*ones(size(bc.back.nodes',1),1) ;		
				
%% 2D dirichlet condition in 3D mesh 
        % if (Flag3D == 1)
			% Sol.u(:,3,1) = zeros(size(crd,1),1) ;
        % end
       
	   % Sol.p(1,1) = fluid.pInf ;
       % Sol.p(bc.top.nodes',1) = fluid.pInf ;
	   % Sol.p(bc.bottom.nodes',1) = fluid.pInf ;
       % Sol.p(bc.right.nodes',1) = fluid.pInf ;
	   % Sol.p(bc.left.nodes',1) = fluid.pInf ;	

       % Sol.vapFrac(bc.top.nodes',1) = 1.0 ;
	   % Sol.vapFrac(bc.bottom.nodes',1) = 1.0 ;
       % Sol.vapFrac(bc.right.nodes',1) = 1.0 ;
	   % Sol.vapFrac(bc.left.nodes',1) = 1.0 ;	  		
end