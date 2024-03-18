function [A_in, F_in, u_D] =  CSM_ApplyBC(A, F, FE_SPACE, MESH, DATA, t, zero_Dirichlet)
if nargin < 6
    t = [];
end

if isempty(A)
    A = sparse(MESH.numNodes*MESH.dim, MESH.numNodes*MESH.dim);
end

if isempty(F)
    F = sparse(MESH.numNodes*MESH.dim, 1);
end

if nargin < 7
    zero_Dirichlet = 0;
end

param = DATA.param;

u_D = [];

        %% Dirichlet condition
        for k = 1 : 2
            x  = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
            y  = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
            u_Dirichlet{k} = DATA.bcDir{k}(x,y,t,param);               
            u_D = [u_D; u_Dirichlet{k}'];
        end
        
    u_D  = u_D * (1 - zero_Dirichlet);

    F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;
    
    A_in = A(MESH.internal_dof,MESH.internal_dof);

end
