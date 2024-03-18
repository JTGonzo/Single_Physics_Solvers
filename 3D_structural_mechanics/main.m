%clc
clear all; close all; clc;
[~,~,~] = mkdir('Figures');
addpath('C_Files/')
addpath('Tools/')

%% Load P1 mesh
load('beam.mat')

dim      =  3;
fem      =  'P2';

t      = [];
param = [];
vtk_filename = 'Figures/Sol_'; % []; % 
sol_history = true;
data_file = 'datafile';

%% Read problem parameters and BCs from data_file
eval(data_file);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end

DATA.param = param;

quad_order = 5;

[ MESH, FE_SPACE ] = build( dim, elements, vertices, boundaries, fem, quad_order, DATA, 'CSM' );


%% Gather Time Setting
t0  = DATA.time.t0;
dt  = DATA.time.dt;
tf  = DATA.time.tf;
t   = DATA.time.t0;
k_t = 0;

TimeAdvance = GeneralizedAlpha_TimeAdvance( DATA.time.beta, DATA.time.gamma, DATA.time.alpha_m, DATA.time.alpha_f, dt );

u0  = [];
du0 = [];

for k = 1 : FE_SPACE.numComponents
    
    u0  = [u0; DATA.u0{k}(  MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];
    du0 = [du0; DATA.du0{k}( MESH.nodes(1,:), MESH.nodes(2,:), MESH.nodes(3,:), t0, param )'];

end

u = u0;
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, 0);
end

Coef_Mass = TimeAdvance.MassCoefficient( );

%% Generate Domain Decomposition (if required)
SolidModel = CSM_Assembler( MESH, DATA, FE_SPACE );

%% Assemble mass matrix
M    =  SolidModel.compute_mass()* DATA.Density;

% Assemble Robin BC (if it's the case)
A_robin = SolidModel.assemble_ElasticRobinBC();

%% Initial Acceleration
F_ext_0  = SolidModel.compute_volumetric_forces(t0);

F_in_0   =  SolidModel.compute_internal_forces( u0 ) + A_robin * u0;

d2u0 = M \ (F_ext_0 - F_in_0);

TimeAdvance.Initialize( u0, du0, d2u0 );

U_n = u0;
step = 0;

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );

    % Newton Method
    tol        = DATA.NonLinearSolver.tol;
    resRelNorm = tol + 1;
    incrNorm   = tol + 1;
    maxIter    = DATA.NonLinearSolver.maxit;
    k          = 1;

    [~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE, MESH, DATA, t);
    dU             = zeros(MESH.numNodes*MESH.dim,1);
    U_k            = u(:,end);
    U_k(MESH.Dirichlet_dof) = u_D;
    
    Csi = TimeAdvance.RhsContribute( );

    % Assemble matrix and right-hand side
    F_ext      = SolidModel.compute_volumetric_forces( (1 - TimeAdvance.M_alpha_f) * t + TimeAdvance.M_alpha_f * (t-dt) );
    
    F_in      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
    
    Residual  = Coef_Mass * M * U_k + F_in - F_ext - M * Csi ...
                + A_robin * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n);
            
    dF_in     = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );      
            
    Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in + A_robin * (1 - TimeAdvance.M_alpha_f);
    
    % Apply boundary conditions
    [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1);
    
    res0Norm = norm(b);
    
    fprintf('\n============ Start Newton Iterations ============\n\n');
    while (k <= maxIter && incrNorm > tol && resRelNorm > tol)
        
        % Solve
        dU(MESH.internal_dof) = A \ b;
        
        % update solution
        U_k        = U_k + dU;
        incrNorm   = norm(dU)/norm(U_k);
        
                % Assemble matrix and right-hand side
        F_in      = SolidModel.compute_internal_forces( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        dF_in     = SolidModel.compute_jacobian( (1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n );
        
        Residual  = Coef_Mass * M * U_k + F_in - F_ext - M * Csi ...
                    + A_robin * ((1 - TimeAdvance.M_alpha_f) * U_k + TimeAdvance.M_alpha_f * U_n);
            
        Jacobian  = Coef_Mass * M + (1 - TimeAdvance.M_alpha_f) * dF_in + A_robin * (1 - TimeAdvance.M_alpha_f);

        % Apply boundary conditions
        [A, b]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE, MESH, DATA, t, 1);
        
        resRelNorm = norm(b) / res0Norm;
        
        fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;
        
    end
    
    if sol_history
        u = [u U_k];
    else
        u = U_k;
    end
    
%% Export to VTK   
    if (mod(step,10)==0)
        if ~isempty(vtk_filename)
            CSM_export_solution(MESH.dim, U_k, MESH.vertices, MESH.elements, MESH.numNodes, vtk_filename, k_t);
        end
    end
    
    TimeAdvance.Update( U_k );
    
    U_n = U_k;
    
    iter_time = toc(iter_time);
    fprintf('\n-------------- Iteration time: %3.2f s -----------------',iter_time);
    
end