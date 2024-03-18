%DATAFILE for linear elasticity problem

% Source term
%data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{1} = @(x, y, z, t, param)(-3*1 + 0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);
%data.force{3} = @(x, y, z, t, param)(-0.05*1 + 0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y.*z);
%data.bcNeu{3} = @(x, y, z, t, param)(0 + 0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0*x.*y.*z);
%data.bcPrex   = @(x, y, z, t, param)(-2*(t<10) + 0*x.*y.*z);

% BC flag
data.flag_dirichlet{1} = [1];
data.flag_neumann{1}   = [2 3];
data.flag_pressure{1}  = [ ];
data.flag_robin{1}     = [];
data.flag_clamp_points{1} = [];

data.flag_dirichlet{2} = [1];
data.flag_neumann{2}   = [2 3];
data.flag_pressure{2}  = [ ];
data.flag_robin{2}     = [];
data.flag_clamp_points{2} = [];

data.flag_dirichlet{3} = [1 3];
data.flag_neumann{3}   = [2];
data.flag_pressure{3}  = [ ];
data.flag_robin{3}     = [];
data.flag_clamp_points{3} = [];

%initial condition
data.u0{1} = @(x, y, z, t, param)(0.*x.*y);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y);

% material parameters
data.Material_Model = 'StVenantKirchhoff'; %'StVenantKirchhoff', 'Linear', 'NeoHookean'
data.model   = 'CSM';
data.Young   = 2100.5;
data.Poisson = 0.29;
data.Density = 10;
data.ElasticCoefRobin = 1e+6;
data.flag_dirichletNormal = [];

% NonLinear Solver
data.NonlinearSolver.first       = 'newton';
data.NonLinearSolver.tol               = 1e-7;
data.NonLinearSolver.maxit             = 25;

% Linear Solver
data.LinearSolver.type           = 'backslash'; 
data.LinearSolver.mumps_reordering  = 5;
    
% Preconditioner
data.Preconditioner.type         = 'None'; 

% Time options
data.time.t0         = 0;
data.time.dt         = 0.01;
data.time.tf         = 8;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;
