/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

#include "Tools.h"
#include "LinearElasticMaterial.h"
#include "StVenantKirchhoffMaterial.h"

#ifdef _OPENMP
#include <omp.h>
#else
#warning "OpenMP not enabled. Compile with mex CSM_assembler_C_omp.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp""
#endif

/*************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    char *Material_Model = mxArrayToString(prhs[1]);
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
            
    if (strcmp(Material_Model, "Linear_forces")==0)
    {
            LinearElasticMaterial_forces(plhs, prhs);
    }
    
    
    if (strcmp(Material_Model, "Linear_jacobian")==0)
    {
        if (dim == 2)
        {
            LinearElasticMaterial_jacobianFast2D(plhs, prhs);
        }
        
    }

    if (strcmp(Material_Model, "StVenantKirchhoff_forces")==0)
    {
            StVenantKirchhoffMaterial_forces(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "StVenantKirchhoff_jacobianSlow")==0)
    {
            StVenantKirchhoffMaterial_jacobian(plhs, prhs);
    }
    
    if (strcmp(Material_Model, "StVenantKirchhoff_jacobian")==0)
    {
        if (dim == 2)
        {
            StVenantKirchhoffMaterial_jacobianFast2D(plhs, prhs);
        }        
    }
    
    if (strcmp(Material_Model, "StVenantKirchhoff_stress")==0)
    {
            StVenantKirchhoffMaterial_stress(plhs, prhs);
    }
    
    mxFree(Material_Model);

}
/*************************************************************************/

