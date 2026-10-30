/*---------------------------------------------------------------------------*\
    cutCellIBM: turbulence models that see the immersed walls, for the
    incompressible solvers (ibmSimpleFoam, ibmPimpleFoam). Load with
    libs ("libcutCellIBMTurbulence.so"); in system/controlDict.
\*---------------------------------------------------------------------------*/

#include "turbulentTransportModels.H"
#include "ibmSpalartAllmaras.H"

makeRASModel(ibmSpalartAllmaras);
