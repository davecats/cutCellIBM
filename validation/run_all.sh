#!/bin/bash
# Reproduce the validation of the cut-cell IBM solvers against the notebook.
# Needs: OpenFOAM-v2106 sourced, the project built (./Allwmake), and the
# notebook's python environment (fvm/.venv) for the reference and comparison.
cd "${0%/*}" || exit 1
set -e
PY=${PY:-python3}   # needs numpy, scipy, matplotlib
source ../etc/bashrc
T=../tutorials

echo "== reference solutions (notebook numerics, ~40 s each)"
[ -f ref_nx81.npz ]          || $PY reference_notebook.py --nIter 600 --out ref_nx81.npz > ref_nx81.log
[ -f ref_nx81_freeSlip.npz ] || $PY reference_notebook.py --nIter 600 --bc freeSlip --out ref_nx81_freeSlip.npz > ref_nx81_freeSlip.log

run() {   # run <case> <solver> [mpi ranks]
    local c=$1 app=$2 np=${3:-1}
    echo "== $c ($app)"
    ( cd $T/$c && ./Allclean > /dev/null 2>&1; rm -rf processor* log.*
      blockMesh > log.blockMesh 2>&1
      postProcess -func writeCellCentres -time 0 > log.postProcess 2>&1
      ibmSetGeometry -writeFaceCentres > log.ibmSetGeometry 2>&1
      if [ $np -gt 1 ]; then
          decomposePar > log.decomposePar 2>&1
          mpirun -np $np $app -parallel > log.$app 2>&1
          reconstructPar -latestTime > log.reconstructPar 2>&1
      else
          $app > log.$app 2>&1
      fi )
}

run cylinderRe1                ibmSimpleFoam
run cylinderRe1_pimple         ibmPimpleFoam
run cylinderRe1_stl            ibmSimpleFoam
run cylinderRe1_alphaField     ibmSimpleFoam
run cylinderRe1_pimple_freeSlip ibmPimpleFoam
run twoCylinders               ibmSimpleFoam
run cylinderRe1_parallel       ibmSimpleFoam 4

echo; echo "== comparisons with the notebook"
$PY compare.py $T/cylinderRe1                 --fig fig_simple_vs_notebook.png
$PY compare.py $T/cylinderRe1_pimple          --fig fig_pimple_vs_notebook.png
$PY compare.py $T/cylinderRe1_stl
$PY compare.py $T/cylinderRe1_alphaField
$PY compare.py $T/cylinderRe1_parallel
$PY compare.py $T/cylinderRe1_pimple_freeSlip --ref ref_nx81_freeSlip.npz --fig fig_pimple_freeSlip_vs_notebook.png
echo; echo "== two bodies: forces must sum to gradP * fluid volume"
grep -E "IBM force|pressure gradient" $T/twoCylinders/log.ibmSimpleFoam | tail -3
grep "fluid volume" $T/twoCylinders/log.ibmSimpleFoam | head -1
