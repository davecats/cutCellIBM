#!/bin/bash
# Richardson study of the periodic-cylinder drag (Re_D = 1) for no-slip and
# free-slip walls, with and without the second-order wall treatment, on
# N = 41 81 161 [321]. Cases are built in $1 (default: runs) from
# tutorials/cylinderRe1; prints the drag per unit depth and the observed order.
cd "${0%/*}" || exit 1
W=${1:-runs}; NS=${NS:-"41 81 161"}
T=../tutorials
for so in false true; do for wall in noSlip freeSlip; do for N in $NS; do
    d=$W/rich_${wall}_so${so}_N$N; mkdir -p $d/0 $d/constant $d/system
    cp $T/cylinderRe1/0/U $T/cylinderRe1/0/p $d/0/
    cp $T/cylinderRe1/constant/{transportProperties,turbulenceProperties,fvOptions,ibmProperties} $d/constant/
    cp $T/cylinderRe1/system/{blockMeshDict,controlDict,fvSchemes,fvSolution} $d/system/
    sed -i "s/nx 81; ny 81;/nx $N; ny $N;/" $d/system/blockMeshDict
    sed -i "s/wallType  noSlip;.*/wallType  $wall;/; s/^writeFields .*/writeFields   false;\nsecondOrder   $so;/" $d/constant/ibmProperties
    sed -i "s/default Gauss linear corrected;/default Gauss linear cutCellCorrected;/; s/^snGradSchemes .*/snGradSchemes   { default cutCellCorrected; }/" $d/system/fvSchemes
    sed -i "s/^endTime .*/endTime         8000;/; s/^writeInterval .*/writeInterval   8000;/" $d/system/controlDict
    ( cd $d && blockMesh > log.blockMesh 2>&1 && ibmSimpleFoam > log.ibmSimpleFoam 2>&1 ) &
done; done; done
wait
for so in false true; do for wall in noSlip freeSlip; do
    echo -n "$wall secondOrder=$so:"
    prev=""; prevd=""
    for N in $NS; do
        d=$W/rich_${wall}_so${so}_N$N
        F=$(grep "IBM force" $d/log.ibmSimpleFoam | tail -1 | awk -v N=$N '{gsub(/[()]/,""); printf "%.5f", ($6+$10)/(10/N)}')
        c=$(grep -c "converged in" $d/log.ibmSimpleFoam)
        echo -n "  N=$N F=$F$( [ $c -eq 0 ] && echo '(not converged)')"
        if [ -n "$prevd" ]; then p=$(python3 -c "import math; print(f'{math.log(abs($prevd)/abs($prev-$F))/math.log(2):.2f}')"); echo -n " (p=$p)"; fi
        [ -n "$prev" ] && prevd=$(python3 -c "print($prev-$F)")
        prev=$F
    done; echo
done; done
