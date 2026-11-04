# hexRef: 2-D and 3-D hexahedral refinement engine (third-party code)

The files in this directory are taken unchanged from **multiDimAMR**
(Henning Scheufler, <https://github.com/HenningScheufler/multiDimAMR>,
branch v2012, commit 346b4aab0b7b6b2b5fa45e048991f65890ba3e90), as patched
for OpenFOAM-v2106 in the ISTM `ibm_secondOrder` project
(`code/multiDimAMR/code/patches/ThirdParty`: `normal` -> `areaNormal`,
`Hash<labelPair>`, `autoPtr` validity tests).

| file | content |
|---|---|
| `hexRef.C/H` | common base: cell/point levels, 2:1 balance, unrefinement bookkeeping |
| `hexRef3D/` | 3-D split of a hexahedron into 8 (the OpenFOAM `hexRef8` algorithm) |
| `hexRef4/` | 2-D split into 4 for meshes one cell thick (Baniabedalruhman 2015, Ouedraogo 2017) |
| `hexRefData.C/H`, `refinementHistoryMultiDim.C/H` | refinement history for unrefinement and restarts |
| `hexRefNew.C` | selection of hexRef3D/hexRef4 from the mesh dimensions |

The axisymmetric variant (`hexRef4Axi`), the load balancing and the adaptation
criteria of multiDimAMR are not included; the refinement criteria of the
cut-cell IBM are in `../ibmRefineFvMesh`.

The load balancing of `../ibmRefineFvMesh` copies no file from multiDimAMR. It
uses `hexRef::distribute` and `refinementHistoryMultiDim::distribute` from this
directory, OpenFOAM's `fvMeshDistribute` (unmodified; multiDimAMR's
`fvMeshDistributeAddPatch` is not needed) and OpenFOAM's decomposition methods.
Two ideas are taken from multiDimAMR's `dynamicMultiDimRefineBalancedFvMesh`:
decompose the graph of the refinement trees (cells grouped by their top
ancestor in the history, weighted by their number of cells) instead of adding a
decomposition constraint (`refinementHistoryMultiDimConstraint`), and evaluate
the coupled patches of the volume fields after the redistribution.

Reference: D. Rettenmaier et al., Load balanced 2D and 3D adaptive mesh
refinement in OpenFOAM, SoftwareX 10 (2019) 100317.

## License

The code is derived from OpenFOAM and is distributed, like OpenFOAM and
multiDimAMR, under the GNU General Public License v3 or later; the original
copyright headers are kept in every file. The project as a whole is distributed under the
same license (`COPYING` at the top of the repository).
