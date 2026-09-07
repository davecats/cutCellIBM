/*---------------------------------------------------------------------------*\
    cutCellIBM: concrete immersed bodies

    cylinderBody   infinite circular cylinder: centre, axis, radius
    sphereBody     sphere: centre, radius
    boxBody        axis-aligned box: min, max
    planeBody      half space: point, normal (normal points into the fluid)
    searchableSurfaceBody
                   any OpenFOAM searchableSurface with inside/outside
                   information (triSurfaceMesh from an STL, searchableSphere,
                   searchableCylinder, ...). The level set is the signed
                   distance to the surface.
    alphaFieldBody the fluid volume fraction is read from a volScalarField;
                   the face fractions are reconstructed from the 0.5
                   iso-surface of the point-interpolated field.
\*---------------------------------------------------------------------------*/

#ifndef ibmBodies_H
#define ibmBodies_H

#include "ibmBody.H"
#include "volFields.H"

namespace Foam
{

class searchableSurface;

class cylinderBody : public ibmBody
{
    point centre_;
    vector axis_;
    scalar radius_;
public:
    TypeName("cylinder");
    cylinderBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
};


class sphereBody : public ibmBody
{
    point centre_;
    scalar radius_;
public:
    TypeName("sphere");
    sphereBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
};


class boxBody : public ibmBody
{
    point min_;
    point max_;
public:
    TypeName("box");
    boxBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
};


class planeBody : public ibmBody
{
    point point_;
    vector normal_;
public:
    TypeName("plane");
    planeBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
};


class searchableSurfaceBody : public ibmBody
{
    autoPtr<searchableSurface> surfPtr_;
public:
    TypeName("searchableSurface");
    searchableSurfaceBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
};


class alphaFieldBody : public ibmBody
{
    word fieldName_;
    bool solidFraction_;
public:
    TypeName("alphaField");
    alphaFieldBody(const word&, const fvMesh&, const dictionary&);
    virtual bool hasLevelSet() const { return false; }

    //- Fluid volume fraction read from the field (clipped to [0,1])
    tmp<volScalarField> alpha() const;
};

} // End namespace Foam

#endif
