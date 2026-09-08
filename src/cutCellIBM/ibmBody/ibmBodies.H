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
    alphaFieldBody the fluid volume fraction is read from a volScalarField
                   (registered on the mesh, or from the time directory); the
                   face fractions are reconstructed from the 0.5 iso-surface
                   of the point-interpolated field. With "moving true" the
                   field is re-fetched at every time step and the body
                   velocity is "velocity (ux uy uz)" or "velocityField Ub".
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
    virtual point centre() const { return centre_; }
};


class sphereBody : public ibmBody
{
    point centre_;
    scalar radius_;
public:
    TypeName("sphere");
    sphereBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
    virtual point centre() const { return centre_; }
};


class boxBody : public ibmBody
{
    point min_;
    point max_;
public:
    TypeName("box");
    boxBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
    virtual point centre() const { return 0.5*(min_ + max_); }
};


class planeBody : public ibmBody
{
    point point_;
    vector normal_;
public:
    TypeName("plane");
    planeBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
    virtual point centre() const { return point_; }
};


class searchableSurfaceBody : public ibmBody
{
    autoPtr<searchableSurface> surfPtr_;
public:
    TypeName("searchableSurface");
    searchableSurfaceBody(const word&, const fvMesh&, const dictionary&);
    virtual void levelSet(const pointField&, scalarField&) const;
    virtual point centre() const;
};


class alphaFieldBody : public ibmBody
{
    word fieldName_;
    bool solidFraction_;

    //- True if the field changes in time (re-fetched at every update)
    bool moving_;

    //- Uniform body velocity, or the name of a velocity field
    vector velocity_;
    word velocityFieldName_;

public:
    TypeName("alphaField");
    alphaFieldBody(const word&, const fvMesh&, const dictionary&);
    virtual bool hasLevelSet() const { return false; }
    virtual point centre() const;
    virtual bool moving() const { return moving_; }
    virtual void velocity(const pointField& pts, vectorField& U) const;

    //- Fluid volume fraction (clipped to [0,1]): looked up from the object
    //  registry if a field of that name is registered (e.g. solved for by
    //  the solver), otherwise read from the current time directory
    tmp<volScalarField> alpha() const;
};

} // End namespace Foam

#endif
