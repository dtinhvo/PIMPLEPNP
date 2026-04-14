/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "subgridZeroIonicFlux.H"

// constructions are copied from rheotool
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subgridZeroIonicFluxFvPatchScalarField::subgridZeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    mesh(iF.mesh()),
    electricDict_(
        IOobject
            (
                "electricProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ) 
    ), 
    zi_(0),
    Psi0_(0) ,
    specieName_(iF.name().substr(1)),
    fieldTypeName_(iF.name()[0])
{}


Foam::subgridZeroIonicFluxFvPatchScalarField::subgridZeroIonicFluxFvPatchScalarField
(
    const subgridZeroIonicFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    mesh(iF.mesh()),
    electricDict_(ptf.electricDict_), 
    zi_(ptf.zi_),
    Psi0_(ptf.Psi0_),
    specieName_(iF.name().substr(1)),
    fieldTypeName_(iF.name()[0])
{}


Foam::subgridZeroIonicFluxFvPatchScalarField::subgridZeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    mesh(iF.mesh()),
    electricDict_(
        IOobject
            (
                "electricProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ) 
    ), 
    zi_(0),
    Psi0_(0),
    specieName_(iF.name().substr(1)),
    fieldTypeName_(iF.name()[0])
{
    // re-read some vars, because improper dict init 
    zi_ = readScalar(mesh.lookupObject<IOdictionary>("electricProperties").lookup("z")); 
    fvPatchField<scalar>::operator=(patchInternalField());
}


Foam::subgridZeroIonicFluxFvPatchScalarField::subgridZeroIonicFluxFvPatchScalarField
(
    const subgridZeroIonicFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    specieName_(wbppsf.specieName_),
    fieldTypeName_(wbppsf.fieldTypeName_),
    mesh(wbppsf.mesh),
    electricDict_ (wbppsf.electricDict_ ),
    Psi0_(wbppsf.Psi0_),
    zi_(wbppsf.zi_)
{}


Foam::subgridZeroIonicFluxFvPatchScalarField::subgridZeroIonicFluxFvPatchScalarField
(
    const subgridZeroIonicFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    mesh(wbppsf.mesh),
    electricDict_ (wbppsf.electricDict_ ), 
    zi_(wbppsf.zi_),
    Psi0_(wbppsf.Psi0_),
    specieName_(iF.name().substr(1)),
    fieldTypeName_(iF.name()[0])
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        subgridZeroIonicFluxFvPatchScalarField
    );
}

// ************************************************************************* //
