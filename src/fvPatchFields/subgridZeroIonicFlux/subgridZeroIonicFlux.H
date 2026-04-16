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

#ifndef subgridZeroIonicFlux_H
#define subgridZeroIonicFlux_H


#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H" 
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
    namespace constants
    {
        namespace chemical
        {
            namespace electric
            {
                const dimensionedScalar e  // [C] = [kg0⋅m0⋅s1⋅K0⋅mol0⋅A1⋅cd0] 
                ("elementaryCharge",
                dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),
                1.6021766208e-19);
                const dimensionedScalar eps_0  // [F/m] = [kg-1⋅m-3⋅s4⋅K0⋅mol0⋅A2⋅cd0]
                ( "vacuumPermittivity", 
                dimensionSet(-1, -3, 4, 0, 0, 2, 0),
                8.8541878176e-12 ) ;
                const dimensionedScalar A // [ 1/mol ] = [kg0⋅m0⋅s0⋅K0⋅mol-1⋅A0⋅cd0] 
                ( "AvogradoNumber", 
                dimensionSet(0, 0, 0, 0, -1, 0, 0),
                6.022140857e+23);
                const dimensionedScalar F // [C/mol] = [kg0⋅m0⋅s0⋅K0⋅mol-1⋅A1⋅cd0]
                ( "FaradayConstant", 
                dimensionSet(0, 0, 0, 0,-1, 1, 0),
                96485.33212331001); 
                const dimensionedScalar R // [J/molK] = [kg1⋅m2⋅s-2⋅K-1⋅mol-1⋅A0⋅cd0]
                ( "gasConstant", 
                dimensionSet(1, 2, -2, -1,-1, 0, 0),
                8.31446261815324); 
                const dimensionedScalar k_B  // [J/K] = [kg1⋅m2⋅s-2⋅K-1⋅mol0⋅A0⋅cd0] 
                ( "BoltzmannConstant", 
                dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),
                1.38064852e-23); 
            }
        }
    }
} // End namespace Foam
 
namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class subgridZeroIonicFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class subgridZeroIonicFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data
        //- mesh
        const fvMesh& mesh; 
        //- Ionic valenceNumber
        scalar zi_; // const not possible: see constructor  
        //- Reference potential
        scalar  Psi0_ ;
        //- Species information
        const word specieName_; 
        //- Field type: c is derived from Y
        const char fieldTypeName_;  

        const IOdictionary electricDict_;

public:

    //- Runtime type information
    TypeName("subgridZeroIonicFlux");

    // Constructors

        //- Construct from patch and internal field
        subgridZeroIonicFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subgridZeroIonicFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given subgridZeroIonicFluxFvPatchScalarField onto
        //  a new patch
        subgridZeroIonicFluxFvPatchScalarField
        (
            const subgridZeroIonicFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subgridZeroIonicFluxFvPatchScalarField
        (
            const subgridZeroIonicFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new subgridZeroIonicFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subgridZeroIonicFluxFvPatchScalarField
        (
            const subgridZeroIonicFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new subgridZeroIonicFluxFvPatchScalarField(*this, iF)
            );
        }

        //- Update the coefficients associated with the patch field
            virtual void updateCoeffs()
        {
            if (updated())
            {
                return;       
            }
            
            const volScalarField& Psi_E_field = db().lookupObject<volScalarField>("Psi_E"); 

                
            scalarField Epatch( Psi_E_field.boundaryField()[patch().index()].snGrad() );
             
                gradient() = -(*this ) 
                    * zi_ * Foam::constants::chemical::electric::F.value() * Epatch / (Foam::constants::chemical::electric::R.value()* 300)/1000;      
             
            fixedGradientFvPatchScalarField::updateCoeffs();
        };
         
        //- Evaluate the patch field
        //      This is equivalent to c_i,f 
        virtual void evaluate(const Pstream::commsTypes commsType=Pstream::blocking)
        {
            const fvPatchField<scalar>& Psi_E_field =patch().lookupPatchField<volScalarField, scalar>("Psi_E") ; 

            scalarField deltaPsi = Psi_E_field - Psi_E_field.patchInternalField();

                    Field<Foam::scalar>::operator=
                    (
                        ( this->patchInternalField() /* c_i,P  */)
                        * Foam::exp( -Foam::constants::chemical::electric::F.value() * zi_* deltaPsi/ (Foam::constants::chemical::electric::R.value() * 300) ) /1000 
                    );
            
            // This will set update_ = false and force updateCoeffs(). Gradient
            // should be updated always we call operator= in order to ensure
            // that both use the same deltaPsi.  
            fvPatchField<Foam::scalar>::evaluate();
            
            if (!this->updated())
            {
                this->updateCoeffs();
            }

            fvPatchField<Foam::scalar>::evaluate();
        };

        // Render the BC fully-explicit
         
        Foam::tmp<Field<Foam::scalar> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const
        {
            return tmp<Field<Foam::scalar> >(new Field<Foam::scalar>(this->size(), pTraits<Foam::scalar>::zero));
        }
         
        Foam::tmp<Field<Foam::scalar> > 
        gradientInternalCoeffs() const
        {
            return tmp<Field<Foam::scalar>>
            (
                new Field<Foam::scalar>(this->size(), Zero)
            );
        }
         
        Foam::tmp<Field<Foam::scalar> > 
        gradientBoundaryCoeffs() const
        {
            return gradient();    
        }
                
        //- Write
        virtual void write(Ostream& os) const
        {
            fixedGradientFvPatchScalarField::write(os);
            writeEntry("value", os); 
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 
