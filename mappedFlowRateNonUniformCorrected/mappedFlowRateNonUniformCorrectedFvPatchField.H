/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::mappedFlowRateNonUniformCorrectedFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    Variant of mappedFlowRateCorrected that enforces a target water volumetric
    flow rate by a *multiplicative* rescaling of the normal velocity component
    instead of an additive uniform shift.

    This preserves the donor's normal-component profile shape (and therefore
    the relative size of normal fluctuations vs. mean) — the entire normal
    field is rescaled by a single scalar gamma. Tangential donor structure is
    passed through unchanged, exactly as in mappedFlowRateCorrected.

    Algorithm in updateCoeffs() each step:
      1) map donor U onto the inlet
      2) compute current water flux  Qw = sum( alpha * (U . n) * Sf )
      3) compute  gamma = -Qtarget / Qw   (outward normals -> minus sign)
      4) clip gamma to [gammaMin, gammaMax]
      5) split U into normal/tangential, write back
              U = Ut + gammaEff * Un * n,
         where  gammaEff = (1 - a) + a*gamma
         so cells with alpha = 0 receive the donor field untouched.

    Compared to the additive mappedFlowRateCorrected:
      - mean and fluctuations on the normal component are rescaled
        by the same factor gamma (ratios preserved)
      - no uniform offset is applied, so faces that already had small Un are
        not jolted by a large dUn when QwRec drifts
      - dry cells (alpha < alphaMin) keep the donor field

Usage
    \verbatim
    inlet
    {
        type                mappedFlowRateNonUniformCorrected;

        field               U;
        setAverage          no;
        average             uniform (0 0 0);
        interpolationScheme cell;

        Qtarget             0.089;       // target water volumetric flow rate
        alphaMin            1e-3;        // wet-cell threshold
        gammaMin            0.5;         // lower clip on rescale factor
        gammaMax            2.0;         // upper clip on rescale factor

        sampleMode          nearestCell;
        sampleRegion        region0;
        samplePatch         none;
        offsetMode          uniform;
        offset              (0 2 0);

        value               uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    mappedFlowRateNonUniformCorrectedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFlowRateNonUniformCorrectedFvPatchField_H
#define mappedFlowRateNonUniformCorrectedFvPatchField_H

#include "mappedPatchBase.H"
#include "mappedFlowRateCorrectedPatchFieldBase.H"
#include "fixedValueFvPatchFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class mappedFlowRateNonUniformCorrectedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedFlowRateNonUniformCorrectedFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public mappedPatchBase,
    public mappedFlowRateCorrectedPatchFieldBase<Type>
{
    // Private Data

        //- Target water volumetric flow rate through the patch
        scalar Qtarget_;

        //- Below this alpha a face cell is treated as dry
        scalar alphaMin_;

        //- Lower clip on the multiplicative rescale factor gamma
        scalar gammaMin_;

        //- Upper clip on the multiplicative rescale factor gamma
        scalar gammaMax_;

public:

    //- Runtime type information
    TypeName("mappedFlowRateNonUniformCorrected");


    // Constructors

        //- Construct from patch and internal field
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct from patch, internal field and distance for normal type
        //  sampling
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,

            // mappedPatchBase
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance,

            // My settings
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct by mapping given field onto a new patch
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const mappedFlowRateNonUniformCorrectedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const mappedFlowRateNonUniformCorrectedFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        mappedFlowRateNonUniformCorrectedFvPatchField
        (
            const mappedFlowRateNonUniformCorrectedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedFlowRateNonUniformCorrectedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
