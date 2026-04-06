/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mappedFlowRateCorrectedFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this),
    Qtarget_(0.0),
    alphaMin_(1e-8)
{}


template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, dict, *this),
    Qtarget_(dict.lookupOrDefault<scalar>("Qtarget", 0.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 1e-8))
{
    if (!dict.found("value"))
    {
        this->evaluate();
    }
}


template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const mappedFlowRateCorrectedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_)
{}


template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,

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
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase
    (
        p.patch(),
        sampleRegion,
        sampleMode,
        samplePatch,
        distance
    ),
    mappedFlowRateCorrectedPatchFieldBase<Type>
    (
        *this,
        *this,
        fieldName,
        setAverage,
        average,
        interpolationScheme
    ),
    Qtarget_(0.0),
    alphaMin_(1e-8)
{}


template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const mappedFlowRateCorrectedFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_)
{}


template<class Type>
Foam::mappedFlowRateCorrectedFvPatchField<Type>::mappedFlowRateCorrectedFvPatchField
(
    const mappedFlowRateCorrectedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFlowRateCorrectedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFlowRateCorrectedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFlowRateCorrectedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // First get the mapped donor field from the parent machinery
    this->operator==(this->mappedField());

    // This BC is compiled only for vector, so it is safe to treat as vectorField
    vectorField& U = static_cast<vectorField&>(*this);

    const fvPatch& p = this->patch();
    const vectorField nHat(p.nf());
    const scalarField magSf(p.magSf());

    const volScalarField& alpha =
        this->db().template lookupObject<volScalarField>("alpha.water");

    const labelUList& faceCells = p.faceCells();

    scalar QwRec = 0.0;
    scalar Awet  = 0.0;

    forAll(U, facei)
    {
        scalar a = alpha[faceCells[facei]];
        a = min(max(a, scalar(0)), scalar(1));

        const scalar Un = (U[facei] & nHat[facei]);

        QwRec += a * Un * magSf[facei];
        Awet  += a * magSf[facei];
    }

    QwRec = returnReduce(QwRec, sumOp<scalar>());
    Awet  = returnReduce(Awet, sumOp<scalar>());

    scalar dUn = 0.0;

    if (Awet > SMALL)
    {
        // Patch normals point outward, so inflow is negative flux
        dUn = (-Qtarget_ - QwRec)/Awet;
    }

    forAll(U, facei)
    {
        scalar a = alpha[faceCells[facei]];
        a = min(max(a, scalar(0)), scalar(1));

        const vector n = nHat[facei];
        const scalar Un = (U[facei] & n);
        const vector Ut = U[facei] - Un*n;

        // Preserve tangential donor structure, correct only normal component
        U[facei] = Ut + (Un + a*dUn)*n;
    }

    if (debug)
    {
        auto limits = gMinMax(U);
        auto avg = gAverage(U);

        Info<< "mappedFlowRateCorrected on field:" << this->internalField().name()
            << " patch:" << this->patch().name()
            << " Qtarget:" << Qtarget_
            << " QwRec:" << QwRec
            << " Awet:" << Awet
            << " dUn:" << dUn
            << " avg:" << avg
            << " min:" << limits.min()
            << " max:" << limits.max()
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void Foam::mappedFlowRateCorrectedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    mappedFlowRateCorrectedPatchFieldBase<Type>::write(os);
    os.writeKeyword("Qtarget") << Qtarget_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaMin") << alphaMin_ << token::END_STATEMENT << nl;
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //