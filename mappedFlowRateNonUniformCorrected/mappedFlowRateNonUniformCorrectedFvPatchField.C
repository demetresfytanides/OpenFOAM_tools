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

\*---------------------------------------------------------------------------*/

#include "mappedFlowRateNonUniformCorrectedFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this),
    Qtarget_(0.0),
    alphaMin_(1e-8),
    gammaMin_(0.5),
    gammaMax_(2.0)
{}


template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, dict, *this),
    Qtarget_(dict.lookupOrDefault<scalar>("Qtarget",  0.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 1e-8)),
    gammaMin_(dict.lookupOrDefault<scalar>("gammaMin", 0.5)),
    gammaMax_(dict.lookupOrDefault<scalar>("gammaMax", 2.0))
{
    if (!dict.found("value"))
    {
        this->evaluate();
    }
}


template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
(
    const mappedFlowRateNonUniformCorrectedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_),
    gammaMin_(ptf.gammaMin_),
    gammaMax_(ptf.gammaMax_)
{}


template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
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
    alphaMin_(1e-8),
    gammaMin_(0.5),
    gammaMax_(2.0)
{}


template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
(
    const mappedFlowRateNonUniformCorrectedFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_),
    gammaMin_(ptf.gammaMin_),
    gammaMax_(ptf.gammaMax_)
{}


template<class Type>
Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::
mappedFlowRateNonUniformCorrectedFvPatchField
(
    const mappedFlowRateNonUniformCorrectedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedFlowRateCorrectedPatchFieldBase<Type>(*this, *this, ptf),
    Qtarget_(ptf.Qtarget_),
    alphaMin_(ptf.alphaMin_),
    gammaMin_(ptf.gammaMin_),
    gammaMax_(ptf.gammaMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::updateCoeffs()
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

    // First pass: compute current water flux through the inlet using the
    // donor velocity that was just mapped in.
    scalar QwRec = 0.0;
    scalar Awet  = 0.0;

    forAll(U, facei)
    {
        scalar a = alpha[faceCells[facei]];
        a = min(max(a, scalar(0)), scalar(1));

        // Faces below the wet threshold do not contribute to Qw or Awet
        const scalar aw = (a >= alphaMin_) ? a : scalar(0);

        const scalar Un = (U[facei] & nHat[facei]);

        QwRec += aw * Un * magSf[facei];
        Awet  += aw * magSf[facei];
    }

    QwRec = returnReduce(QwRec, sumOp<scalar>());
    Awet  = returnReduce(Awet,  sumOp<scalar>());

    // Multiplicative scale factor on the normal component.
    // Patch normals point outward, so inflow corresponds to QwRec < 0,
    // and the target is gamma*QwRec = -Qtarget.
    scalar gamma = 1.0;

    if (mag(QwRec) > SMALL)
    {
        gamma = -Qtarget_/QwRec;
    }

    // Clip to keep the rescale physical even when QwRec is small or has
    // the wrong sign on a transient.
    gamma = max(min(gamma, gammaMax_), gammaMin_);

    // Second pass: rebuild U as
    //   U = Ut + gammaEff * Un * n,
    // with gammaEff = (1 - a) + a*gamma so that fully-dry faces (a=0) are
    // left at the donor field and fully-wet faces (a=1) get the full
    // multiplicative correction.
    forAll(U, facei)
    {
        scalar a = alpha[faceCells[facei]];
        a = min(max(a, scalar(0)), scalar(1));

        const vector n  = nHat[facei];
        const scalar Un = (U[facei] & n);
        const vector Ut = U[facei] - Un*n;

        const scalar gammaEff = (scalar(1) - a) + a*gamma;

        // Preserve tangential donor structure, rescale only normal component
        U[facei] = Ut + gammaEff*Un*n;
    }

    if (debug)
    {
        // Optional: report the post-correction water flux
        scalar QwAfter = 0.0;
        forAll(U, facei)
        {
            scalar a = alpha[faceCells[facei]];
            a = min(max(a, scalar(0)), scalar(1));
            const scalar aw = (a >= alphaMin_) ? a : scalar(0);
            QwAfter += aw * (U[facei] & nHat[facei]) * magSf[facei];
        }
        QwAfter = returnReduce(QwAfter, sumOp<scalar>());

        auto limits = gMinMax(U);
        auto avg = gAverage(U);

        Info<< "mappedFlowRateNonUniformCorrected on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << " Qtarget:" << Qtarget_
            << " QwRec:"   << QwRec
            << " Awet:"    << Awet
            << " gamma:"   << gamma
            << " QwAfter:" << QwAfter
            << " avg:"     << avg
            << " min:"     << limits.min()
            << " max:"     << limits.max()
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedFlowRateNonUniformCorrectedFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    mappedPatchBase::write(os);
    mappedFlowRateCorrectedPatchFieldBase<Type>::write(os);
    os.writeKeyword("Qtarget")  << Qtarget_  << token::END_STATEMENT << nl;
    os.writeKeyword("alphaMin") << alphaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("gammaMin") << gammaMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("gammaMax") << gammaMax_ << token::END_STATEMENT << nl;
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
