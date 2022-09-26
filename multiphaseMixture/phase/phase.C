/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "phase.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phase::phase
(
    const word& phaseName,
    const dictionary& phaseDict,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    name_(phaseName),
    phaseDict_(phaseDict),
    nuModel_
    (
        viscosityModel::New
        (
            IOobject::groupName("nu", phaseName),
            phaseDict_,
            U,
            phi
        )
    ),
    rho_("rho", dimDensity, phaseDict_),
    gFlag_("gFlag", dimless, phaseDict_),
	U_(U),
	magGradAlpha_
    (
        IOobject
        (
            IOobject::groupName("magGradAlpha", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		U.mesh(),
		dimless
    ),
	specificStrainRate_
    (
        IOobject
        (
            IOobject::groupName("specificStrainRate", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		U.mesh(),
		dimless/dimTime
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    NotImplemented;
    return nullptr;
}

void Foam::phase::calcMagGradAlpha()
{
	volScalarField& alpha = *this;
	magGradAlpha_ = mag(fvc::grad(alpha))*
		dimensionedScalar(dimensionSet(0,1,0,0,0),1);
	//magGradAlpha_ /= magGradAlpha_.weightedAverage(U_.mesh().V());//average();
	//magGradAlpha_ /= magGradAlpha_.average();
	//magGradAlpha_ *= (scalar(1) - alpha);
	//magGradAlpha_.clip(0, 1);
	//magGradAlpha_ -= scalar(min(magGradAlpha_));
	//magGradAlpha_ /= scalar(max(magGradAlpha_));
}

void Foam::phase::calcSpecificStrainRate()
{
	calcMagGradAlpha();
	volScalarField& alpha = *this;
	//invariantII(strainRateTensor2Inv_, symm(fvc::grad(U_))*dimensionedScalar(dimensionSet(0,-1,0,0,0),1));
	//invariantII(strainRateTensor2Inv_, symm(fvc::grad(U_)));
	specificStrainRate_ = nuModel_->strainRate();//*dimensionedScalar(dimless*dimTime, 1);//strainRate();
	//specificStrainRate_ -= Foam::min(strainRateTensor2Inv_);
	//specificStrainRate_ -= scalar(1);
	//specificStrainRate_.clip(0, 1);
	specificStrainRate_ *= magGradAlpha_ * alpha;
	//specificStrainRate_ /= dimensionedScalar(dimless, 200);
	//specificStrainRate_ -= dimensionedScalar(dimless/dimTime, 1);
    //specificStrainRate_.clip(0, 1);
	//specificStrainRate_ /= Foam::max(strainRateTensor2Inv_);
}

void Foam::phase::correct()
{
    nuModel_->correct();
}


bool Foam::phase::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    if (nuModel_->read(phaseDict_))
    {
        phaseDict_.readEntry("rho", rho_);
        phaseDict_.readEntry("gFlag", gFlag_);

        return true;
    }

    return false;
}


// ************************************************************************* //
