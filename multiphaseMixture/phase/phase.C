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
	U_(U),
	strainRateTensor2Inv_
    (
        IOobject
        (
            IOobject::groupName("strainRateTensor2Inv", phaseName),
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		U.mesh(),
		dimensionSet(0,0,-2,0,0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phase> Foam::phase::clone() const
{
    NotImplemented;
    return nullptr;
}

void Foam::phase::calcStrainRateTensor2Inv()
{
    Info<< "phaseModel.calcStrainRateTensor2Inv()" << nl << endl;
	volScalarField& alpha = *this;
	invariantII(strainRateTensor2Inv_, symm(fvc::grad(U_)));
	strainRateTensor2Inv_ *= alpha;
    strainRateTensor2Inv_.clip(0, 1);
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

        return true;
    }

    return false;
}


// ************************************************************************* //
