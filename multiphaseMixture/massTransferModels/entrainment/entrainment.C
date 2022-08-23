/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "entrainment.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferModels
{
    defineTypeNameAndDebug(entrainment, 0);

    addToRunTimeSelectionTable
    (
        massTransferModel,
        entrainment,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferModels::entrainment::entrainment
(
    const dictionary& interfaceDict,
    const phase& phase1,
    const phase& phase2
)
:
    massTransferModel(interfaceDict, phase1, phase2),
	dict_(interfaceDict.optionalSubDict("entrainmentCoeffs")),
	entrCoeff_("entrCoeff", dimless, dict_)
{
	Info << "entrCoeff = " << entrCoeff_ << endl; 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//phase1 - from, phase2 - to
Foam::tmp<Foam::volScalarField> Foam::massTransferModels::entrainment::K() const
{
	//volScalarField& alpha2 = phase2_;
	tmp<volScalarField> Kvalue(
		//min(
		//	dimensionedScalar(dimensionSet(0,0,-2,0,0), 1),
		//	phase1_.strainRateTensor2Inv()
		//	)*
		phase1_.strainRateTensor2Inv()*
		//phase1_.magGradAlpha()*
		dimensionedScalar(dimensionSet(0,0,-1,0,0),-1)*
		entrCoeff_
		);
    
	return Kvalue;
}


// ************************************************************************* /
