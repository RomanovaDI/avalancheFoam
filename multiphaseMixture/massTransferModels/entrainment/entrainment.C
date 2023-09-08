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
    const phase& phase1,//soil
    const phase& phase2 //entrained
)
:
    massTransferModel(interfaceDict, phase1, phase2),
	dict_(interfaceDict.optionalSubDict("entrainmentCoeffs")),
	entrCoeff_("entrCoeff", dimless, dict_),
	breakingPoint_("breakingPoint", dimless, dict_)
{
	Info << "entrCoeff = " << entrCoeff_ << endl; 
	Info << "breakingPoint = " << breakingPoint_ << endl; 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//phase1 - from, phase2 - to
Foam::tmp<Foam::volScalarField> Foam::massTransferModels::entrainment::K() const
{
	/*tmp<volScalarField> limitedSpecificStrainRate = phase1_.specificStrainRate();
	//limitedSpecificStrainRate = limitedSpecificStrainRate.ref().clip(0,1);
	limitedSpecificStrainRate = Foam::min(limitedSpecificStrainRate / breakingPoint_, dimensionedScalar(dimless/dimTime, 1.0));
	forAll(limitedSpecificStrainRate.ref(), i)
	{
		limitedSpecificStrainRate.ref()[i] = floor(limitedSpecificStrainRate.ref()[i]+1e-04);
	}
	//const volScalarField& alpha1 = phase1_();
	tmp<volScalarField> alpha1 = phase1_;
	limitedSpecificStrainRate = limitedSpecificStrainRate * residualPhaseFraction_ / (Foam::min(alpha1, residualPhaseFraction_) + SMALL);
	//forAll (limitedSpecificStrainRate, i)
	//{
	//	if (alpha1[i] < residualPhaseFraction())
	//		limitedSpecificStrainRate[i] *= VGREAT;
	//}
	//limitedSpecificStrainRate -= dimencionedScalar(dimless/dimTime, 1);
	//limitedSpecificStrainRate.clip(0,1);*/
	tmp<volScalarField> limitedSpecificStrainRate = phase1_.specificStrainRate();
	tmp<volScalarField> alpha1 = phase1_;
	tmp<volScalarField> erDepth = phase1_.magU() / (dimensionedScalar(dimless/dimTime, SMALL) + limitedSpecificStrainRate);
	limitedSpecificStrainRate = \
		Foam::min(\
			Foam::max(\
				(limitedSpecificStrainRate - breakingPoint_*dimensionedScalar(dimless/dimTime, 1.0)) * GREAT,\
				dimensionedScalar(dimless/dimTime, 0.0)\
			),\
			dimensionedScalar(dimless/dimTime, 1.0)\
		);
	//		/ (phase1_.magGradStrainRate()+SMALL*dimensionedScalar(dimless/(dimLength*dimTime), 1.0));
	/////////////////////////////////////////
	tmp<volScalarField> Kvalue(
			//Foam::min(
				erDepth * limitedSpecificStrainRate * entrCoeff_ * dimensionedScalar(dimless/dimLength, 1.0)//,
			//	alpha1
			//)
			//phase1_.specificStrainRate()/
			//breakingPoint_*
			//dimensionedScalar(dimensionSet(0,0,0,0,0),-1)*
		);
	return Kvalue;
}


// ************************************************************************* /
