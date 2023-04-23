/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 33 "/home/davide/OpenFOAM/davide-6/OpenFOAM-nrg/tutorials/melting/conduction/2D_melting/0/T.boundaryField.north"
#include "fvCFD.H"
        	#include <cmath>
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void codedFixedValue_18297ff26d8fbb5f6cf3736e40246ccc7ac84fce(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    codedFixedValueFixedValueFvPatchScalarField
);


const char* const codedFixedValueFixedValueFvPatchScalarField::SHA1sum =
    "18297ff26d8fbb5f6cf3736e40246ccc7ac84fce";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

codedFixedValueFixedValueFvPatchScalarField::
codedFixedValueFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce"
            " from patch/DimensionedField\n";
    }
}


codedFixedValueFixedValueFvPatchScalarField::
codedFixedValueFixedValueFvPatchScalarField
(
    const codedFixedValueFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce"
            " from patch/DimensionedField/mapper\n";
    }
}


codedFixedValueFixedValueFvPatchScalarField::
codedFixedValueFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce"
            " from patch/dictionary\n";
    }
}


codedFixedValueFixedValueFvPatchScalarField::
codedFixedValueFixedValueFvPatchScalarField
(
    const codedFixedValueFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce"
            " as copy\n";
    }
}


codedFixedValueFixedValueFvPatchScalarField::
codedFixedValueFixedValueFvPatchScalarField
(
    const codedFixedValueFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

codedFixedValueFixedValueFvPatchScalarField::
~codedFixedValueFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void codedFixedValueFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs codedFixedValue sha1: 18297ff26d8fbb5f6cf3736e40246ccc7ac84fce\n";
    }

//{{{ begin code
    #line 46 "/home/davide/OpenFOAM/davide-6/OpenFOAM-nrg/tutorials/melting/conduction/2D_melting/0/T.boundaryField.north"
const fvPatch& boundaryPatch = patch();
        	const vectorField& Cf = boundaryPatch.Cf();
        	scalarField& field = *this;
        	
        	
        	// geometrical parameters of the plane of melting
        	const double x0		= 0.10;
        	const double y0		= 0.10;
        	const double z0		= 0.;
        	const double Theta	= 45*3.14159/180;;
        	const double Phi        = 90*3.14159/180;;
        	
        	// thermophysical parameters of the PCM (water)
        	const double Cp_l	= 4.186e3;
        	const double Cp_s 	= 2.04e3;
        	const double rho 	= 1000;
        	const double k 		= 0.6;
        	const double Tm 	= 273.15;

        	const double T_l 	= 293.15;
        	const double T_s 	= 260;

        	const double alfa_l 	= k/(rho*Cp_l);
        	const double alfa_s 	= k/(rho*Cp_s);
        	const double nu		= std::sqrt(alfa_l/alfa_s);
        	
        	const scalar lambda	= 0.30492;
        	        	
        	
        	const scalar t = this->db().time().value(); //time
        	
        	forAll(Cf, faceI)
        	{
        		const scalar x = Cf[faceI].x(); // x coordinate of the face i
        		const scalar y = Cf[faceI].y(); // y coordinate of the face i
        		const scalar z = Cf[faceI].z(); // y coordinate of the face i
        		const scalar r = std::abs( (x-x0)*std::cos(Theta)*std::sin(Phi) + (y-y0)*std::sin(Theta)*std::sin(Phi) + (z-z0)*std::cos(Phi) ); // distance from the plane
        					

        		const scalar front = 2*lambda*std::sqrt(alfa_l*t);  // liquid front position
        		
        		if (r<front)
        		{
        		  field[faceI] = T_l - (T_l-Tm)*std::erf(r/2/std::sqrt(alfa_l*t))/std::erf(lambda);
        		}
        		
        		else if (r>=front)
        		{
        		  field[faceI] = T_s + (Tm-T_s)*std::erfc(r/2/std::sqrt(alfa_s*t))/std::erfc(lambda*nu);
        		}        		  
        		
        	}
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

