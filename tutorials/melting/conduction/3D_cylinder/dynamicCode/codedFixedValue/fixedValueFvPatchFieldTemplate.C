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
#line 120437 "/home/davide/OpenFOAM/davide-6/OpenFOAM-nrg/tutorials/melting/conduction/3D_cylinder/0/T.boundaryField.top"
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
    // SHA1 = f1e9ce8774add6c652c3e9315f24083f60df52b2
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void codedFixedValue_f1e9ce8774add6c652c3e9315f24083f60df52b2(bool load)
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
    "f1e9ce8774add6c652c3e9315f24083f60df52b2";


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
        Info<<"construct codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2"
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
        Info<<"construct codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2"
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
        Info<<"construct codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2"
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
        Info<<"construct codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2"
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
        Info<<"construct codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

codedFixedValueFixedValueFvPatchScalarField::
~codedFixedValueFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2\n";
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
        Info<<"updateCoeffs codedFixedValue sha1: f1e9ce8774add6c652c3e9315f24083f60df52b2\n";
    }

//{{{ begin code
    #line 120450 "/home/davide/OpenFOAM/davide-6/OpenFOAM-nrg/tutorials/melting/conduction/3D_cylinder/0/T.boundaryField.top"
const fvPatch& boundaryPatch = patch();
        	const vectorField& Cf = boundaryPatch.Cf();
        	scalarField& field = *this;
        	
        	// thermophysical parameters of the PCM (water)
        	const double Cp	= 578;
        	const double rho = 8860;
        	const double k = 23.6;

        	const double L 	= 0.2;

        	const double T_H = 1273;
        	const double T_0 = 836;

        	const double alfa = k/(rho*Cp);
        	
        	const scalar t = this->db().time().value(); //time
        	
        	forAll(Cf, faceI)
        	{
        		const scalar z = Cf[faceI].z(); // z coordinate of the face i
        		const scalar r = L - z;

        		const scalar denominator = 2*std::sqrt(alfa*t);  // liquid front position
        		
        		field[faceI] = T_H - (T_H-T_0)*std::erf(r/2/std::sqrt(alfa*t));
        	}
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

