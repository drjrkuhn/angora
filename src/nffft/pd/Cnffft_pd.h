/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CNFFFT_PD_H
#define CNFFFT_PD_H

//Declaration of the PHASOR_DOMAIN near-field-to-far-field transformer (NFFFT) object "Cnffft_pd", which represents a collection of transformer objects "Ctr_pd".

//only the declaration of Ctr_pd needed: use forward declaration
class Ctr_pd;
//only the declaration of TrDataType_pd needed: use forward declaration
class TrDataType_pd;

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>


class Cnffft_pd
{
 public:
	 int AddTransformer(const TrDataType_pd& MyData,
		 const string& InputFarFieldFileName = "");	//adds a transformer

	 void UpdateFarField(const int& n);		//Update far-field arrays using the current E,H on the virtual surface
	 void ConstructFarField();	//Construct far-field from the far-field arrays

	 int NumberOfTransformers()
	 {
		 return Transformers.size();
	 }

 private:
	 vector<boost::shared_ptr<Ctr_pd> > Transformers;	//array of pointers to transformer objects
};

#endif
