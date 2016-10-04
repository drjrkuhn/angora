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

#ifndef CMAT_EXCP_H
#define CMAT_EXCP_H

//the exceptions thrown by Cmat

//base Angora exception class
#include "angora_excp.h"


class MaterialPropertyDoesNotExist: public AngoraException
{// exception raised when the material does not have the specified constitutive parameter
public:
  MaterialPropertyDoesNotExist(const string& matprop) :_matprop(matprop){};
  virtual ~MaterialPropertyDoesNotExist() throw() {};

  virtual const string getError() const
  {//error message
  	ostringstream _msgstr;
	_msgstr << "Error in Cmat (material class): material property \"" << _matprop << "\" does not exist in material";
  	return _msgstr.str();
  }
protected:
 const string _matprop;
};

#endif // CMAT_EXCP_H
