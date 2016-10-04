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

//Includes routines that carry out the main-grid updates.

#include "update_no_dispersion.h"
#include "update_full.h"

extern bool dispersion_exists;


//**********************************
//	UPDATE THE ELECTRIC FIELD
//**********************************
void updateE(const int& n)
{
	if (dispersion_exists)
	{
		updateE_full(n);
	}
	else
	{
		updateE_no_dispersion(n);
	}
}

//**********************************
//	UPDATE THE MAGNETIC FIELD
//**********************************
void updateH(const int& n)
{
	if (dispersion_exists)
	{
		updateH_full(n);
	}
	else
	{
		updateH_no_dispersion(n);
	}
}
