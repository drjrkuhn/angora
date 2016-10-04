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

//function that checks for dispersion in the grid

#include "headers.h"

extern bool dispersion_exists;
extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;

extern int Ex_min_index_in_x,Ex_max_index_in_x,Ex_min_index_in_y,Ex_max_index_in_y,Ex_min_index_in_z,Ex_max_index_in_z;
extern int Ey_min_index_in_x,Ey_max_index_in_x,Ey_min_index_in_y,Ey_max_index_in_y,Ey_min_index_in_z,Ey_max_index_in_z;
extern int Ez_min_index_in_x,Ez_max_index_in_x,Ez_min_index_in_y,Ez_max_index_in_y,Ez_min_index_in_z,Ez_max_index_in_z;
extern int Hx_min_index_in_x,Hx_max_index_in_x,Hx_min_index_in_y,Hx_max_index_in_y,Hx_min_index_in_z,Hx_max_index_in_z;
extern int Hy_min_index_in_x,Hy_max_index_in_x,Hy_min_index_in_y,Hy_max_index_in_y,Hy_min_index_in_z,Hy_max_index_in_z;
extern int Hz_min_index_in_x,Hz_max_index_in_x,Hz_min_index_in_y,Hz_max_index_in_y,Hz_min_index_in_z,Hz_max_index_in_z;

namespace{
	int i,j,k;
};


void check_dispersion()
{
	dispersion_exists = false;

	for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
		for (j=Ex_min_index_in_y; j<=Ex_max_index_in_y; j++){
			for (k=Ex_min_index_in_z; k<=Ex_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ex_position(i,j,k))
				{
					dispersion_exists = true;
					return;
				}
			}
		}
	}

	for (i=Ey_min_index_in_x; i<=Ey_max_index_in_x; i++){
		for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
			for (k=Ey_min_index_in_z; k<=Ey_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ey_position(i,j,k))
				{
					dispersion_exists = true;
					return;
				}
			}
		}
	}

	for (i=Ez_min_index_in_x; i<=Ez_max_index_in_x; i++){
		for (j=Ez_min_index_in_y; j<=Ez_max_index_in_y; j++){
			for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ez_position(i,j,k))
				{
					dispersion_exists = true;
					return;
				}
			}
		}
	}
}
