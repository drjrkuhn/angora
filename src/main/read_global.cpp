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

//Function that reads some global variables from the config file

#include "headers.h"

#include "read_global.h"

extern string config_filename;

extern double courant,dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML,NSTEPS;
extern int OriginX,OriginY,OriginZ;

extern double max_field_value;	//maximum electric field value encountered in grid
extern bool max_field_value_set_in_configfile;		//is the max. field value fixed in the config file?
extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;


void read_global(const Config& fdtdconfig, const Config& validsettings)
//Reads some global variables
{
	//Read and override maximum field value, if specified in the config file
	if (read_optional_value_from_group<double>(fdtdconfig.getRoot(),"max_field_value",max_field_value))
	{
		max_field_value_set_in_configfile = true;	//max. field value set in config file
	}
	//Read and override dB accuracy, if specified in the config file
	if (read_optional_value_from_group<double>(fdtdconfig.getRoot(),"dB_accuracy",dB_accuracy))
	{
		dB_accuracy_set_in_configfile = true;	//dB accuracy set in config file
	}

	//Read courant stability factor
	read_value_from_group<double>(fdtdconfig.getRoot(),"courant",courant);

	//Read grid spacing
	read_value_from_group<double>(fdtdconfig.getRoot(),"dx",dx);

	dt = courant*dx/c/1.73205;		//Time step

	//Read number of grid cells in the x,y,z directions
	double grid_dimension_x,grid_dimension_y,grid_dimension_z;
	if (read_length_from_group<double>(fdtdconfig.getRoot(),"grid_dimension_x",grid_dimension_x))
	{
		NCELLS_X = (int)round(grid_dimension_x/dx);
	}
	else
	{
		NCELLS_X = (int)round(grid_dimension_x);
	}
	if (read_length_from_group<double>(fdtdconfig.getRoot(),"grid_dimension_y",grid_dimension_y))
	{
		NCELLS_Y = (int)round(grid_dimension_y/dx);
	}
	else
	{
		NCELLS_Y = (int)round(grid_dimension_y);
	}
	if (read_length_from_group<double>(fdtdconfig.getRoot(),"grid_dimension_z",grid_dimension_z))
	{
		NCELLS_Z = (int)round(grid_dimension_z/dx);
	}
	else
	{
		NCELLS_Z = (int)round(grid_dimension_z);
	}

	//Read PML thickness
	double pml_thickness;
	if (read_length_from_group<double>(fdtdconfig.getRoot(),"pml_thickness",pml_thickness))
	{
		NPML = (int)round(pml_thickness/dx);
	}
	else
	{
		NPML = (int)round(pml_thickness);
	}


	//Read number of time steps
	read_value_from_group<int>(fdtdconfig.getRoot(),"num_of_time_steps",NSTEPS);

	//Read x index of the origin cell
	double origin_x,origin_y,origin_z;
	try{
		if (read_length_from_group<double>(fdtdconfig.getRoot(),"origin_x",origin_x))
		{
			OriginX += (int)(origin_x/dx)+1; //actual grid cell index is (coordinate+1)
		}
		else
		{
			OriginX += (int)origin_x+1; //actual grid cell index is (coordinate+1)
		}
	}
	catch (AngoraSettingNotFoundException&)
	{
		OriginX = (NCELLS_X+2*NPML)/2+1;
	}

	try{
		if (read_length_from_group<double>(fdtdconfig.getRoot(),"origin_y",origin_y))
		{
			OriginY += (int)(origin_y/dx)+1; //actual grid cell index is (coordinate+1)
		}
		else
		{
			OriginY += (int)origin_y+1; //actual grid cell index is (coordinate+1)
		}
	}
	catch (AngoraSettingNotFoundException&)
	{
		OriginY = (NCELLS_Y+2*NPML)/2+1;
	}

	try{
		if (read_length_from_group<double>(fdtdconfig.getRoot(),"origin_z",origin_z))
		{
			OriginZ += (int)(origin_z/dx)+1; //actual grid cell index is (coordinate+1)
		}
		else
		{
			OriginZ += (int)origin_z+1; //actual grid cell index is (coordinate+1)
		}
	}
	catch (AngoraSettingNotFoundException&)
	{
		OriginZ = (NCELLS_Z+2*NPML)/2+1;
	}
}
