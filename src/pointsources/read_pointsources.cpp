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

//Includes the routine that reads the point-source definitions

#include "headers.h"

#include "read_pointsources.h"

//definition of Cpointsources needed
#include "Cpointsources.h"

//definition of Cwfs needed
#include "waveforms/Cwfs.h"

extern bool check_mode;

extern double dx;

extern string config_filename;

extern int OriginX,OriginY,OriginZ;


void read_pointsources(Cpointsources &PointSources, const Config& fdtdconfig, const Config& validsettings, const Cwfs& Waveforms)
{
	//Read point-source settings
	string point_sources_setting_path = "PointSources";
	if (fdtdconfig.exists(point_sources_setting_path))
	{
		Setting& PointSourceCollec = read_list_from_group(fdtdconfig.getRoot(),point_sources_setting_path);

		int num_of_pointsources = PointSourceCollec.getLength();
		for (int pointsourceindex=0; pointsourceindex<num_of_pointsources; pointsourceindex++)
		{
			Setting& Pointsourcesettings = PointSourceCollec[pointsourceindex];		//go to the pointsourceindex'th point-source setting
			//check group for invalid settings
			CheckAngoraGroupSetting(Pointsourcesettings,validsettings);

			if (SettingEnabledForGrid(Pointsourcesettings))		//apply only if enabled for this grid
			{
				string source_orientation;
				read_value_from_group<string>(Pointsourcesettings,"source_orientation",source_orientation);

				int cell_index_x,cell_index_y,cell_index_z;
				double coord_x,coord_y,coord_z;
				if (read_length_from_group<double>(Pointsourcesettings,"coord_x",coord_x))
				{//returned true: given in meters
					cell_index_x = (int)round(coord_x/dx)+OriginX; //relative to grid origin
				}
				else
				{//returned false: given in grid cells
					cell_index_x = (int)round(coord_x)+OriginX; //relative to grid origin
				}
				if (read_length_from_group<double>(Pointsourcesettings,"coord_y",coord_y))
				{//returned true: given in meters
					cell_index_y = (int)round(coord_y/dx)+OriginY; //relative to grid origin
				}
				else
				{//returned false: given in grid cells
					cell_index_y = (int)round(coord_y)+OriginY; //relative to grid origin
				}
				if (read_length_from_group<double>(Pointsourcesettings,"coord_z",coord_z))
				{//returned true: given in meters
					cell_index_z = (int)round(coord_z/dx)+OriginZ; //relative to grid origin
				}
				else
				{//returned false: given in grid cells
					cell_index_z = (int)round(coord_z)+OriginZ; //relative to grid origin
				}

				//read the extra current moment
				double j0;
				if (!read_optional_value_from_group<double>(Pointsourcesettings,"j_0",j0))
				{
					j0 = 1.0;
				}

				//read the index of the waveform in the Waveforms object that will represent the current moment (dimensions:current)
				string waveform_tag;
				read_value_from_group<string>(Pointsourcesettings,"waveform_tag",waveform_tag);
				const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
				//If not in check mode, add the point source to the PointSources object
				if (!check_mode)
				{
					PointSources.AddElectricDipole(cell_index_x,cell_index_y,cell_index_z,source_orientation,j0,waveformptr);
				}
			}
		}
	}
}
