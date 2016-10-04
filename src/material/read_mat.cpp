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

// reading the config options that specify the material definitions

#include "headers.h"

#include "read_mat.h"

#include "Cmats.h"


void read_mat(Cmats &Materials, const Config& fdtdconfig, const Config& validsettings)
{
	string material_setting_path = "Materials";
	if (fdtdconfig.exists(material_setting_path))
	{
		Setting& Materiallistsettings = read_list_from_group(fdtdconfig.getRoot(),material_setting_path);

		int num_of_newmaterials = Materiallistsettings.getLength();
		for (int newmaterialindex=0; newmaterialindex<num_of_newmaterials; newmaterialindex++)
		{
			Setting& Newmaterialsettings = Materiallistsettings[newmaterialindex];	//go to the newmaterialindex'th new material setting
			//check group for invalid settings
			CheckAngoraGroupSetting(Newmaterialsettings,validsettings);

			if (SettingEnabledForGrid(Newmaterialsettings))		//apply only if enabled for this grid
			{
				//create the material object
				Cmat NewMaterial;

				//is the material transparent? (are the unspecified parameters assigned default values, or left untouched)
				bool transparent;
				if (!read_optional_value_from_group<bool>(Newmaterialsettings,"transparent",transparent))
				{
					transparent = false; //by default, unspecified parameters are assigned default values
				}

				float rel_permittivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"rel_permittivity",rel_permittivity))
				{
					if (!transparent)
					{
						NewMaterial.set_eps(1.0); //default: eps_r = 1
					}// if transparent, don't set the permittivity
				}
				else
				{
					NewMaterial.set_eps(rel_permittivity);
				}

				float rel_permeability;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"rel_permeability",rel_permeability))
				{
					if (!transparent)
					{
						NewMaterial.set_mu(1.0); //default: mu_r = 1
					}
				}
				else
				{
					NewMaterial.set_mu(rel_permeability);
				}

				float electric_conductivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"electric_conductivity",electric_conductivity))
				{
					if (!transparent)
					{
						NewMaterial.set_cond_e(0); //default: sigma = 0
					}
				}
				else
				{
					NewMaterial.set_cond_e(electric_conductivity);
				}

				float magnetic_conductivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"magnetic_conductivity",magnetic_conductivity))
				{
					if (!transparent)
					{
						NewMaterial.set_cond_h(0); //default: sigma = 0
					}
				}
				else
				{
					NewMaterial.set_cond_h(magnetic_conductivity);
				}

				float drude_pole_frequency;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"drude_pole_frequency",drude_pole_frequency))
				{
					if (!transparent)
					{
						NewMaterial.set_omega_p(0); //default: omega_p = 0
					}
				}
				else
				{
					NewMaterial.set_omega_p(drude_pole_frequency);
				}

				float drude_pole_relaxation_time;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"drude_pole_relaxation_time",drude_pole_relaxation_time))
				{
					if (!transparent)
					{
						NewMaterial.set_tau_r(0); //default: tau_r = 0
					}
				}
				else
				{
					NewMaterial.set_tau_r(drude_pole_relaxation_time);
				}

				string material_tag;
				read_value_from_group<string>(Newmaterialsettings,"material_tag",material_tag);

//				AddIsotropicMaterial(NewMaterialId,rel_permittivity,rel_permeability,electric_conductivity,magnetic_conductivity);
				//Add the new material to the material-collector object (even in check mode, since other objects refer to these materials)
				Materials.CreateMaterial(NewMaterial,material_tag);
			}
		}
	}//if exists
}
