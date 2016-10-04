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

//Definition of the class "Cmats" that holds the named-material information

#include "Cmats.h"


bool Cmats::MaterialTagExists(const string& material_tag)
{//returns true if a material with the given tag exists.
	return (NamedMaterials.find(material_tag)!=NamedMaterials.end());
}

bool Cmats::lookupMaterialWithTag(const string& material_tag, Cmat& material) const
{//copies the material with tag into "material" if it exists.
//returns false if the string tag is not found
	map<string,Cmat>::const_iterator map_it = NamedMaterials.find(material_tag);
	if (map_it==NamedMaterials.end())
	{//string tag does not correspond to any material, return false
		return false;
	}
	else
	{
		material = map_it->second;  //second element of the pair<string,MaterialId> object is the material identifier
		return true;
	}
}

const Cmat  Cmats::operator[] (const string& material_tag) const
{
	map<string,Cmat>::const_iterator map_it = NamedMaterials.find(material_tag);
	if (map_it==NamedMaterials.end())
	{//string tag does not correspond to any material, throw exception
		throw NamedMaterialNotFoundException(material_tag);
	}
	else
	{
		return map_it->second;
	}
}

void Cmats::CreateMaterial(const Cmat& material, const string& material_tag)
{//creates a planar sheet with the given string tag
	if (!MaterialTagExists(material_tag))
	{
		//add the material with the given tag
		NamedMaterials.insert(pair<string,Cmat>(material_tag,material));
	}
	else
	{//string tag already exists, throw exception
		throw NamedMaterialExistsException(material_tag);
	}
}
