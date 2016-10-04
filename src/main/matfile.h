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

#ifndef MATFILE_H
#define MATFILE_H

//template declaration for the function that reads a material region from a file

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

#include "material/Cmat_types.h"

#include "shape/Cshape.h"

#include <fstream>

#ifndef ANGORA_MAX_NEWMAT
#define ANGORA_MAX_NEWMAT 1000
#endif

extern double dx,dt;

extern int rank;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

//extern int num_of_distinct_eps_x,num_of_distinct_eps_y,num_of_distinct_eps_z;
//extern int num_of_distinct_mu_x,num_of_distinct_mu_y,num_of_distinct_mu_z;
//extern int num_of_distinct_cond_e_x,num_of_distinct_cond_e_y,num_of_distinct_cond_e_z;
//extern int num_of_distinct_cond_h_x,num_of_distinct_cond_h_y,num_of_distinct_cond_h_z;

extern Array<update_coeff_type,3> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<eps_x_type,3> eps_x_indices;
extern Array<eps_y_type,3> eps_y_indices;
extern Array<eps_z_type,3> eps_z_indices;
extern Array<mu_x_type,3> mu_x_indices;
extern Array<mu_y_type,3> mu_y_indices;
extern Array<mu_z_type,3> mu_z_indices;
extern Array<cond_e_x_type,3> cond_e_x_indices;
extern Array<cond_e_y_type,3> cond_e_y_indices;
extern Array<cond_e_z_type,3> cond_e_z_indices;
extern Array<cond_h_x_type,3> cond_h_x_indices;
extern Array<cond_h_y_type,3> cond_h_y_indices;
extern Array<cond_h_z_type,3> cond_h_z_indices;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;


template<typename MatType> //data type for the material property
void PlaceMaterialRegionFromFile(const string& MaterialFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor, const string& constitutive_param_type, const int& max_number_of_new_materials = ANGORA_MAX_NEWMAT, const_Cshape_shared_ptr shape_mask = const_Cshape_shared_ptr(new Cuniverse()))
{//Reads rectangular-prism-shaped dielectric region from file and places into grid
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the region (measured in cells from the back-left-lower corner of the grid)

	if ((constitutive_param_type!="rel_permittivity")&&(constitutive_param_type!="rel_permeability")&&(constitutive_param_type!="electric_conductivity")&&(constitutive_param_type!="magnetic_conductivity"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,constitutive_param_type,
			"(valid arguments are \"rel_permittivity\", \"rel_permeability\", \"electric_conductivity\", and \"magnetic_conductivity\")");
	}

	ifstream MaterialFile;	//temporary ifstream object for reading the input
	MaterialFile.open(MaterialFileName.c_str(),ios::binary);	//open file for reading
	if (!MaterialFile)
	{
		cout << "Error opening material input file " << MaterialFileName << "." << endl << endl;
		exit(-1);
	}
	int xExtentInCells,yExtentInCells,zExtentInCells;	//x, y and z extents of the material region (in cells)
	MaterialFile.read((char*)&xExtentInCells,sizeof(xExtentInCells));	//read the x extent
	MaterialFile.read((char*)&yExtentInCells,sizeof(yExtentInCells));	//read the y extent
	MaterialFile.read((char*)&zExtentInCells,sizeof(zExtentInCells));	//read the z extent

	//read through the file to determine the maximum and minumum constitutive parameter values
	int pos_saved = MaterialFile.tellg();	//first, save the current position of the read pointer
	MatType max_param=0;	//maximum constitutive parameter value
	MatType min_param=1e10;	//minimum constitutive parameter value
	MatType param_temp;	//constitutive parameter value that has been read
	double param_lower_limit; //lower limit of the constitutive parameter
	if ((constitutive_param_type=="rel_permittivity")||(constitutive_param_type=="rel_permeability"))
	{
		param_lower_limit = 1;
	}
	else if ((constitutive_param_type=="electric_conductivity")||(constitutive_param_type=="magnetic_conductivity"))
	{
		param_lower_limit = 0;
	}
	else
	{
		throw AngoraDeveloperException("Error in PlaceMaterialRegionFromFile: unknown constitutive parameter type");
	}
	for (int i=1; i<=xExtentInCells; i++)
	{
		for (int j=1; j<=yExtentInCells; j++)
		{
			for (int k=1; k<=zExtentInCells; k++)
			{
				MaterialFile.read((char*)&param_temp,sizeof(param_temp));	//read the constitutive parameter
				if (param_temp>=param_lower_limit)
				{//if the value is nonpositive, don't bother with it at all
					if (param_temp>max_param) max_param=param_temp;	//update the maximum constitutive parameter
					if (param_temp<min_param) min_param=param_temp;	//update the minimum constitutive parameter
				}
			}
		}
	}
	//maximum number of different material types that can be extracted from the region
//	int max_num_of_materials = 1000; 	//pretty random, may have to find a more efficient way in the future
	//minimum difference in constitutive parameter between different materials
	MatType param_step = (max_param-min_param)/(max_number_of_new_materials-1);
	//add the new materials to the material list
	//before increasing the number of materials, save the current maximum material indices
	int material_index_saved_eps_x = eps_x.size()-1;
	int material_index_saved_eps_y = eps_y.size()-1;
	int material_index_saved_eps_z = eps_z.size()-1;
	int material_index_saved_mu_x = mu_x.size()-1;
	int material_index_saved_mu_y = mu_y.size()-1;
	int material_index_saved_mu_z = mu_z.size()-1;
	int material_index_saved_cond_e_x = cond_e_x.size()-1;
	int material_index_saved_cond_e_y = cond_e_y.size()-1;
	int material_index_saved_cond_e_z = cond_e_z.size()-1;
	int material_index_saved_cond_h_x = cond_h_x.size()-1;
	int material_index_saved_cond_h_y = cond_h_y.size()-1;
	int material_index_saved_cond_h_z = cond_h_z.size()-1;
	//dummy material index
	Cmat NewMaterial;
	if ((constitutive_param_type=="rel_permittivity"))
	{//the relative permittivity of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_eps(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="rel_permeability"))
	{//the relative permeability of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_mu(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="electric_conductivity"))
	{//the electric conductivity (in S/m) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_cond_e(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="magnetic_conductivity"))
	{//the magnetic conductivity (in Ohm/m) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_cond_h(min_param+(i-1)*param_step);
		}
	}

	MaterialFile.seekg(pos_saved,ios::beg);	//return to the saved position in the file

	//now, read through the file again, determine the material index for each point, and update the material indices in the main grid
	int material_offset;	//offset of the current material index beginning from the material index that was saved before the creation of new materials
	int material_index_eps,material_index_cond_e,material_index_mu,material_index_cond_h;		//material indices at a given point

	//calculate the coordinates of the back-left-lower corner of the region
	int xCornerPos=xPos;
	int yCornerPos=yPos;
	int zCornerPos=zPos;
	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
	{
		if (rank==0)
		{
			cout << "Invalid anchor point \"" << anchor << "\" for material input file " << MaterialFileName << " in node " << rank << endl << endl;
			exit(-1);
		}
	}

	if (anchor=="center")
	{
		xCornerPos = (int)floor(xPos-xExtentInCells/2.0);	//shift reference to the center
		yCornerPos = (int)floor(yPos-yExtentInCells/2.0);	//shift reference to the center
		zCornerPos = (int)floor(zPos-zExtentInCells/2.0);	//shift reference to the center
	}
	else if (anchor=="BLL")	//back-left-lower
	{
		//reference point is the back-left-lower corner by default
	}
	else if (anchor=="BLU")	//back-left-upper
	{
		zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
	}
	else if (anchor=="BRL")	//back-right-lower
	{
		yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
	}
	else if (anchor=="BRU")	//back-right-upper
	{
		yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
		zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
	}
	else if (anchor=="FLL")	//front-left-lower
	{
		xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
	}
	else if (anchor=="FLU")	//front-left-upper
	{
		xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
		zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
	}
	else if (anchor=="FRL")	//front-right-lower
	{
		xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
		yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
	}
	else if (anchor=="FRU")	//front-right-upper
	{
		xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
		yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
		zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
	}

	//indices of the cell at the back-left-lower corner
	//these are simply equal to the coordinates of the back-left-lower corner plus one
	int xCornerCell = xCornerPos+1;
	int yCornerCell = yCornerPos+1;
	int zCornerCell = zCornerPos+1;

	//every node has to read the file, and update the necessary portions of their grid
	Array<MatType,1> material_row(Range(xCornerCell,xCornerCell+xExtentInCells-1));	//temporary array of size xExtentInCells that will hold the constitutive parameter values belonging to a x-row in the 2-D material property array in the file
	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.

	if (constitutive_param_type=="rel_permittivity")
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtentInCells-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtentInCells-1; j++)
			{
				MaterialFile.read((char*)material_row.data(),xExtentInCells*sizeof(material_row(0)));	//read the constitutive parameter x-row into temporary array
				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-1,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_eps = material_index_saved_eps_x + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								eps_x_indices(i,j,k) = material_index_eps;
								Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
								Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-0.5,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_eps = material_index_saved_eps_y + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								eps_y_indices(i,j,k) = material_index_eps;
								Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
								Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-1,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_eps = material_index_saved_eps_z + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								eps_z_indices(i,j,k) = material_index_eps;
								Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
								Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}
			}
		}
	}

	if (constitutive_param_type=="electric_conductivity")
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtentInCells-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtentInCells-1; j++)
			{
				MaterialFile.read((char*)material_row.data(),xExtentInCells*sizeof(material_row(0)));	//read the constitutive parameter x-row into temporary array
				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-1,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_e = material_index_saved_cond_e_x + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								cond_e_x_indices(i,j,k) = material_index_cond_e;
								Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
								Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-0.5,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_e = material_index_saved_cond_e_y + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								cond_e_y_indices(i,j,k) = material_index_cond_e;
								Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
								Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-1,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_e = material_index_saved_cond_e_z + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								cond_e_z_indices(i,j,k) = material_index_cond_e;
								Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
								Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
							}
						}
					}
				}
			}
		}
	}

	if (constitutive_param_type=="rel_permeability")
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtentInCells-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtentInCells-1; j++)
			{
				MaterialFile.read((char*)material_row.data(),xExtentInCells*sizeof(material_row(0)));	//read the constitutive parameter x-row into temporary array
				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-0.5,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_mu = material_index_saved_mu_x + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								mu_x_indices(i,j,k) = material_index_mu;
								Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
								Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-1,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_mu = material_index_saved_mu_y + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								mu_y_indices(i,j,k) = material_index_mu;
								Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
								Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-0.5,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_mu = material_index_saved_mu_z + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								mu_z_indices(i,j,k) = material_index_mu;
								Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
								Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
							}
						}
					}
				}
			}
		}
	}

	if (constitutive_param_type=="magnetic_conductivity")
	{
		for (int k=zCornerCell; k<=zCornerCell+zExtentInCells-1; k++)
		{
			for (int j=yCornerCell; j<=yCornerCell+yExtentInCells-1; j++)
			{
				MaterialFile.read((char*)material_row.data(),xExtentInCells*sizeof(material_row(0)));	//read the constitutive parameter x-row into temporary array
				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-1,j-0.5,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_h = material_index_saved_cond_h_x + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material at the center of the cell
								cond_h_x_indices(i,j,k) = material_index_cond_h;
								Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
								Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper)&&(j>=jleft)&&(j<=jright+1))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-1,k-0.5))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_h = material_index_saved_cond_h_y + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the lower side of the cell
								cond_h_y_indices(i,j,k) = material_index_cond_h;
								Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
								Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
							}
						}
					}
				}

				if ((k>=klower)&&(k<=kupper+1)&&(j>=jleft)&&(j<=jright))
				{
					for (int i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtentInCells-1); i++)
					{
						param_temp = material_row(i);
						if (shape_mask->IsInside(i-0.5,j-0.5,k-1))
						{
							if (param_temp>=min_param)
							{//if the value is nonpositive, don't bother with it at all
								material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
								//material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
								material_index_cond_h = material_index_saved_cond_h_z + (material_offset + 1);	//this is the absolute index in the current material list
																									// +1 because of the range of material_offset above
								//place material on the left side of the cell
								cond_h_z_indices(i,j,k) = material_index_cond_h;
								Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
								Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
							}
						}
					}
				}
			}
		}
	}

	//finally, close the file
	MaterialFile.close();
//	if (rank==0)
//	{
//		cout << "Material region read." << endl << endl;
//	}
}

#endif // MATFILE_H
