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

//Defines functions that place materials with different geometries in the grid.

#include "headers.h"

#include "placegeom.h"

#include "material/Cmat.h"

#include <fstream>

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int rank;

extern const int PEC;
extern const int vacuum;

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

extern Array<eps_x_type,1> eps_x_indices_on_z_axis;
extern Array<eps_y_type,1> eps_y_indices_on_z_axis;
extern Array<eps_z_type,1> eps_z_indices_on_z_axis;
extern Array<mu_x_type,1> mu_x_indices_on_z_axis;
extern Array<mu_y_type,1> mu_y_indices_on_z_axis;
extern Array<mu_z_type,1> mu_z_indices_on_z_axis;
extern Array<cond_e_x_type,1> cond_e_x_indices_on_z_axis;
extern Array<cond_e_y_type,1> cond_e_y_indices_on_z_axis;
extern Array<cond_e_z_type,1> cond_e_z_indices_on_z_axis;
extern Array<cond_h_x_type,1> cond_h_x_indices_on_z_axis;
extern Array<cond_h_y_type,1> cond_h_y_indices_on_z_axis;
extern Array<cond_h_z_type,1> cond_h_z_indices_on_z_axis;

extern int number_of_layers;
extern Array<float,1> eps_x_values_in_layers;
extern Array<float,1> eps_y_values_in_layers;
extern Array<float,1> eps_z_values_in_layers;
extern Array<float,1> mu_x_values_in_layers;
extern Array<float,1> mu_y_values_in_layers;
extern Array<float,1> mu_z_values_in_layers;
extern Array<float,1> cond_e_x_values_in_layers;
extern Array<float,1> cond_e_y_values_in_layers;
extern Array<float,1> cond_e_z_values_in_layers;
extern Array<float,1> cond_h_x_values_in_layers;
extern Array<float,1> cond_h_y_values_in_layers;
extern Array<float,1> cond_h_z_values_in_layers;
extern Array<int,1> LayerLowerZIndices;
extern Array<int,1> LayerThicknesses;
extern Array<bool,1> IsLayerGrounded;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;

/** should be removed eventually **/
extern double epsilon_r_max_x,epsilon_r_min_x,epsilon_r_max_y,epsilon_r_min_y,epsilon_r_max_z,epsilon_r_min_z;
extern double mu_r_max_x,mu_r_min_x,mu_r_max_y,mu_r_min_y,mu_r_max_z,mu_r_min_z;
extern double epsilon_r_max,epsilon_r_min;
extern double mu_r_max,mu_r_min;
extern double epsilon_r_upper,mu_r_upper,cond_e_upper,cond_h_upper,epsilon_r_lower,mu_r_lower,cond_e_lower,cond_h_lower;
extern double c_upper,c_lower;
bool MaterialIsSame(const eps_x_type&, const eps_y_type&, const eps_z_type&,
				   const mu_x_type&, const mu_y_type&, const mu_z_type&,
				   const cond_e_x_type&, const cond_e_y_type&, const cond_e_z_type&,
				   const cond_h_x_type&, const cond_h_y_type&, const cond_h_z_type&,
				   const eps_x_type&, const eps_y_type&, const eps_z_type&,
				   const mu_x_type&, const mu_y_type&, const mu_z_type&,
				   const cond_e_x_type&, const cond_e_y_type&, const cond_e_z_type&,
				   const cond_h_x_type&, const cond_h_y_type&, const cond_h_z_type&);
/** should be removed eventually **/


void PlaceSlab(const Cmat& material, const double& start_position, const double& end_position)
//Places an infinite material slab at the specified position
//start_position is the z coordinate (w.r.t. the origin, in grid cells) of the lower interface
//end_position is the z coordinate (w.r.t. the origin, in grid cells) of the upper interface
//Both start_position and end_position can be non-integer. Effective permittivities and permeabilities are used at interfaces for second order accuracy.
// The formulas for the permittivity and permeability are from
// K.-P. Hwang and A. C. Cangellaris, “Effective permittivities for second-order accurate FDTD equations at dielectric interfaces,” IEEE Microw. Wireless Compon. Lett., vol. 11, no. 4, pp. 158–60, Apr. 2001.
// The arithmetic averaging of the electric and magnetic conductivities is more or less empirical.
{
	//Define the "uppermost" and "lowermost" grid voxel indices
	int SlabLower = (int)round(start_position)+1;//cell index is 1 more than the z coordinate of the x-y interface layer
	int SlabUpper = (int)round(end_position);//cell index is equal to the z coordinate of the x-y interface layer

	/** These are not used yet!!! start_position, end_position assumed integer! **/
	//the distance (in cell units) of the interfaces to the closest x-y grid interface layer
	double d_upper = end_position-SlabUpper;
	double d_lower = start_position-(SlabLower-1);
	int tang_comp_upper_cell_index = SlabUpper+1;
	int norm_comp_upper_cell_index = (d_upper>=0?SlabUpper+1:SlabUpper);
	int tang_comp_lower_cell_index = SlabLower;
	int norm_comp_lower_cell_index = (d_lower>=0?SlabLower:SlabLower-1);
	/** These are not used yet!!! start_position, end_position assumed integer! **/

	if (SlabUpper<SlabLower)
	{
//		if (rank==0)
//		{
//			cout << "Error: Uppermost cell position (" << SlabUpper << ") in material slab is smaller than the lowermost cell position (" << SlabLower << ")." << endl;
//		}
//		exit(-1);
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//			InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<int>(func_name,SlabUpper,"(thickness of slab should be positive)");
	}

	{
	int i,j,k; //counters
	/** Place the bulk of the slab **/
	//full-integer positions
	for (k=max(klower,SlabLower); k<=min(kupper+1,SlabUpper+1); k++)
	{
		for (i=Ca_X.lbound(firstDim);i<=Ca_X.ubound(firstDim);i++) {
			for (j=Ca_X.lbound(secondDim);j<=Ca_X.ubound(secondDim);j++)
			{
				if (material.eps_x_exists()) eps_x_indices(i,j,k)=material.eps_x_index();
				if (material.cond_e_x_exists()) cond_e_x_indices(i,j,k)=material.cond_e_x_index();
				Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
				Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
			}
		}

		for (i=Ca_Y.lbound(firstDim);i<=Ca_Y.ubound(firstDim);i++) {
			for (j=Ca_Y.lbound(secondDim);j<=Ca_Y.ubound(secondDim);j++)
			{
				if (material.eps_y_exists()) eps_y_indices(i,j,k)=material.eps_y_index();
				if (material.cond_e_y_exists()) cond_e_y_indices(i,j,k)=material.cond_e_y_index();
				Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
				Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
			}
		}

		for (i=Da_Z.lbound(firstDim);i<=Da_Z.ubound(firstDim);i++) {
			for (j=Da_Z.lbound(secondDim);j<=Da_Z.ubound(secondDim);j++)
			{
				if (material.mu_z_exists()) mu_z_indices(i,j,k)=material.mu_z_index();
				if (material.cond_h_z_exists()) cond_h_z_indices(i,j,k)=material.cond_h_z_index();
				Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
				Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
			}
		}
	}
	//half-integer positions
	for (k=max(klower,SlabLower); k<=min(kupper,SlabUpper); k++)
	{
		for (i=Da_X.lbound(firstDim);i<=Da_X.ubound(firstDim);i++) {
			for (j=Da_X.lbound(secondDim);j<=Da_X.ubound(secondDim);j++)
			{
				if (material.mu_x_exists()) mu_x_indices(i,j,k)=material.mu_x_index();
				if (material.cond_h_x_exists()) cond_h_x_indices(i,j,k)=material.cond_h_x_index();
				Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
				Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
			}
		}

		for (i=Da_Y.lbound(firstDim);i<=Da_Y.ubound(firstDim);i++) {
			for (j=Da_Y.lbound(secondDim);j<=Da_Y.ubound(secondDim);j++)
			{
				if (material.mu_y_exists()) mu_y_indices(i,j,k)=material.mu_y_index();
				if (material.cond_h_y_exists()) cond_h_y_indices(i,j,k)=material.cond_h_y_index();
				Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
				Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
			}
		}

		for (i=Ca_Z.lbound(firstDim);i<=Ca_Z.ubound(firstDim);i++) {
			for (j=Ca_Z.lbound(secondDim);j<=Ca_Z.ubound(secondDim);j++)
			{
				if (material.eps_z_exists()) eps_z_indices(i,j,k)=material.eps_z_index();
				if (material.cond_e_z_exists()) cond_e_z_indices(i,j,k)=material.cond_e_z_index();
				Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
				Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
			}
		}
	}
	}//loop-variable protection block

	/** Place interface layers **/
	//create materials for the two planar interfaces with averaged constitutive parameters
	//Material indices for upper and lower interfaces
	Cmat upper_interf_material,lower_interf_material;
	/** Place upper interface layer **/
	if ((SlabUpper<NCELLS_Z+2*NPML)&&(SlabUpper>=0))	//if the uppermost cell is outside the upper edge of the grid
														//then this is a half space, so do not place upper interface
														// (also check if the uppermost cell is beyond the lower edge)
	{
		//second permittivities in the permittivity averages are taken from the next higher z position
		if (material.eps_x_exists()) upper_interf_material.set_eps_x((material.eps_x_value() + eps_x(eps_x_indices_on_z_axis(SlabUpper+2)))/2);
		if (material.eps_y_exists()) upper_interf_material.set_eps_y((material.eps_y_value() + eps_y(eps_y_indices_on_z_axis(SlabUpper+2)))/2);
		//second permeability in both the denominator and numerator below are taken from the next higher z position
		if (material.mu_z_exists()) upper_interf_material.set_mu_z(2*material.mu_z_value()*mu_z(mu_z_indices_on_z_axis(SlabUpper+2))/(material.mu_z_value()+mu_z(mu_z_indices_on_z_axis(SlabUpper+2))));
		//the following conductivity averages are empirical
		if (material.cond_e_x_exists()) upper_interf_material.set_cond_e_x((material.cond_e_x_value() + cond_e_x(cond_e_x_indices_on_z_axis(SlabUpper+2)))/2); //this is empirical
		if (material.cond_e_y_exists()) upper_interf_material.set_cond_e_y((material.cond_e_y_value() + cond_e_y(cond_e_y_indices_on_z_axis(SlabUpper+2)))/2); //this is empirical
		if (material.cond_h_z_exists()) upper_interf_material.set_cond_h_z((material.cond_h_z_value() + cond_h_z(cond_h_z_indices_on_z_axis(SlabUpper+2)))/2); //this is empirical

		{
		int i,j,k; //counters
		//Then, place the material index pointers
		if ((klower<=SlabUpper+1)&&(kupper>=SlabUpper)) //don't do anything if this interface layer does not belong to this node
		{
			for (i=Ca_X.lbound(firstDim);i<=Ca_X.ubound(firstDim);i++) {
				for (j=Ca_X.lbound(secondDim);j<=Ca_X.ubound(secondDim);j++)
				{
					if (upper_interf_material.eps_x_exists()) eps_x_indices(i,j,SlabUpper+1)=upper_interf_material.eps_x_index();
					if (upper_interf_material.cond_e_x_exists()) cond_e_x_indices(i,j,SlabUpper+1)=upper_interf_material.cond_e_x_index();
					Ca_X(i,j,SlabUpper+1)=(1-dt*cond_e_x(cond_e_x_indices(i,j,SlabUpper+1))/(2.0*eps_x(eps_x_indices(i,j,SlabUpper+1))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,SlabUpper+1))/(2.0*eps_x(eps_x_indices(i,j,SlabUpper+1))*epsilon_0));
					Cb_X(i,j,SlabUpper+1)=dt/eps_x(eps_x_indices(i,j,SlabUpper+1))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,SlabUpper+1))/(2.0*eps_x(eps_x_indices(i,j,SlabUpper+1))*epsilon_0));
				}
			}
			for (i=Ca_Y.lbound(firstDim);i<=Ca_Y.ubound(firstDim);i++) {
				for (j=Ca_Y.lbound(secondDim);j<=Ca_Y.ubound(secondDim);j++)
				{
					if (upper_interf_material.eps_y_exists()) eps_y_indices(i,j,SlabUpper+1)=upper_interf_material.eps_y_index();
					if (upper_interf_material.cond_e_y_exists()) cond_e_y_indices(i,j,SlabUpper+1)=upper_interf_material.cond_e_y_index();
					Ca_Y(i,j,SlabUpper+1)=(1-dt*cond_e_y(cond_e_y_indices(i,j,SlabUpper+1))/(2.0*eps_y(eps_y_indices(i,j,SlabUpper+1))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,SlabUpper+1))/(2.0*eps_y(eps_y_indices(i,j,SlabUpper+1))*epsilon_0));
					Cb_Y(i,j,SlabUpper+1)=dt/eps_y(eps_y_indices(i,j,SlabUpper+1))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,SlabUpper+1))/(2.0*eps_y(eps_y_indices(i,j,SlabUpper+1))*epsilon_0));
				}
			}
			for (i=Da_Z.lbound(firstDim);i<=Da_Z.ubound(firstDim);i++) {
				for (j=Da_Z.lbound(secondDim);j<=Da_Z.ubound(secondDim);j++)
				{
					if (upper_interf_material.mu_z_exists()) mu_z_indices(i,j,SlabUpper+1)=upper_interf_material.mu_z_index();
					if (upper_interf_material.cond_h_z_exists()) cond_h_z_indices(i,j,SlabUpper+1)=upper_interf_material.cond_h_z_index();
					Da_Z(i,j,SlabUpper+1)=(1-dt*cond_h_z(cond_h_z_indices(i,j,SlabUpper+1))/(2.0*mu_z(mu_z_indices(i,j,SlabUpper+1))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,SlabUpper+1))/(2.0*mu_z(mu_z_indices(i,j,SlabUpper+1))*mu_0));
					Db_Z(i,j,SlabUpper+1)=dt/mu_z(mu_z_indices(i,j,SlabUpper+1))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,SlabUpper+1))/(2.0*mu_z(mu_z_indices(i,j,SlabUpper+1))*mu_0));
				}
			}
		}
		}//loop-variable protection block
	}

	/** Place lower interface layer **/
	if ((SlabLower>1)&&(SlabLower<=NCELLS_Z+2*NPML+1))	//if the lowermost cell is outside the lower edge of the grid
														//then this is a half space, so do not place lower interface
														// (also check if the lowermost cell is beyond the upper edge)
	{
		//second permittivities in the permittivity averages are taken from the next lower z position
		if (material.eps_x_exists()) lower_interf_material.set_eps_x((material.eps_x_value() + eps_x(eps_x_indices_on_z_axis(SlabLower-1)))/2);
		if (material.eps_y_exists()) lower_interf_material.set_eps_y((material.eps_y_value() + eps_y(eps_y_indices_on_z_axis(SlabLower-1)))/2);
		//second permeability in both the denominator and numerator below are taken from the next lower z position
		if (material.mu_z_exists()) lower_interf_material.set_mu_z(2*material.mu_z_value()*mu_z(mu_z_indices_on_z_axis(SlabLower-1))/(material.mu_z_value()+mu_z(mu_z_indices_on_z_axis(SlabLower-1))));
		//the following conductivity averages are empirical
		if (material.cond_e_x_exists()) lower_interf_material.set_cond_e_x((material.cond_e_x_value() + cond_e_x(cond_e_x_indices_on_z_axis(SlabLower-1)))/2);
		if (material.cond_e_y_exists()) lower_interf_material.set_cond_e_y((material.cond_e_y_value() + cond_e_y(cond_e_y_indices_on_z_axis(SlabLower-1)))/2);
		if (material.cond_h_z_exists()) lower_interf_material.set_cond_h_z((material.cond_h_z_value() + cond_h_z(cond_h_z_indices_on_z_axis(SlabLower-1)))/2);

		{
		int i,j,k; //counters
		//Then, place the material index pointers
		if ((klower<=SlabLower)&&(kupper>=SlabLower-1)) //don't do anything if this interface layer does not belong to this node
		{
			for (i=Ca_X.lbound(firstDim);i<=Ca_X.ubound(firstDim);i++) {
				for (j=Ca_X.lbound(secondDim);j<=Ca_X.ubound(secondDim);j++)
				{
					if (lower_interf_material.eps_x_exists()) eps_x_indices(i,j,SlabLower)=lower_interf_material.eps_x_index();
					if (lower_interf_material.cond_e_x_exists()) cond_e_x_indices(i,j,SlabLower)=lower_interf_material.cond_e_x_index();
					Ca_X(i,j,SlabLower)=(1-dt*cond_e_x(cond_e_x_indices(i,j,SlabLower))/(2.0*eps_x(eps_x_indices(i,j,SlabLower))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,SlabLower))/(2.0*eps_x(eps_x_indices(i,j,SlabLower))*epsilon_0));
					Cb_X(i,j,SlabLower)=dt/eps_x(eps_x_indices(i,j,SlabLower))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,SlabLower))/(2.0*eps_x(eps_x_indices(i,j,SlabLower))*epsilon_0));
				}
			}

			for (i=Ca_Y.lbound(firstDim);i<=Ca_Y.ubound(firstDim);i++) {
				for (j=Ca_Y.lbound(secondDim);j<=Ca_Y.ubound(secondDim);j++)
				{
					if (lower_interf_material.eps_y_exists()) eps_y_indices(i,j,SlabLower)=lower_interf_material.eps_y_index();
					if (lower_interf_material.cond_e_y_exists()) cond_e_y_indices(i,j,SlabLower)=lower_interf_material.cond_e_y_index();
					Ca_Y(i,j,SlabLower)=(1-dt*cond_e_y(cond_e_y_indices(i,j,SlabLower))/(2.0*eps_y(eps_y_indices(i,j,SlabLower))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,SlabLower))/(2.0*eps_y(eps_y_indices(i,j,SlabLower))*epsilon_0));
					Cb_Y(i,j,SlabLower)=dt/eps_y(eps_y_indices(i,j,SlabLower))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,SlabLower))/(2.0*eps_y(eps_y_indices(i,j,SlabLower))*epsilon_0));
				}
			}

			for (i=Da_Z.lbound(firstDim);i<=Da_Z.ubound(firstDim);i++) {
				for (j=Da_Z.lbound(secondDim);j<=Da_Z.ubound(secondDim);j++)
				{
					if (lower_interf_material.mu_z_exists()) mu_z_indices(i,j,SlabLower)=lower_interf_material.mu_z_index();
					if (lower_interf_material.cond_h_z_exists()) cond_h_z_indices(i,j,SlabLower)=lower_interf_material.cond_h_z_index();
					Da_Z(i,j,SlabLower)=(1-dt*cond_h_z(cond_h_z_indices(i,j,SlabLower))/(2.0*mu_z(mu_z_indices(i,j,SlabLower))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,SlabLower))/(2.0*mu_z(mu_z_indices(i,j,SlabLower))*mu_0));
					Db_Z(i,j,SlabLower)=dt/mu_z(mu_z_indices(i,j,SlabLower))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,SlabLower))/(2.0*mu_z(mu_z_indices(i,j,SlabLower))*mu_0));
				}
			}
		}
		}//loop-variable protection block
	}

	/** Update layering info **/
	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
	//x and y components
	for (int k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML+1,SlabUpper+1); k++)
	{
		if (material.eps_x_exists()) eps_x_indices_on_z_axis(k)=material.eps_x_index();
		if (material.eps_y_exists()) eps_y_indices_on_z_axis(k)=material.eps_y_index();
		if (material.mu_z_exists()) mu_z_indices_on_z_axis(k)=material.mu_z_index();
		if (material.cond_e_x_exists()) cond_e_x_indices_on_z_axis(k)=material.cond_e_x_index();
		if (material.cond_e_y_exists()) cond_e_y_indices_on_z_axis(k)=material.cond_e_y_index();
		if (material.cond_h_z_exists()) cond_h_z_indices_on_z_axis(k)=material.cond_h_z_index();
	}
	//z component
	for (int k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML,SlabUpper); k++)
	{
		if (material.eps_z_exists()) eps_z_indices_on_z_axis(k)=material.eps_z_index();
		if (material.mu_x_exists()) mu_x_indices_on_z_axis(k)=material.mu_x_index();
		if (material.mu_y_exists()) mu_y_indices_on_z_axis(k)=material.mu_y_index();
		if (material.cond_e_z_exists()) cond_e_z_indices_on_z_axis(k)=material.cond_e_z_index();
		if (material.cond_h_x_exists()) cond_h_x_indices_on_z_axis(k)=material.cond_h_x_index();
		if (material.cond_h_y_exists()) cond_h_y_indices_on_z_axis(k)=material.cond_h_y_index();
	}
	if ((SlabLower>1)&&(SlabLower<=NCELLS_Z+2*NPML+1))
	{
		if (lower_interf_material.eps_x_exists()) eps_x_indices_on_z_axis(SlabLower)=lower_interf_material.eps_x_index();
		if (lower_interf_material.eps_y_exists()) eps_y_indices_on_z_axis(SlabLower)=lower_interf_material.eps_y_index();
		if (lower_interf_material.mu_z_exists()) mu_z_indices_on_z_axis(SlabLower)=lower_interf_material.mu_z_index();
		if (lower_interf_material.cond_e_x_exists()) cond_e_x_indices_on_z_axis(SlabLower)=lower_interf_material.cond_e_x_index();
		if (lower_interf_material.cond_e_y_exists()) cond_e_y_indices_on_z_axis(SlabLower)=lower_interf_material.cond_e_y_index();
		if (lower_interf_material.cond_h_z_exists()) cond_h_z_indices_on_z_axis(SlabLower)=lower_interf_material.cond_h_z_index();
	}
	if ((SlabUpper<NCELLS_Z+2*NPML)&&(SlabUpper>=0))
	{
		if (upper_interf_material.eps_x_exists()) eps_x_indices_on_z_axis(SlabUpper+1)=upper_interf_material.eps_x_index();
		if (upper_interf_material.eps_y_exists()) eps_y_indices_on_z_axis(SlabUpper+1)=upper_interf_material.eps_y_index();
		if (upper_interf_material.mu_z_exists()) mu_z_indices_on_z_axis(SlabUpper+1)=upper_interf_material.mu_z_index();
		if (upper_interf_material.cond_e_x_exists()) cond_e_x_indices_on_z_axis(SlabUpper+1)=upper_interf_material.cond_e_x_index();
		if (upper_interf_material.cond_e_y_exists()) cond_e_y_indices_on_z_axis(SlabUpper+1)=upper_interf_material.cond_e_y_index();
		if (upper_interf_material.cond_h_z_exists()) cond_h_z_indices_on_z_axis(SlabUpper+1)=upper_interf_material.cond_h_z_index();
	}
}

void PlaceGround(const int& GroundCoord)
//Places an infinite PEC ground plane at the specified position
{
	//Define the "uppermost" and "lowermost" grid voxel indices
	int GroundCellIndex = GroundCoord+1;  //index of the cell right above the ground layer
//	//Define the "uppermost" and "lowermost" grid voxel indices
//	int GroundCellIndex = (int)round(GroundCoord)+1;  //index of the cell right above the ground layer

//	/** This is not used yet!!! GroundCoord assumed integer! **/
//	//the distance (in cell units) of the PEC ground layer to the next higher closest x-y grid grid interface
//	double d = GroundCoord-GroundCellIndex;
//	/** This is not used yet!!! GroundCoord assumed integer! **/

	{
	int i,j; //counters
	//Place PEC ground
	if ((klower<=GroundCellIndex)&&(kupper>=GroundCellIndex-1))
	{
		for (i=Ca_X.lbound(firstDim);i<=Ca_X.ubound(firstDim);i++) {
			for (j=Ca_X.lbound(secondDim);j<=Ca_X.ubound(secondDim);j++)
			{
				eps_x_indices(i,j,GroundCellIndex)=PEC;
				cond_e_x_indices(i,j,GroundCellIndex)=PEC;
				Ca_X(i,j,GroundCellIndex)=1.0;
				Cb_X(i,j,GroundCellIndex)=0.0;
			}
		}

		for (i=Ca_Y.lbound(firstDim);i<=Ca_Y.ubound(firstDim);i++) {
			for (j=Ca_Y.lbound(secondDim);j<=Ca_Y.ubound(secondDim);j++)
			{
				eps_y_indices(i,j,GroundCellIndex)=PEC;
				cond_e_y_indices(i,j,GroundCellIndex)=PEC;
				Ca_Y(i,j,GroundCellIndex)=1.0;
				Cb_Y(i,j,GroundCellIndex)=0.0;
			}
		}
	}
	}//loop-variable protection block

	//Update layering info
	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
	if ((GroundCellIndex>=1)&&(GroundCellIndex<=NCELLS_Z+2*NPML+1))
	{
		eps_x_indices_on_z_axis(GroundCellIndex)=PEC;
		cond_e_x_indices_on_z_axis(GroundCellIndex)=PEC;

		eps_y_indices_on_z_axis(GroundCellIndex)=PEC;
		cond_e_y_indices_on_z_axis(GroundCellIndex)=PEC;
	}
}

//void PlacePECSlab(const double& start_position, const double& end_position)
////Places an infinite perfect-electric-conductor (PEC) slab at the specified position
//{
//	int k; //counter for looping through z positions
//
//	//Define the "uppermost" and "lowermost" grid voxel indices
//	int SlabLower = (int)round(start_position)+1;//cell index is 1 more than the z coordinate of the x-y interface layer
//	int SlabUpper = (int)round(end_position);//cell index is equal to the z coordinate of the x-y interface layer
//
//	/** These are not used yet!!! start_position, end_position assumed integer! **/
//	//the distance (in cell units) of the interfaces to the closest x-y grid interface layer
//	double d_lower = start_position-SlabLower;
//	double d_upper = end_position-SlabUpper;
//	/** These are not used yet!!! start_position, end_position assumed integer! **/
//
//	//Update eps_x,eps_y
//	for (k=max(klower,SlabLower); k<=min(kupper+1,SlabUpper+1); k++)
//	{
//		eps_x_indices(Range::all(),Range::all(),k)=PEC;
//		eps_y_indices(Range::all(),Range::all(),k)=PEC;
//		cond_e_x_indices(Range::all(),Range::all(),k)=PEC;
//		cond_e_y_indices(Range::all(),Range::all(),k)=PEC;
///** TODO : **/
//Ca_X(i,j,k)=(1-dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0))/(1+dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0));
//Cb_X(i,j,k)=dt/eps_x_indices(i,j,k)/epsilon_0/dx/(1+dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0));
///** TODO : **/
////		//x component
////		Media_Ex(Range::all(),Range::all(),k)=PEC;
////		//y component
////		Media_Ey(Range::all(),Range::all(),k)=PEC;
//	}
//	//Update eps_z
//	for (k=max(klower,SlabLower); k<=min(kupper,SlabUpper); k++)
//	{
//		eps_z_indices(Range::all(),Range::all(),k)=PEC;
//		cond_e_z_indices(Range::all(),Range::all(),k)=PEC;
///** TODO : **/
//Ca_X(i,j,k)=(1-dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0))/(1+dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0));
//Cb_X(i,j,k)=dt/eps_x_indices(i,j,k)/epsilon_0/dx/(1+dt*cond_e_x_indices(i,j,k)/(2.0*eps_x_indices(i,j,k)*epsilon_0));
///** TODO : **/
////		//z component
////		Media_Ez(Range::all(),Range::all(),k)=PEC;
//	}
//
//	//Update layering info
//	//(Note that this is done whether the slab passes through the node or not! The layering info must be known by all nodes!)
//	//x and y components
//	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML+1,SlabUpper+1); k++)
//	{
//		eps_x_indices_on_z_axis(k)=PEC;
//		eps_y_indices_on_z_axis(k)=PEC;
//		cond_e_x_indices_on_z_axis(k)=PEC;
//		cond_e_y_indices_on_z_axis(k)=PEC;
////		Layering_e_x(k)=PEC;
////		Layering_e_y(k)=PEC;
//	}
//	//z component
//	for (k=max(1,SlabLower); k<=min(NCELLS_Z+2*NPML,SlabUpper); k++)
//	{
//		eps_z_indices_on_z_axis(k)=PEC;
//		cond_e_z_indices_on_z_axis(k)=PEC;
////		Layering_e_z(k)=PEC;
//	}
//
//	//In PlaceSlab(), epsilon_r_upper, c_upper, etc. is updated at this point. (04/24/2011: Not anymore!) If the PEC slab extends to the upper or lower limit of the grid, it is the responsibility of the particular class that normally requires epsilon_r_upper, c_upper, etc. to be smart and accomodate the PEC layer at the top or bottom instead of using epsilon_r_upper, c_upper, etc.
//}

//void PlacePECBlock(
//		const int& BlockBack, const int& BlockFront,	//indices of the rearmost and foremost cells in the block
//		const int& BlockLeft, const int& BlockRight,	//indices of the leftmost and rightmost cells in the block
//		const int& BlockLower, const int& BlockUpper)	//indices of the lowermost and uppermost cells in the block
////Places a PEC block at the specified position
//{
//	int i,j,k;
//	//x component
//	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
//			{
//				Media_Ex(i,j,k)=PEC;
//			}
//		}
//	}
//	//y component
//	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
//			{
//				Media_Ey(i,j,k)=PEC;
//			}
//		}
//	}
//	//z component
//	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
//			{
//				Media_Ez(i,j,k)=PEC;
//			}
//		}
//	}
//}
//
//void PlaceMaterialBlock(
//			const MaterialId& MaterialIdentifier,
//			const int& BlockBack, const int& BlockFront,	//indices of the rearmost and foremost cells in the block
//			const int& BlockLeft, const int& BlockRight, 	//indices of the leftmost and rightmost cells in the block
//			const int& BlockLower, const int& BlockUpper)	//indices of the lowermost and uppermost cells in the block
////Places a square block made of material with MaterialIndex at the specified position
//{
//	int i,j,k;
//	//electric-field components
//	//x component
//	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
//			{
//				Media_Ex(i,j,k)=MaterialIdentifier.ElectricIndex_X;
//			}
//		}
//	}
//	//y component
//	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
//			{
//				Media_Ey(i,j,k)=MaterialIdentifier.ElectricIndex_Y;
//			}
//		}
//	}
//	//z component
//	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
//			{
//				Media_Ez(i,j,k)=MaterialIdentifier.ElectricIndex_Z;
//			}
//		}
//	}
//
//	//magnetic-field components
//	//x component
//	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
//			{
//				Media_Hx(i,j,k)=MaterialIdentifier.MagneticIndex_X;
//			}
//		}
//	}
//	//y component
//	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
//			{
//				Media_Hy(i,j,k)=MaterialIdentifier.MagneticIndex_Y;
//			}
//		}
//	}
//	//z component
//	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
//	{
//		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
//		{
//			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
//			{
//				Media_Hz(i,j,k)=MaterialIdentifier.MagneticIndex_Z;
//			}
//		}
//	}
//}
//
//void PlaceRotatedBlock(
//			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
//			const double& ThicknessX, const double& ThicknessY, const double& ThicknessZ,	//thicknesses
//			const double& phi_deg)			//rotation angle (in degrees)
////Places a rotated (in the xy-plane) PEC block at the specified position
//{
//	double x,y,z;	//x,y,z coordinates with respect to the center of the block
//	double x_dist,y_dist,z_dist;	//distances of the point to the orthogonal faces of the cube
//
//	double phi = phi_deg*M_PI/180;
//
//	int i,j,k;
//	//x component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j-CenterY);
//				z = (k-CenterZ);
//				x_dist = x*cos(phi)+y*sin(phi);
//				y_dist = -x*sin(phi)+y*cos(phi);
//				z_dist = z;
//				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
//				{
//					Media_Ex(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//y component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k-CenterZ);
//				x_dist = x*cos(phi)+y*sin(phi);
//				y_dist = -x*sin(phi)+y*cos(phi);
//				z_dist = z;
//				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
//				{
//					Media_Ey(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//z component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i-CenterX);
//				y = (j-CenterY);
//				z = (k+0.5-CenterZ);
//				x_dist = x*cos(phi)+y*sin(phi);
//				y_dist = -x*sin(phi)+y*cos(phi);
//				z_dist = z;
//				if ((abs(x_dist)<=ThicknessX/2.0)&&(abs(y_dist)<=ThicknessY/2.0)&&(abs(z_dist)<=ThicknessZ/2.0))
//				{
//					Media_Ez(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//}
//
//void PlaceCircularBlock(
//			const double& CenterX, const double& CenterY, const double& CenterZ,	//coordinates of the center point
//			const double& Radius, //radius of the bowtie (in cells)
//			const double& Height)	//height of the block (in cells)
////Places a circular (in the xy-plane) PEC block
//{
//	double x,y,z;	//x,y,z coordinates with respect to the center of the block
//	double dist_r;	//polar distance of the point to the center (in the xy plane)
//
//	int i,j,k;
//	//x component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2));
//				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
//				{
//					Media_Ex(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//y component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2));
//				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
//				{
//					Media_Ey(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//z component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i-CenterX);
//				y = (j-CenterY);
//				z = (k+0.5-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2));
//				if ((dist_r<=Radius)&&(abs(z)<=Height/2.0))
//				{
//					Media_Ez(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//}
//
//void PlaceSphere(const double& CenterX, const double& CenterY, const double& CenterZ, const double& Radius)
////Places a PEC sphere
//{
//	double x,y,z;	//x,y,z coordinates with respect to the center of the sphere
//	double dist_r;	//spherical distance of the point to the center of the sphere
//
//	int i,j,k;
//	//x component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ex(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//y component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ey(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//	//z component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i-CenterX);
//				y = (j-CenterY);
//				z = (k+0.5-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ez(i,j,k)=PEC;
//				}
//			}
//		}
//	}
//}
//
//void PlaceMaterialSphere(
//			const MaterialId& MaterialIdentifier,
//			const double& CenterX, const double& CenterY, const double& CenterZ,
//			const double& Radius)
//{//Places a sphere made of material with identifier MaterialIdentifier
//	double x,y,z;	//x,y,z coordinates with respect to the center of the sphere
//	double dist_r;	//spherical distance of the point to the center of the sphere
//
//	int i,j,k;
//	//electric-field components
//	//x component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ex(i,j,k)=MaterialIdentifier.ElectricIndex_X;
//				}
//			}
//		}
//	}
//	//y component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ey(i,j,k)=MaterialIdentifier.ElectricIndex_Y;
//				}
//			}
//		}
//	}
//	//z component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i-CenterX);
//				y = (j-CenterY);
//				z = (k+0.5-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Ez(i,j,k)=MaterialIdentifier.ElectricIndex_Z;
//				}
//			}
//		}
//	}
//
//	//magnetic-field components
//	//x component
//	for (i=iback; i<=ifront+1; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k+0.5-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Hx(i,j,k)=MaterialIdentifier.MagneticIndex_X;
//				}
//			}
//		}
//	}
//	//y component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright+1; j++)
//		{
//			for (k=klower; k<=kupper; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j-CenterY);
//				z = (k+0.5-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Hy(i,j,k)=MaterialIdentifier.MagneticIndex_Y;
//				}
//			}
//		}
//	}
//	//z component
//	for (i=iback; i<=ifront; i++)
//	{
//		for (j=jleft; j<=jright; j++)
//		{
//			for (k=klower; k<=kupper+1; k++)
//			{
//				x = (i+0.5-CenterX);
//				y = (j+0.5-CenterY);
//				z = (k-CenterZ);
//				dist_r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
//				if (dist_r<=Radius)
//				{
//					Media_Hz(i,j,k)=MaterialIdentifier.MagneticIndex_Z;
//				}
//			}
//		}
//	}
//}
//
//void PlacePECMaskFromFile(const string& PECMaskFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor)
//{//Reads a rectangular-prism-shaped boolean PEC mask from file and applies to grid
//// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the applied mask (measured in cells from the back-left-lower corner of the grid)
//	if (rank==0)
//	{
//		cout << endl << "Reading PEC mask from " << PECMaskFileName << " ..." << endl;
//	}
//	ifstream PECMaskFile;	//temporary ifstream object for reading the input
//	PECMaskFile.open(PECMaskFileName.c_str(),ios::binary);	//open file for reading
//	if (!PECMaskFile)
//	{
//		cout << "Error opening PEC mask input file " << PECMaskFileName << "." << endl << endl;
//		exit(-1);
//	}
//	int xExtent,yExtent,zExtent;	//x, y and z extents of the material region (in cells)
//	PECMaskFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
//	PECMaskFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
//	PECMaskFile.read((char*)&zExtent,sizeof(zExtent));	//read the z extent
//
//	//calculate the coordinates of the back-left-lower corner of the region
//	int xCornerPos=xPos;
//	int yCornerPos=yPos;
//	int zCornerPos=zPos;
//	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
//	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
//	{
//		if (rank==0)
//		{
//			cout << "Invalid anchor point \"" << anchor << "\" for PEC mask input file " << PECMaskFileName << " in node " << rank << endl << endl;
//			exit(-1);
//		}
//	}
//
//	if (anchor=="center")
//	{
//		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
//		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
//		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
//	}
//	else if (anchor=="BLL")	//back-left-lower
//	{
//		//reference point is the back-left-lower corner by default
//	}
//	else if (anchor=="BLU")	//back-left-upper
//	{
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="BRL")	//back-right-lower
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="BRU")	//back-right-upper
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FLL")	//front-left-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//	}
//	else if (anchor=="FLU")	//front-left-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FRL")	//front-right-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="FRU")	//front-right-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//
//	//indices of the cell at the back-left-lower corner
//	//these are simply equal to the coordinates of the back-left-lower corner plus one
//	int xCornerCell = xCornerPos+1;
//	int yCornerCell = yCornerPos+1;
//	int zCornerCell = zCornerPos+1;
//
//	//every node has to read the file, and update the necessary portions of their grid
//	Array<bool,1> PEC_mask_row(Range(xCornerCell,xCornerCell+xExtent-1));	//temporary array of size xExtent that will hold the boolean values belonging to a x-row in the 3-D PEC mask array in the file
//	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
//	int i,j,k;
//	for (k=zCornerCell; k<=zCornerCell+zExtent-1; k++)
//	{
//		for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
//		{
//			PECMaskFile.read((char*)PEC_mask_row.data(),xExtent*sizeof(PEC_mask_row(0)));	//read the permittivity x-row into temporary array
//			//Media_Ex
//			if ((k>=klower)&&(k<=kupper+1))
//			{
//				if ((j>=jleft)&&(j<=jright+1))
//				{
//					for (i=max(iback,xCornerCell); i<=min(ifront,xCornerCell+xExtent-1); i++)
//					{
//						if (PEC_mask_row(i)==true)
//						{
//							//mark this position as perfect electric conductor
//							Media_Ex(i,j,k) = PEC;
//						}
//					}
//				}
//			}
//
//			//Media_Ey
//			if ((k>=klower)&&(k<=kupper+1))
//			{
//				if ((j>=jleft)&&(j<=jright))
//				{
//					for (i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
//					{
//						if (PEC_mask_row(i)==true)
//						{
//							//mark this position as perfect electric conductor
//							Media_Ey(i,j,k) = PEC;
//						}
//					}
//				}
//			}
//
//			//Media_Ez
//			if ((k>=klower)&&(k<=kupper))
//			{
//				if ((j>=jleft)&&(j<=jright+1))
//				{
//					for (i=max(iback,xCornerCell); i<=min(ifront+1,xCornerCell+xExtent-1); i++)
//					{
//						if (PEC_mask_row(i)==true)
//						{
//							//mark this position as perfect electric conductor
//							Media_Ez(i,j,k) = PEC;
//						}
//					}
//				}
//			}
//		}
//	}
//
//	//finally, close the file
//	PECMaskFile.close();
//	if (rank==0)
//	{
//		cout << "PEC mask applied to grid." << endl << endl;
//	}
//}
//
//void PlaceSurfaceProfileFromFile(
//			const string& SurfaceProfileFileName,
//			const MaterialId& MaterialIdentifier,
//			const int& xPos, const int& yPos, const int& zPos, const string& anchor)
//{
////Reads a rectangular homogeneous region with surface roughness from a file.
////The region is completely specified by a 2D array indicating the surface fluctuations, the average thickness of the region, and the relative dielectric permittivity of the region (equal to square of the refractive index).
//// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the region (measured in cells from the back-left-lower corner of the grid)
//	if (rank==0)
//	{
//		cout << endl << "Reading surface profile from " << SurfaceProfileFileName << " ..." << endl;
//	}
//	ifstream SurfaceProfileFile;	//temporary ifstream object for reading the input
//	SurfaceProfileFile.open(SurfaceProfileFileName.c_str(),ios::binary);	//open file for reading
//	if (!SurfaceProfileFile)
//	{
//		cout << "Error opening surface profile input file " << SurfaceProfileFileName << "." << endl << endl;
//		exit(-1);
//	}
//	int xExtent,yExtent,zExtent;	//x, y and average z extents of the material region (in cells)
//	SurfaceProfileFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
//	SurfaceProfileFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
//	SurfaceProfileFile.read((char*)&zExtent,sizeof(zExtent));	//read the average z extent
//
//	//Following removed after ver.>= 0.12
///*	double eps_r; //relative dielectric permittivity of the homogeneous region
//	//read the relative dielectric permittivity of the homogeneous region
//	SurfaceProfileFile.read((char*)&eps_r,sizeof(eps_r));*/
//// 	//add the new material to the material-index array
//// 	int surfacematerialindex = AddMaterial(eps_r,0);	//relative permittivity of the new material is eps_r
//
//	//calculate the coordinates of the back-left-lower corner of the region
//	int xCornerPos=xPos;
//	int yCornerPos=yPos;
//	int zCornerPos=zPos;
//	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
//	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
//	{
//		if (rank==0)
//		{
//			cout << "Invalid anchor point \"" << anchor << "\" for surface-profile input file " << SurfaceProfileFileName << " in node " << rank << endl << endl;
//			exit(-1);
//		}
//	}
//
//	if (anchor=="center")
//	{
//		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
//		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
//		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
//	}
//	else if (anchor=="BLL")	//back-left-lower
//	{
//		//reference point is the back-left-lower corner by default
//	}
//	else if (anchor=="BLU")	//back-left-upper
//	{
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="BRL")	//back-right-lower
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="BRU")	//back-right-upper
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FLL")	//front-left-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//	}
//	else if (anchor=="FLU")	//front-left-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FRL")	//front-right-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="FRU")	//front-right-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//
//	//indices of the cell at the back-left-lower corner
//	//these are simply equal to the coordinates of the back-left-lower corner plus one
//	int xCornerCell = xCornerPos+1;
//	int yCornerCell = yCornerPos+1;
//	int zCornerCell = zCornerPos+1;
//
//	//extra and absolute thickness values at each point of the surface
//	int extra_thickness,thickness; //in cells
//
//	//every node has to read the file, and update the necessary portions of their grid
//	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
//	int i,j,k;
//	for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
//	{
//		for (i=xCornerCell; i<=xCornerCell+xExtent-1; i++)
//		{
//			SurfaceProfileFile.read((char*)&extra_thickness,sizeof(extra_thickness));	//read the surface-profile at the (x,y) position
//			//calculate the thickness at this (x,y) position by adding the profile value to the average thickness (zExtent)
//			//(extra_thickness might also be negative)
//			thickness = zExtent+extra_thickness;
//			for (k=zCornerCell; k<=zCornerCell+thickness-1; k++)
//			{
//				//Media_Ex
//				if ((k>=klower)&&(k<=kupper+1))
//				{
//					if ((j>=jleft)&&(j<=jright+1))
//					{
//						if ((i>=iback)&&(i<=ifront))
//						{
//								Media_Ex(i,j,k) = MaterialIdentifier.ElectricIndex_X;
//						}
//					}
//				}
//
//				//Media_Ey
//				if ((k>=klower)&&(k<=kupper+1))
//				{
//					if ((j>=jleft)&&(j<=jright))
//					{
//						if ((i>=iback)&&(i<=ifront+1))
//						{
//								Media_Ey(i,j,k) = MaterialIdentifier.ElectricIndex_Y;
//						}
//					}
//				}
//
//				//Media_Ez
//				if ((k>=klower)&&(k<=kupper))
//				{
//					if ((j>=jleft)&&(j<=jright+1))
//					{
//						if ((i>=iback)&&(i<=ifront+1))
//						{
//								Media_Ez(i,j,k) = MaterialIdentifier.ElectricIndex_Z;
//						}
//					}
//				}
//			}
//		}
//	}
//
//	//finally, close the file
//	SurfaceProfileFile.close();
//	if (rank==0)
//	{
//		cout << "Surface profile read." << endl << endl;
//	}
//}
//
//void PlaceSurfaceEngravingProfileFromFile(
//			const string& SurfaceEngravingProfileFileName,
//			const MaterialId& MaterialIdentifier,
//			const int& xPos, const int& yPos, const int& zPos, const string& anchor)
//{
////Reads a 2D array that specifies an engraving profile. The surface is engraved at the specified height values in the array (in grid cells), and filled with the material specified by "MaterialIndex".
//// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the base of the engraving profile (measured in cells from the back-left-lower corner of the grid)
//	if (rank==0)
//	{
//		cout << endl << "Reading surface-engraving profile from " << SurfaceEngravingProfileFileName << " ..." << endl;
//	}
//	ifstream SurfaceEngravingProfileFile;	//temporary ifstream object for reading the input
//	SurfaceEngravingProfileFile.open(SurfaceEngravingProfileFileName.c_str(),ios::binary);	//open file for reading
//	if (!SurfaceEngravingProfileFile)
//	{
//		cout << "Error opening surface-engraving profile input file " << SurfaceEngravingProfileFileName << "." << endl << endl;
//		exit(-1);
//	}
//	int xExtent,yExtent,zExtent;	//x, y and average z extents of the material region (in cells)
//	SurfaceEngravingProfileFile.read((char*)&xExtent,sizeof(xExtent));	//read the x extent
//	SurfaceEngravingProfileFile.read((char*)&yExtent,sizeof(yExtent));	//read the y extent
//	SurfaceEngravingProfileFile.read((char*)&zExtent,sizeof(zExtent));	//read the average z extent
//
///*	double eps_r; //relative dielectric permittivity of the homogeneous region
//	//read the relative dielectric permittivity of the homogeneous region
//	SurfaceEngravingProfileFile.read((char*)&eps_r,sizeof(eps_r));*/
//
//// 	//add the new material to the material-index array
//// 	int surfacematerialindex = AddMaterial(eps_r,0);	//relative permittivity of the new material is eps_r
//
//	//calculate the coordinates of the back-left-lower corner of the region
//	int xCornerPos=xPos;
//	int yCornerPos=yPos;
//	int zCornerPos=zPos;
//	if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
//	   					  &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
//	{
//		if (rank==0)
//		{
//			cout << "Invalid anchor point \"" << anchor << "\" for surface-engraving profile input file " << SurfaceEngravingProfileFileName << " in node " << rank << endl << endl;
//			exit(-1);
//		}
//	}
//
//	if (anchor=="center")
//	{
//		xCornerPos = int(xPos - xExtent/2.0);	//shift reference to the center
//		yCornerPos = int(yPos - yExtent/2.0);	//shift reference to the center
//		zCornerPos = int(zPos - zExtent/2.0);	//shift reference to the center
//	}
//	else if (anchor=="BLL")	//back-left-lower
//	{
//		//reference point is the back-left-lower corner by default
//	}
//	else if (anchor=="BLU")	//back-left-upper
//	{
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="BRL")	//back-right-lower
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="BRU")	//back-right-upper
//	{
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FLL")	//front-left-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//	}
//	else if (anchor=="FLU")	//front-left-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//	else if (anchor=="FRL")	//front-right-lower
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//	}
//	else if (anchor=="FRU")	//front-right-upper
//	{
//		xCornerPos = xPos - xExtent;	//shift reference to the front corner
//		yCornerPos = yPos - yExtent;	//shift reference to the right corner
//		zCornerPos = zPos - zExtent;	//shift reference to the upper corner
//	}
//
//	//indices of the cell at the back-left-lower corner
//	//these are simply equal to the coordinates of the back-left-lower corner plus one
//	int xCornerCell = xCornerPos+1;
//	int yCornerCell = yCornerPos+1;
//	int zCornerCell = zCornerPos+1;
//
//	//extra and absolute thickness values at each point of the surface
//	int extra_thickness,thickness; //in cells
//// cout << xCornerCell << "," << yCornerCell << "," << zCornerCell << endl;
//	//every node has to read the file, and update the necessary portions of their grid
//	// Note that the x-y-z dimensions are written in this order in the file. Therefore, the x-rows are read first.
//	int i,j,k;
//	for (j=yCornerCell; j<=yCornerCell+yExtent-1; j++)
//	{
//		for (i=xCornerCell; i<=xCornerCell+xExtent-1; i++)
//		{
//			SurfaceEngravingProfileFile.read((char*)&extra_thickness,sizeof(extra_thickness));	//read the surface-profile at the (x,y) position
//			//calculate the thickness at this (x,y) position by adding the profile value to the average thickness (zExtent)
//			//(extra_thickness might also be negative)
//			thickness = zExtent+extra_thickness;
//// 			for (int k=zCornerCell; k<=zCornerCell+thickness-1; k++)
//// 			cout << zExtent << "," << extra_thickness << endl;
//// 			cout << xCornerCell << "," << yCornerCell << "," << zCornerCell << endl;
//// 			exit(-1);
//// cout << zCornerCell+thickness << endl;
//			for (k=zCornerCell+thickness; k<=kupper+1; k++)  //from bottom of the engraving profile to the top of the grid
//			{
//// if ((j>=yCornerCell)&&(j<=yCornerCell+1)&&(i>=xCornerCell)&&(i<=xCornerCell+1)&&(k==(zCornerCell+thickness))){
//// if ((j>=yCornerCell)&&(j<=yCornerCell+20)&&(i>=xCornerCell)&&(i<=xCornerCell)&&(k==(zCornerCell+thickness))){
//// cout << i << "," << j << "," << k << endl;
//// }
//				//Media_Ex
//				if ((k>=klower)&&(k<=kupper+1))
//				{
//					if ((j>=jleft)&&(j<=jright+1))
//					{
//						if ((i>=iback)&&(i<=ifront))
//						{
//								Media_Ex(i,j,k) = MaterialIdentifier.ElectricIndex_X; //engrave with the specified material
//						}
//					}
//				}
//
//				//Media_Ey
//				if ((k>=klower)&&(k<=kupper+1))
//				{
//					if ((j>=jleft)&&(j<=jright))
//					{
//						if ((i>=iback)&&(i<=ifront+1))
//						{
//								Media_Ey(i,j,k) = MaterialIdentifier.ElectricIndex_Y; //engrave with the specified material
//						}
//					}
//				}
//
//				//Media_Ez
//				if ((k>=klower)&&(k<=kupper))
//				{
//					if ((j>=jleft)&&(j<=jright+1))
//					{
//						if ((i>=iback)&&(i<=ifront+1))
//						{
//								Media_Ez(i,j,k) = MaterialIdentifier.ElectricIndex_Z; //engrave with the specified material
//						}
//					}
//				}
//			}
//		}
//	}
//
//	//finally, close the file
//	SurfaceEngravingProfileFile.close();
//	if (rank==0)
//	{
//		cout << "Surface-engraving profile read." << endl << endl;
//	}
//}

/****************************************************************************************************************/
/****************************************************************************************************************/
/** following functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/
/****************************************************************************************************************/
/****************************************************************************************************************/

void analyze_layering()
//Analyzes the layering structure, and builds the layering information arrays
{
	/** NOTE: **/
	/** This routine can only analyze isotropic layers in its current form. **/
	/** When the classes that use this routine are generalized to handle anisotropic layers, this routine will also be modified. **/
	/** Furthermore, the routine can only handle permittivity variations. Permeability support will be added in the near future. **/

	//the layers are numerically ordered from 0 to 'number_of_layers-1' beginning from the bottom (k=1) of the grid
	//initialize the number of layers and the layering arrays
	number_of_layers=1;

	eps_x_values_in_layers.resize(1);
	eps_y_values_in_layers.resize(1);
	eps_z_values_in_layers.resize(1);
	mu_x_values_in_layers.resize(1);
	mu_y_values_in_layers.resize(1);
	mu_z_values_in_layers.resize(1);
	cond_e_x_values_in_layers.resize(1);
	cond_e_y_values_in_layers.resize(1);
	cond_e_z_values_in_layers.resize(1);
	cond_h_x_values_in_layers.resize(1);
	cond_h_y_values_in_layers.resize(1);
	cond_h_z_values_in_layers.resize(1);

	LayerLowerZIndices.resize(1);
	LayerThicknesses.resize(1);
	IsLayerGrounded.resize(1);
	//first layer is made of whatever is at the bottom of the grid
	eps_x_values_in_layers(0) = eps_z(eps_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	eps_y_values_in_layers(0) = eps_z(eps_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	eps_z_values_in_layers(0) = eps_z(eps_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	mu_x_values_in_layers(0) = mu_x(mu_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	mu_y_values_in_layers(0) = mu_x(mu_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	mu_z_values_in_layers(0) = mu_x(mu_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_e_x_values_in_layers(0) = cond_e_z(cond_e_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_e_y_values_in_layers(0) = cond_e_z(cond_e_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_e_z_values_in_layers(0) = cond_e_z(cond_e_z_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_h_x_values_in_layers(0) = cond_h_x(cond_h_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_h_y_values_in_layers(0) = cond_h_x(cond_h_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	cond_h_z_values_in_layers(0) = cond_h_x(cond_h_x_indices_on_z_axis(1)); /** isotropic layer assumed! **/
	//lowermost index of the first layer is naturally 1
	LayerLowerZIndices(0) = 1;
	//thickness of the first layer is initialized to 1, and updated inside the loop below
	LayerThicknesses(0) = 1;
	//first layer is not grounded
	IsLayerGrounded(0)=false;

	if (NCELLS_Z+2*NPML>=2)
	{//sort of a silly check, but...
	for (int k=2; k<=NCELLS_Z+2*NPML; k++)
	{
		//is the material at k different than the one at k-1?
		//only look at the constitutive parameters located at half-integer positions in z
		if (!MaterialIsSame(1,1,eps_z_indices_on_z_axis(k),
							mu_x_indices_on_z_axis(k),mu_y_indices_on_z_axis(k),1,
							1,1,cond_e_z_indices_on_z_axis(k),
							cond_h_x_indices_on_z_axis(k),cond_h_y_indices_on_z_axis(k),1,
							1,1,eps_z_indices_on_z_axis(k-1),
							mu_x_indices_on_z_axis(k-1),mu_y_indices_on_z_axis(k-1),1,
							1,1,cond_e_z_indices_on_z_axis(k-1),
							cond_h_x_indices_on_z_axis(k-1),cond_h_y_indices_on_z_axis(k-1),1))
		{//new material is encountered
			number_of_layers++;
			//increase the array sizes by 1
			eps_x_values_in_layers.resizeAndPreserve(number_of_layers);
			eps_y_values_in_layers.resizeAndPreserve(number_of_layers);
			eps_z_values_in_layers.resizeAndPreserve(number_of_layers);
			mu_x_values_in_layers.resizeAndPreserve(number_of_layers);
			mu_y_values_in_layers.resizeAndPreserve(number_of_layers);
			mu_z_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_e_x_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_e_y_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_e_z_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_h_x_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_h_y_values_in_layers.resizeAndPreserve(number_of_layers);
			cond_h_z_values_in_layers.resizeAndPreserve(number_of_layers);
//			LayerElectricMaterial_X.resizeAndPreserve(number_of_layers);
//			LayerElectricMaterial_Y.resizeAndPreserve(number_of_layers);
//			LayerElectricMaterial_Z.resizeAndPreserve(number_of_layers);
//			LayerMagneticMaterial_X.resizeAndPreserve(number_of_layers);
//			LayerMagneticMaterial_Y.resizeAndPreserve(number_of_layers);
//			LayerMagneticMaterial_Z.resizeAndPreserve(number_of_layers);
			LayerLowerZIndices.resizeAndPreserve(number_of_layers);
			LayerThicknesses.resizeAndPreserve(number_of_layers);
			IsLayerGrounded.resizeAndPreserve(number_of_layers);
			//fill the info for the new layer
			eps_x_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			eps_y_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			eps_z_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			mu_x_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			mu_y_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			mu_z_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_e_x_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_e_y_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_e_z_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_h_x_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_h_y_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
			cond_h_z_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
//			LayerElectricMaterial_X(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//			LayerElectricMaterial_Y(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//			LayerElectricMaterial_Z(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//			LayerMagneticMaterial_X(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
//			LayerMagneticMaterial_Y(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
//			LayerMagneticMaterial_Z(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
			LayerLowerZIndices(number_of_layers-1) = k;
			LayerThicknesses(number_of_layers-1) = 1; //initialize the thickness to 1
			//check if the new layer is grounded or not
			if (eps_x_indices_on_z_axis(k)==PEC)  /** isotropy assumed! **/
			{
				IsLayerGrounded(number_of_layers-1) = true;
			}
			else
			{
				IsLayerGrounded(number_of_layers-1) = false;
			}
		}
		else
		{//no new material at this z position, but there can be a PEC sheet between k and k-1
			if (eps_x_indices_on_z_axis(k)==PEC)  /** isotropy assumed! **/
			{//there is a PEC sheet between k and k-1
				if (eps_z_indices_on_z_axis(k-1)!=PEC)  /** isotropy assumed! **/
				{//if the previous z position was not PEC, start a new grounded layer
					number_of_layers++;
					//increase the array sizes by 1
					eps_x_values_in_layers.resizeAndPreserve(number_of_layers);
					eps_y_values_in_layers.resizeAndPreserve(number_of_layers);
					eps_z_values_in_layers.resizeAndPreserve(number_of_layers);
					mu_x_values_in_layers.resizeAndPreserve(number_of_layers);
					mu_y_values_in_layers.resizeAndPreserve(number_of_layers);
					mu_z_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_e_x_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_e_y_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_e_z_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_h_x_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_h_y_values_in_layers.resizeAndPreserve(number_of_layers);
					cond_h_z_values_in_layers.resizeAndPreserve(number_of_layers);
//					LayerElectricMaterial_X.resizeAndPreserve(number_of_layers);
//					LayerElectricMaterial_Y.resizeAndPreserve(number_of_layers);
//					LayerElectricMaterial_Z.resizeAndPreserve(number_of_layers);
//					LayerMagneticMaterial_X.resizeAndPreserve(number_of_layers);
//					LayerMagneticMaterial_Y.resizeAndPreserve(number_of_layers);
//					LayerMagneticMaterial_Z.resizeAndPreserve(number_of_layers);
					LayerLowerZIndices.resizeAndPreserve(number_of_layers);
					LayerThicknesses.resizeAndPreserve(number_of_layers);
					IsLayerGrounded.resizeAndPreserve(number_of_layers);
					//fill the info for the new layer
					eps_x_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					eps_y_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					eps_z_values_in_layers(number_of_layers-1) = eps_z(eps_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					mu_x_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					mu_y_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					mu_z_values_in_layers(number_of_layers-1) = mu_x(mu_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_e_x_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_e_y_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_e_z_values_in_layers(number_of_layers-1) = cond_e_z(cond_e_z_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_h_x_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_h_y_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
					cond_h_z_values_in_layers(number_of_layers-1) = cond_h_x(cond_h_x_indices_on_z_axis(k)); /** isotropic layer assumed! **/
//					LayerElectricMaterial_X(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//					LayerElectricMaterial_Y(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//					LayerElectricMaterial_Z(number_of_layers-1) = Layering_e_z(k); /** isotropy assumed! **/
//					LayerMagneticMaterial_X(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
//					LayerMagneticMaterial_Y(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
//					LayerMagneticMaterial_Z(number_of_layers-1) = Layering_h_x(k); /** isotropy assumed! **/
					LayerLowerZIndices(number_of_layers-1) = k;
					LayerThicknesses(number_of_layers-1) = 1; //initialize the thickness to 1
					IsLayerGrounded(number_of_layers-1) = true;
				}
				else
				{//the previous z position was also PEC, so never mind...
					LayerThicknesses(number_of_layers-1) += 1; //increase the thickness of the current layer by 1
				}
			}
			else
			{//no new material at this z position, AND there is no PEC sheet between k and k-1, so we are still in the same layer
				LayerThicknesses(number_of_layers-1) += 1; //increase the thickness of the current layer by 1
			}
		}
	}
	}
// cout << "Layering_e_z is " << Layering_e_z << endl;
// cout << "Layering_h_x is " << Layering_h_x << endl;
// cout << "number_of_layers is " << number_of_layers << endl;
//if (rank==0)  cout << "eps_x_indices_in_layers is " << eps_x_indices_in_layers << endl;
//if (rank==0) cout << "eps_y_indices_in_layers is " << eps_y_indices_in_layers << endl;
//if (rank==0)  cout << "eps_z_indices_in_layers is " << eps_z_indices_in_layers << endl;
//if (rank==0)  cout << "mu_x_indices_in_layers is " << mu_x_indices_in_layers << endl;
//if (rank==0)  cout << "mu_y_indices_in_layers is " << mu_y_indices_in_layers << endl;
//if (rank==0)  cout << "mu_z_indices_in_layers is " << mu_z_indices_in_layers << endl;
//if (rank==0)  cout << "cond_e_x_indices_in_layers is " << cond_e_x_indices_in_layers << endl;
//if (rank==0)  cout << "cond_e_y_indices_in_layers is " << cond_e_y_indices_in_layers << endl;
//if (rank==0)  cout << "cond_e_z_indices_in_layers is " << cond_e_z_indices_in_layers << endl;
//if (rank==0)  cout << "cond_h_x_indices_in_layers is " << cond_h_x_indices_in_layers << endl;
//if (rank==0)  cout << "cond_h_y_indices_in_layers is " << cond_h_y_indices_in_layers << endl;
//if (rank==0)  cout << "cond_h_z_indices_in_layers is " << cond_h_z_indices_in_layers << endl;

// cout << "LayerLowerZIndices is " << LayerLowerZIndices << endl;
// cout << "LayerThicknesses is " << LayerThicknesses << endl;
// cout << "IsLayerGrounded is " << IsLayerGrounded << endl;
}

void find_extremal_constitutive_params()
{
	//if the slab extends to the upper limit of the grid, then update the relative permittivity and velocity in the uppermost layer
	/** NOTE: **/
	/** epsilon_r_upper, mu_r_upper, c_upper, epsilon_r_lower, mu_r_lower, c_lower are only defined for an isotropic medium!!**/
	/** Using the x component, as the y and z components are assumed to be the same **/
	/** This can be changed when the classes using these variables are generalized to handle anisotropic media **/
	epsilon_r_upper = eps_x_values_in_layers(number_of_layers-1); /** isotropy assumed! **/
	mu_r_upper = mu_x_values_in_layers(number_of_layers-1); /** isotropy assumed! **/
	cond_e_upper = cond_e_x_values_in_layers(number_of_layers-1); /** isotropy assumed! **/
	cond_h_upper = cond_h_x_values_in_layers(number_of_layers-1); /** isotropy assumed! **/
	c_upper = c/sqrt(epsilon_r_upper*mu_r_upper); /** assumed lossless **/

	//if the slab extends to the lower limit of the grid, then update the relative permittivity and velocity in the lowermost layer
	epsilon_r_lower = eps_x_values_in_layers(0); /** isotropy assumed! **/
	mu_r_lower = mu_x_values_in_layers(0); /** isotropy assumed! **/
	cond_e_lower = cond_e_x_values_in_layers(0); /** isotropy assumed! **/
	cond_h_lower = cond_h_x_values_in_layers(0); /** isotropy assumed! **/
	c_lower = c/sqrt(epsilon_r_lower*mu_r_lower); /** assumed lossless **/

	//Find the maximum and minimum permittivity/permeability values in the grid
	//max/min permittivities
	Array<float,1> eps_x_noPEC = eps_x(Range(1,eps_x.size()-1)); //exclude the PEC material
	Array<float,1> eps_y_noPEC = eps_y(Range(1,eps_y.size()-1)); //exclude the PEC material
	Array<float,1> eps_z_noPEC = eps_z(Range(1,eps_z.size()-1)); //exclude the PEC material
	epsilon_r_max_x = max(eps_x_noPEC); epsilon_r_min_x = min(eps_x_noPEC);
	epsilon_r_max_y = max(eps_y_noPEC); epsilon_r_min_y = min(eps_y_noPEC);
	epsilon_r_max_z = max(eps_z_noPEC); epsilon_r_min_z = min(eps_z_noPEC);
	epsilon_r_max = max(max(epsilon_r_max_x,epsilon_r_max_y),epsilon_r_max_z);
	epsilon_r_min = min(min(epsilon_r_min_x,epsilon_r_min_y),epsilon_r_min_z);
	//max/min permeabilities
	Array<float,1> mu_x_noPEC = mu_x(Range(1,mu_x.size()-1)); //exclude the PEC material
	Array<float,1> mu_y_noPEC = mu_y(Range(1,mu_y.size()-1)); //exclude the PEC material
	Array<float,1> mu_z_noPEC = mu_z(Range(1,mu_z.size()-1)); //exclude the PEC material
	mu_r_max_x = max(mu_x_noPEC); mu_r_min_x = min(mu_x_noPEC);
	mu_r_max_y = max(mu_y_noPEC); mu_r_min_y = min(mu_y_noPEC);
	mu_r_max_z = max(mu_z_noPEC); mu_r_min_z = min(mu_z_noPEC);
	mu_r_max = max(max(mu_r_max_x,mu_r_max_y),mu_r_max_z);
	mu_r_min = min(min(mu_r_min_x,mu_r_min_y),mu_r_min_z);
}
/****************************************************************************************************************/
/****************************************************************************************************************/
/** the above functions should eventually be removed (their job should be done by PlaceSlab, AddMaterial and similar functions) **/
/****************************************************************************************************************/
/****************************************************************************************************************/

inline bool MaterialIsSame(const eps_x_type& eps_x_index_1,
				   const eps_y_type& eps_y_index_1,
				   const eps_z_type& eps_z_index_1,
				   const mu_x_type& mu_x_index_1,
				   const mu_y_type& mu_y_index_1,
				   const mu_z_type& mu_z_index_1,
				   const cond_e_x_type& cond_e_x_index_1,
				   const cond_e_y_type& cond_e_y_index_1,
				   const cond_e_z_type& cond_e_z_index_1,
				   const cond_h_x_type& cond_h_x_index_1,
				   const cond_h_y_type& cond_h_y_index_1,
				   const cond_h_z_type& cond_h_z_index_1,
				   const eps_x_type& eps_x_index_2,
				   const eps_y_type& eps_y_index_2,
				   const eps_z_type& eps_z_index_2,
				   const mu_x_type& mu_x_index_2,
				   const mu_y_type& mu_y_index_2,
				   const mu_z_type& mu_z_index_2,
				   const cond_e_x_type& cond_e_x_index_2,
				   const cond_e_y_type& cond_e_y_index_2,
				   const cond_e_z_type& cond_e_z_index_2,
				   const cond_h_x_type& cond_h_x_index_2,
				   const cond_h_y_type& cond_h_y_index_2,
				   const cond_h_z_type& cond_h_z_index_2)
{
	//simplest version, but ignores the possibility of two material indices representing same material properties
// 	return Mat1==Mat2;
	//more stringent version, which requires the matching of the material properties to machine precision (in double format)
	const int relaxation_factor = 100;
	return ((abs(eps_x(eps_x_index_1)-eps_x(eps_x_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(eps_y(eps_y_index_1)-eps_y(eps_y_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(eps_z(eps_z_index_1)-eps_z(eps_z_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_x(mu_x_index_1)-mu_x(mu_x_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_y(mu_y_index_1)-mu_y(mu_y_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(mu_z(mu_z_index_1)-mu_z(mu_z_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_x(cond_e_x_index_1)-cond_e_x(cond_e_x_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_y(cond_e_y_index_1)-cond_e_y(cond_e_y_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_e_z(cond_e_z_index_1)-cond_e_z(cond_e_z_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_x(cond_h_x_index_1)-cond_h_x(cond_h_x_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_y(cond_h_y_index_1)-cond_h_y(cond_h_y_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor)
			&&(abs(cond_h_z(cond_h_z_index_1)-cond_h_z(cond_h_z_index_2))<LIBSTD_DBL_EPSILON*relaxation_factor));
}
/****************************************************************************************************************/
