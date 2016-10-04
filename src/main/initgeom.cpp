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

//Defines functions that initialize the geometry in the grid.

#include "headers.h"

#include "initgeom.h"

#include "time_axis.h"

#include "material/Cmat_types.h"

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

//extern int num_of_distinct_eps_x,num_of_distinct_eps_y,num_of_distinct_eps_z;
//extern int num_of_distinct_mu_x,num_of_distinct_mu_y,num_of_distinct_mu_z;
//extern int num_of_distinct_cond_e_x,num_of_distinct_cond_e_y,num_of_distinct_cond_e_z;
//extern int num_of_distinct_cond_h_x,num_of_distinct_cond_h_y,num_of_distinct_cond_h_z;

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
extern Array<omega_p_x_type,3> omega_p_x_indices;
extern Array<omega_p_y_type,3> omega_p_y_indices;
extern Array<omega_p_z_type,3> omega_p_z_indices;
extern Array<tau_r_x_type,3> tau_r_x_indices;
extern Array<tau_r_y_type,3> tau_r_y_indices;
extern Array<tau_r_z_type,3> tau_r_z_indices;

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

//extern Array<double,1> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
//extern Array<double,1> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,tau_r_x,tau_r_y,tau_r_z;

extern Array<double,1> inv_kappa_e_x,inv_kappa_e_y,inv_kappa_e_z,inv_kappa_h_x,inv_kappa_h_y,inv_kappa_h_z;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;
extern double c_upper,c_lower;

extern bool dispersion_exists;
extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<double,3> J_p_x,J_p_y,J_p_z;
extern Array<update_coeff_type,3> Pa_X,Pb_X,Pa_Y,Pb_Y,Pa_Z,Pb_Z;

extern double max_field_value;
extern bool max_field_value_set_in_configfile;
extern double dB_accuracy;
extern bool dB_accuracy_set_in_configfile;


void init_geom()
{//resets the material information and assigns vacuum to the entire grid
	//the following is not automatic, since named materials are handled in the object Cmats
//	/** FIXME: This is not the most ideal solution **/
//	//initialize the array of named (tagged) materials to zero length
//	InitializeNamedMaterialArray();
//	/** FIXME: This is not the most ideal solution **/

	//set the initial time value to 0
//	allow_setting_initial_time_value();
	set_initial_time_value(0);

	//Initialize some global variables
	//initially, there are only 2 materials: vacuum and PEC
//	NumOfElectricMaterials_X = 2;
//	NumOfElectricMaterials_Y = 2;
//	NumOfElectricMaterials_Z = 2;
//	NumOfMagneticMaterials_X = 2;
//	NumOfMagneticMaterials_Y = 2;
//	NumOfMagneticMaterials_Z = 2;
//	num_of_distinct_eps_x = 2;
//	num_of_distinct_eps_y = 2;
//	num_of_distinct_eps_z = 2;
//	num_of_distinct_mu_x = 2;
//	num_of_distinct_mu_y = 2;
//	num_of_distinct_mu_z = 2;
//	num_of_distinct_cond_e_x = 2;
//	num_of_distinct_cond_e_y = 2;
//	num_of_distinct_cond_e_z = 2;
//	num_of_distinct_cond_h_x = 2;
//	num_of_distinct_cond_h_y = 2;
//	num_of_distinct_cond_h_z = 2;

	epsilon_r_upper=1; 	//relative permittivity of the uppermost layer (changed only through PlaceSlab in geometry.cpp)
	epsilon_r_lower=1; 	//relative permittivity of the lowermost layer (changed only through PlaceSlab in geometry.cpp)
	c_upper=299792458.;	//velocity of propagation in the uppermost layer (changed only through PlaceSlab in geometry.cpp)
	c_lower=299792458.;	//velocity of propagation in the lowermost layer (changed only through PlaceSlab in geometry.cpp)

	//set the maximum electric field value encountered in grid to 1, if not set in the config file already
	if (!max_field_value_set_in_configfile) max_field_value = 1;

	if (!dB_accuracy_set_in_configfile) dB_accuracy = -60;		//useful accuracy range of the grid (in -dB), default to -60 dB

	//Initialize field arrays
	Ex=0;
	Ey=0;
	Ez=0;
	Hx=0;
	Hy=0;
	Hz=0;
	//Initialize polarization current arrays
	J_p_x=0;
	J_p_y=0;
	J_p_z=0;
	//Place materials
	eps_x_indices=vacuum;
	eps_y_indices=vacuum;
	eps_z_indices=vacuum;
	mu_x_indices=vacuum;
	mu_y_indices=vacuum;
	mu_z_indices=vacuum;
	cond_e_x_indices=vacuum;
	cond_e_y_indices=vacuum;
	cond_e_z_indices=vacuum;
	cond_h_x_indices=vacuum;
	cond_h_y_indices=vacuum;
	cond_h_z_indices=vacuum;
	omega_p_x_indices=vacuum;
	omega_p_y_indices=vacuum;
	omega_p_z_indices=vacuum;
	tau_r_x_indices=vacuum;
	tau_r_y_indices=vacuum;
	tau_r_z_indices=vacuum;

	//Allocate and initialize constitutive parameter arrays
	eps_x.resize(2);
	eps_y.resize(2);
	eps_z.resize(2);
	mu_x.resize(2);
	mu_y.resize(2);
	mu_z.resize(2);
	cond_e_x.resize(2);
	cond_e_y.resize(2);
	cond_e_z.resize(2);
	cond_h_x.resize(2);
	cond_h_y.resize(2);
	cond_h_z.resize(2);
	omega_p_x.resize(2);
	omega_p_y.resize(2);
	omega_p_z.resize(2);
	tau_r_x.resize(2);
	tau_r_y.resize(2);
	tau_r_z.resize(2);

	//Initialize constitutive parameter vectors:
	//relative permittivities
	eps_x = 1;
	eps_y = 1;
	eps_z = 1;
	//relative permeabilities
	mu_x = 1;
	mu_y = 1;
	mu_z = 1;
	//electric conductivities
	cond_e_x = 0;
	cond_e_y = 0;
	cond_e_z = 0;
	//magnetic conductivities
	cond_h_x = 0;
	cond_h_y = 0;
	cond_h_z = 0;
	//Drude pole frequencies
	omega_p_x = 0;
	omega_p_y = 0;
	omega_p_z = 0;
	//Drude pole relaxation times
	tau_r_x = 0;
	tau_r_y = 0;
	tau_r_z = 0;

	//electric permittivity for PEC material (does not matter unless permittivity profile is visualized)
	//This is supposed to represent "infinity". Better and more portable value for the "maximum double value" can be found later.
	eps_x(PEC) = LIBSTD_DBL_MAX;
	eps_y(PEC) = LIBSTD_DBL_MAX;
	eps_z(PEC) = LIBSTD_DBL_MAX;
	mu_x(PEC) = 1;
	mu_y(PEC) = 1;
	mu_z(PEC) = 1;
	//electric conductivity for PEC material (does not matter unless conductivity profile is visualized)
	//This is supposed to represent "infinity". Better and more portable value for the "maximum double value" can be found later.
	cond_e_x(PEC) = LIBSTD_DBL_MAX;
	cond_e_y(PEC) = LIBSTD_DBL_MAX;
	cond_e_z(PEC) = LIBSTD_DBL_MAX;
	cond_h_x(PEC) = 1;
	cond_h_y(PEC) = 1;
	cond_h_z(PEC) = 1;

	//Allocate and initialize PML kappa parameters
	inv_kappa_e_x.resize(Range(1,NCELLS_X+2*NPML+1));
	inv_kappa_e_y.resize(Range(1,NCELLS_Y+2*NPML+1));
	inv_kappa_e_z.resize(Range(1,NCELLS_Z+2*NPML+1));
	inv_kappa_h_x.resize(Range(1,NCELLS_X+2*NPML));
	inv_kappa_h_y.resize(Range(1,NCELLS_Y+2*NPML));
	inv_kappa_h_z.resize(Range(1,NCELLS_Z+2*NPML));
	inv_kappa_e_x=1.0;
	inv_kappa_h_x=1.0;
	inv_kappa_e_y=1.0;
	inv_kappa_h_y=1.0;
	inv_kappa_e_z=1.0;
	inv_kappa_h_z=1.0;

	//Initialize the layering arrays (layering defined in the z-direction)
	eps_x_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));
	eps_y_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));
	eps_z_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	mu_x_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	mu_y_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	mu_z_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));
	cond_e_x_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));
	cond_e_y_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));
	cond_e_z_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	cond_h_x_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	cond_h_y_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML));
	cond_h_z_indices_on_z_axis.resize(Range(1,NCELLS_Z+2*NPML+1));

	eps_x_indices_on_z_axis = vacuum;
	eps_y_indices_on_z_axis = vacuum;
	eps_z_indices_on_z_axis = vacuum;
	mu_x_indices_on_z_axis = vacuum;
	mu_y_indices_on_z_axis = vacuum;
	mu_z_indices_on_z_axis = vacuum;
	cond_e_x_indices_on_z_axis = vacuum;
	cond_e_y_indices_on_z_axis = vacuum;
	cond_e_z_indices_on_z_axis = vacuum;
	cond_h_x_indices_on_z_axis = vacuum;
	cond_h_y_indices_on_z_axis = vacuum;
	cond_h_z_indices_on_z_axis = vacuum;

	// Update coefficients for default material (vacuum)
	double eaf_x = dt*cond_e_x(vacuum)/(2.0*eps_x(vacuum)*epsilon_0);
	Ca_X=(1-eaf_x)/(1+eaf_x);
	Cb_X=dt/eps_x(vacuum)/epsilon_0/dx/(1+eaf_x);
	double eaf_y = dt*cond_e_y(vacuum)/(2.0*eps_y(vacuum)*epsilon_0);
	Ca_Y=(1-eaf_y)/(1+eaf_y);
	Cb_Y=dt/eps_y(vacuum)/epsilon_0/dx/(1+eaf_y);
	double eaf_z = dt*cond_e_z(vacuum)/(2.0*eps_z(vacuum)*epsilon_0);
	Ca_Z=(1-eaf_z)/(1+eaf_z);
	Cb_Z=dt/eps_z(vacuum)/epsilon_0/dx/(1+eaf_z);
	double haf_x = dt*cond_h_x(vacuum)/(2.0*mu_x(vacuum)*mu_0);
	Da_X=(1-haf_x)/(1+haf_x);
	Db_X=dt/mu_x(vacuum)/mu_0/dx/(1+haf_x);
	double haf_y = dt*cond_h_y(vacuum)/(2.0*mu_y(vacuum)*mu_0);
	Da_Y=(1-haf_y)/(1+haf_y);
	Db_Y=dt/mu_y(vacuum)/mu_0/dx/(1+haf_y);
	double haf_z = dt*cond_h_z(vacuum)/(2.0*mu_z(vacuum)*mu_0);
	Da_Z=(1-haf_z)/(1+haf_z);
	Db_Z=dt/mu_z(vacuum)/mu_0/dx/(1+haf_z);

	//initialize the dispersion variables
	dispersion_exists = false;
	dispersion_exists_at_Ex_position = false;
	dispersion_exists_at_Ey_position = false;
	dispersion_exists_at_Ez_position = false;

	// Polarization-current update coefficients for default material (vacuum)
	Pa_X=1.0;
	Pb_X=0.0;
	Pa_Y=1.0;
	Pb_Y=0.0;
	Pa_Z=1.0;
	Pb_Z=0.0;

////if ((i==31)&&(j==31)&&(k==31))
////{
//	cout << Ca_X(33,33,33) << endl;
//	cout << Cb_X(33,33,33) << endl;
//	cout << Ca_Y(33,33,33) << endl;
//	cout << Cb_Y(33,33,33) << endl;
//	cout << Ca_Z(33,33,33) << endl;
//	cout << Cb_Z(33,33,33) << endl;
//	cout << Da_X(33,33,33) << endl;
//	cout << Db_X(33,33,33) << endl;
//	cout << Da_Y(33,33,33) << endl;
//	cout << Db_Y(33,33,33) << endl;
//	cout << Da_Z(33,33,33) << endl;
//	cout << Db_Z(33,33,33) << endl;
//	exit(-1);
////}
} //init_geom
