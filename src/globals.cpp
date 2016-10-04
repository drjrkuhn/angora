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

#include "headers.h"

#include "material/Cmat_types.h"

//Use MPI if not disabled
#ifndef MPI_DISABLE
#include <mpi.h>
#endif

//global grid variables

/*********************/
/*	ANGORA VERSION **/
/*********************/
//Full version number: MAJOR.MINOR.REVISION
int angora_version_major;		//major version number of the Angora package
int angora_version_minor;		//minor version number of the Angora package
int angora_version_revision;		//revision number of the Angora package


/*******************************/
/*	FILE AND DIRECTORY NAMES ***/
/*******************************/
//current working directory
string currentworkdir;
//base path for all other paths  (defaults to working directory)
string angora_basepath;
//output path
string OutputDir;
//input path
string InputDir;
//log output path
string LogOutputDir;
//recorder output path
string RecorderOutputDir;
//movie-recorder output path
string MovieRecorderOutputDir;
//line-recorder output path
string LineRecorderOutputDir;
//field-value-recorder output path
string FieldValueRecorderOutputDir;
//time-domain near-field-to-far-field transform output path
string TimeDomainNFFFTOutputDir;
//phasor-domain near-field-to-far-field transform output path
string PhasorDomainNFFFTOutputDir;


/***************************************/
/*	DEFAULT FILE AND DIRECTORY NAMES ***/
/***************************************/
//default output path
extern const string default_OutputDir = "output/";
//default input path
extern const string default_InputDir = "input/";
//default log output path
extern const string default_LogOutputDir = "log/";
//default log file name
extern const string default_LogFileName = "angora.log";
//default recorder output path
extern const string default_RecorderOutputDir = "recorder/";
//default movie-recorder output path
extern const string default_MovieRecorderOutputDir = "";
//default line-recorder output path
extern const string default_LineRecorderOutputDir = "";
//default field-value-recorder output path
extern const string default_FieldValueRecorderOutputDir = "";
//default movie filename
extern const string default_movie_filename = "MovieFile";
//default movie file extension
extern const string default_movie_fileextension = "amv";
//default line filename
extern const string default_line_filename = "LineFile";
//default line file extension
extern const string default_line_fileextension = "aln";
//default field-value filename
extern const string default_fieldvalue_filename = "FieldValueFile";
//default field-value file extension
extern const string default_fieldvalue_fileextension = "hd5";
//default time-domain near-field-to-far-field transform output path
extern const string default_TimeDomainNFFFTOutputDir = "nffft/td";
//default phasor-domain near-field-to-far-field transform output path
extern const string default_PhasorDomainNFFFTOutputDir = "nffft/pd";
//default time-domain NFFFT output filename
extern const string default_TimeDomainNFFFT_filename = "FarField_td";
//default time-domain NFFFT output file extension
extern const string default_td_nffft_fileextension = "hd5";
//default phasor-domain NFFFT output filename
extern const string default_PhasorDomainNFFFT_filename = "FarField_pd";
//default phasor-domain NFFFT output file extension
extern const string default_pd_nffft_fileextension = "hd5";
//default optical image output path
extern const string default_OpticalImagingOutputDir = "imaging/";
//default optical image output filename
extern const string default_image_filename = "Image";
//default optical image output file extension
extern const string default_image_fileextension = "hd5";


/********************************/
/*	GRID DIMENSION VARIABLES	*/
/********************************/
double courant,dx,dt;
int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;
int NSTEPS;
//	Indices of the origin cell
int OriginX,OriginY,OriginZ;

/********************/
/*	FIELD ARRAYS	*/
/********************/
//	E-field
Array<double,3> Ex,Ey,Ez;
//	H-field
Array<double,3> Hx,Hy,Hz;


/****************************************/
/*	LAYERING AND MATERIAL INFORMATION	*/
/****************************************/
//Material type information
extern const int PEC=0;
extern const int vacuum=1;

//Update coefficient arrays
//	Update coefficients for the E-field
Array<update_coeff_type,3> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
//	Update coefficients for the H-field
Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;
//Constitutive parameter index arrays
Array<eps_x_type,3> eps_x_indices;
Array<eps_y_type,3> eps_y_indices;
Array<eps_z_type,3> eps_z_indices;
Array<mu_x_type,3> mu_x_indices;
Array<mu_y_type,3> mu_y_indices;
Array<mu_z_type,3> mu_z_indices;
Array<cond_e_x_type,3> cond_e_x_indices;
Array<cond_e_y_type,3> cond_e_y_indices;
Array<cond_e_z_type,3> cond_e_z_indices;
Array<cond_h_x_type,3> cond_h_x_indices;
Array<cond_h_y_type,3> cond_h_y_indices;
Array<cond_h_z_type,3> cond_h_z_indices;
//Constitutive parameter value arrays
Array<float,1> eps_x,eps_y,eps_z;
Array<float,1> mu_x,mu_y,mu_z;
Array<float,1> cond_e_x,cond_e_y,cond_e_z;
Array<float,1> cond_h_x,cond_h_y,cond_h_z;

//Layering information
Array<eps_x_type,1> eps_x_indices_on_z_axis;
Array<eps_y_type,1> eps_y_indices_on_z_axis;
Array<eps_z_type,1> eps_z_indices_on_z_axis;
Array<mu_x_type,1> mu_x_indices_on_z_axis;
Array<mu_y_type,1> mu_y_indices_on_z_axis;
Array<mu_z_type,1> mu_z_indices_on_z_axis;
Array<cond_e_x_type,1> cond_e_x_indices_on_z_axis;
Array<cond_e_y_type,1> cond_e_y_indices_on_z_axis;
Array<cond_e_z_type,1> cond_e_z_indices_on_z_axis;
Array<cond_h_x_type,1> cond_h_x_indices_on_z_axis;
Array<cond_h_y_type,1> cond_h_y_indices_on_z_axis;
Array<cond_h_z_type,1> cond_h_z_indices_on_z_axis;
// Layers are counted from the bottom up, i.e., the lowermost layer is layer #0, and the uppermost layer is layer #N-1 (where N is the number of layers)
//Number of different layers in the grid
int number_of_layers;
//Array of material indices for each layer
//electric properties:
Array<float,1> eps_x_values_in_layers;
Array<float,1> eps_y_values_in_layers;
Array<float,1> eps_z_values_in_layers;
Array<float,1> mu_x_values_in_layers;
Array<float,1> mu_y_values_in_layers;
Array<float,1> mu_z_values_in_layers;
Array<float,1> cond_e_x_values_in_layers;
Array<float,1> cond_e_y_values_in_layers;
Array<float,1> cond_e_z_values_in_layers;
Array<float,1> cond_h_x_values_in_layers;
Array<float,1> cond_h_y_values_in_layers;
Array<float,1> cond_h_z_values_in_layers;
//Array of cell indices that define the lower boundary of each layer (index of the lowest cell that is included in each layer)
Array<int,1> LayerLowerZIndices;
//thicknesses (in grid cells) of each layer
Array<int,1> LayerThicknesses;
//boolean array denoting whether there is a perfect-electric-conductor (PEC) sheet below each layer
Array<bool,1> IsLayerGrounded;

//Coordinate stretching (kappa) arrays
Array<double,1> inv_kappa_e_x,inv_kappa_e_y,inv_kappa_e_z,inv_kappa_h_x,inv_kappa_h_y,inv_kappa_h_z;

//maximum and minimum constitutive parameters in the grid
double epsilon_r_max_x,epsilon_r_min_x,epsilon_r_max_y,epsilon_r_min_y,epsilon_r_max_z,epsilon_r_min_z;
double mu_r_max_x,mu_r_min_x,mu_r_max_y,mu_r_min_y,mu_r_max_z,mu_r_min_z;
double epsilon_r_max,epsilon_r_min;
double mu_r_max,mu_r_min;

double epsilon_r_upper,mu_r_upper,cond_e_upper,cond_h_upper; 	//constitutive params of the uppermost layer (changed only in load_geometry())
double epsilon_r_lower,mu_r_lower,cond_e_lower,cond_h_lower; 	//constitutive params of the lowermost layer (changed only in load_geometry())
double c_upper;	//velocity of propagation in the uppermost layer (changed only through PlaceSlab in geometry.cpp)
double c_lower;	//velocity of propagation in the lowermost layer (changed only through PlaceSlab in geometry.cpp)

//Dispersion-related variables
bool dispersion_exists;
Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
//Drude polarization current arrays
Array<double,3> J_p_x,J_p_y,J_p_z;
//update coefficients for the polarization current updates
//Pa corresponds to k_p, Pb corresponds to beta_p (Taflove&Hagness, 3rd ed., pg. 366)
Array<update_coeff_type,3> Pa_X,Pb_X,Pa_Y,Pb_Y,Pa_Z,Pb_Z;
//Drude pole frequency and relaxation time index arrays
Array<omega_p_x_type,3> omega_p_x_indices;
Array<omega_p_y_type,3> omega_p_y_indices;
Array<omega_p_z_type,3> omega_p_z_indices;
Array<tau_r_x_type,3> tau_r_x_indices;
Array<tau_r_y_type,3> tau_r_y_indices;
Array<tau_r_z_type,3> tau_r_z_indices;
//distinct Drude pole frequency and relaxation time values
Array<float,1> omega_p_x,omega_p_y,omega_p_z,tau_r_x,tau_r_y,tau_r_z;

/******************************/
/*	MULTIPLE-GRID VARIABLES   */
/******************************/
extern const int default_number_of_runs = 1;	//Number of FDTD simulation runs is 1 by default
int number_of_runs;	//Number of FDTD simulation runs
int GridIndex;	//index of the current simulation run
Array<bool,1> grid_is_enabled;	//indices of the grids that are enabled (default: all grids from 0 to NumOfGrids-1)

/********************/
/*	MPI VARIABLES   */
/********************/
#ifndef MPI_DISABLE
//MPI Communicators
MPI_Comm MPI_SubComm;		//subcommunicator for the grid
MPI_Comm MPI_CartSubComm;	//cartesian subcommunicator for the grid
MPI_Status Status;
#endif
//Basic MPI variables
int rank, nodes;
//Ranks of adjacent nodes
int rank_behind,rank_front,rank_left,rank_right,rank_below,rank_above;
//Grid limits determined by node_sectioning() and node_limits()
int rank_x, rank_y, rank_z;		//section rank (or coordinate) in the x, y, z directions
int nodes_x, nodes_y, nodes_z;	//number of sections in the x, y, z directions
int iback,ifront;	//x-indices of rearmost and foremost cells in the grid
					// that is maintained in this node
int jleft,jright;	//y-indices of leftmost and rightmost cells in the grid
					// that is maintained in this node
int klower,kupper;	//z-indices of lowermost and uppermost cells in the grid
					// that is maintained in this node
//min. and max. indices of field components at full-integer grid positions that are updated at the node
//int FullIntPosMin_x,FullIntPosMax_x;
//int FullIntPosMin_y,FullIntPosMax_y;
//int FullIntPosMin_z,FullIntPosMax_z;
int Ex_min_index_in_x,Ex_max_index_in_x,Ex_min_index_in_y,Ex_max_index_in_y,Ex_min_index_in_z,Ex_max_index_in_z;
int Ey_min_index_in_x,Ey_max_index_in_x,Ey_min_index_in_y,Ey_max_index_in_y,Ey_min_index_in_z,Ey_max_index_in_z;
int Ez_min_index_in_x,Ez_max_index_in_x,Ez_min_index_in_y,Ez_max_index_in_y,Ez_min_index_in_z,Ez_max_index_in_z;
int Hx_min_index_in_x,Hx_max_index_in_x,Hx_min_index_in_y,Hx_max_index_in_y,Hx_min_index_in_z,Hx_max_index_in_z;
int Hy_min_index_in_x,Hy_max_index_in_x,Hy_min_index_in_y,Hy_max_index_in_y,Hy_min_index_in_z,Hy_max_index_in_z;
int Hz_min_index_in_x,Hz_max_index_in_x,Hz_min_index_in_y,Hz_max_index_in_y,Hz_min_index_in_z,Hz_max_index_in_z;

int gridwidth_x, gridwidth_y, gridwidth_z;		//width of the x, y, z sections
//MPI send and receive buffer arrays:
Array<double,2> SendBuf_Hx_y,SendBuf_Hz_y,RecvBuf_Hx_y,RecvBuf_Hz_y,
 SendBuf_Hy_x,SendBuf_Hz_x,RecvBuf_Hy_x,RecvBuf_Hz_x,
 SendBuf_Hx_z,SendBuf_Hy_z,RecvBuf_Hx_z,RecvBuf_Hy_z;


/********************************************************/
/*	MAXIMUM FIELD VALUE IN GRID, AND USEFUL ACCURACY	*/
/********************************************************/
double max_field_value;	//maximum electric field value encountered in grid
double dB_accuracy;		//useful accuracy range of the grid (in -dB)


/********************************/
/*	GLOBAL ITERATION INDICES	*/
/********************************/
int i,j,k;
