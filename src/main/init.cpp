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

//Defines several routines that initialize some of the global variables used throughout the simulation.

#include "headers.h"

#include "init.h"

#include "material/Cmat_types.h"

extern bool check_mode;

extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML,NSTEPS;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

//	Update coefficients for the E-field
extern Array<update_coeff_type,3> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
//	Update coefficients for the H-field
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;
//Constitutive parameter index arrays
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

extern bool dispersion_exists;
extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<double,3> J_p_x,J_p_y,J_p_z;
extern Array<update_coeff_type,3> Pa_X,Pb_X,Pa_Y,Pb_Y,Pa_Z,Pb_Z;
extern Array<omega_p_x_type,3> omega_p_x_indices;
extern Array<omega_p_y_type,3> omega_p_y_indices;
extern Array<omega_p_z_type,3> omega_p_z_indices;
extern Array<tau_r_x_type,3> tau_r_x_indices;
extern Array<tau_r_y_type,3> tau_r_y_indices;
extern Array<tau_r_z_type,3> tau_r_z_indices;

#ifndef MPI_DISABLE
extern MPI_Comm MPI_SubComm,MPI_CartSubComm;
#endif
extern int rank, nodes;
extern int rank_x, rank_y, rank_z;
extern int rank_behind,rank_front,rank_left,rank_right,rank_below,rank_above;
extern int nodes_x, nodes_y, nodes_z;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern Array<double,2> SendBuf_Hx_y,SendBuf_Hz_y,RecvBuf_Hx_y,RecvBuf_Hz_y,
 SendBuf_Hy_x,SendBuf_Hz_x,RecvBuf_Hy_x,RecvBuf_Hz_x,
 SendBuf_Hx_z,SendBuf_Hy_z,RecvBuf_Hx_z,RecvBuf_Hy_z;

extern int Ex_min_index_in_x,Ex_max_index_in_x,Ex_min_index_in_y,Ex_max_index_in_y,Ex_min_index_in_z,Ex_max_index_in_z;
extern int Ey_min_index_in_x,Ey_max_index_in_x,Ey_min_index_in_y,Ey_max_index_in_y,Ey_min_index_in_z,Ey_max_index_in_z;
extern int Ez_min_index_in_x,Ez_max_index_in_x,Ez_min_index_in_y,Ez_max_index_in_y,Ez_min_index_in_z,Ez_max_index_in_z;
extern int Hx_min_index_in_x,Hx_max_index_in_x,Hx_min_index_in_y,Hx_max_index_in_y,Hx_min_index_in_z,Hx_max_index_in_z;
extern int Hy_min_index_in_x,Hy_max_index_in_x,Hy_min_index_in_y,Hy_max_index_in_y,Hy_min_index_in_z,Hy_max_index_in_z;
extern int Hz_min_index_in_x,Hz_max_index_in_x,Hz_min_index_in_y,Hz_max_index_in_y,Hz_min_index_in_z,Hz_max_index_in_z;

extern int gridwidth_x, gridwidth_y, gridwidth_z;
extern void MPI_exit(const int& exitcode);

extern void init_geom();


void init_grid()
{
	/** First, initialize the sectioning info **/

	//First, determine the optimum nodes_x,nodes_y,nodes_z (number of sections in the x,y,z directions)
	//Scans through the possible factorizations of nodes such that the each node is assigned a section that is as close to a cube as possible (which is equivalent to minimum sum of side lengths, given the constant volume of each section).
	Array<TinyVector<int,3>,1> Factors;	//list of different factorizations [(nodes_x),(nodes_y),(nodes_z)] of (nodes).
	Factors.resize(1);
	Factors(0) = nodes,1,1;	//default factorization (inserted for safety, it will be re-inserted below)
	for (int i=1; i<=nodes; i++)
	{
		if (nodes==(nodes/i)*i)	//if nodes is divisible by i
		{
			for (int j=1; j<=nodes/i; j++)
			{
				if (nodes==(nodes/(i*j))*i*j)	//if nodes is divisible by (i*j)
				{
					int k=nodes/(i*j);
					Factors.resizeAndPreserve(Factors.size()+1);
					Factors(Factors.size()-1) = i,j,k;	//add the newly found factorization [i,j,k]
				}
			}
		}
	}
	double minsum = (NCELLS_X+NCELLS_Y+NCELLS_Z)+1.0;	//minimum sum of side lengths for the sections
	int minfactor;	//the position in Factors that holds the optimum factorization of nodes
	for (int n=0; n<Factors.size(); n++)
	{
		if (((double)NCELLS_X/Factors(n)(0)+(double)NCELLS_Y/Factors(n)(1)+(double)NCELLS_Z/Factors(n)(2))<=minsum)
		{
			minsum = (double)NCELLS_X/Factors(n)(0)+(double)NCELLS_Y/Factors(n)(1)+(double)NCELLS_Z/Factors(n)(2);
			minfactor = n;
		}
	}
	nodes_x = Factors(minfactor)(0);
	nodes_y = Factors(minfactor)(1);
	nodes_z = Factors(minfactor)(2);

///** REMOVE LATER **/
//	nodes_x = 5;
//	nodes_y = 1;
//	nodes_z = 1;
///** REMOVE LATER **/

	//Using nodes_x,nodes_y,nodes_z, create global communicator with cartesian topology
	int ndims = 3;
	int dims[3] = {nodes_x,nodes_y,nodes_z};
	int periods[3] = {0,0,0};	//non-periodic
	int reorder = 1;	//permit reorder
#ifndef MPI_DISABLE
	MPI_Cart_create(MPI_SubComm, ndims, dims, periods, reorder, &MPI_CartSubComm);
#endif

	//Store the ranks of the adjacent nodes (store MPI_PROC_NULL if the node is at the boundary of the global grid)
	int next = 1; //return the nodes in the immediate adjacency
#ifndef MPI_DISABLE
	MPI_Cart_shift(MPI_CartSubComm, 0, next, &rank_behind, &rank_front);
	MPI_Cart_shift(MPI_CartSubComm, 1, next, &rank_left, &rank_right);
	MPI_Cart_shift(MPI_CartSubComm, 2, next, &rank_below, &rank_above);
#else
	//these are not used if MPI is disabled, but they're defined anyway
	rank_behind = 0;
	rank_front = 0;
	rank_left = 0;
	rank_right = 0;
	rank_below = 0;
	rank_above = 0;
#endif


	//Now, determine node limits:
	// Ranking is z-major (then y, then x), which is default in MPI
	//	z
	//	|
	//	|---|---|---|
	//	| 1	| 3	| 5 |
	//	|---|---|---|
	//	| 0 | 2 | 4 |
	//	 --------------y
	//
	// For the above, nodes_x=1, nodes_y=3, nodes_z=2.
	// For the node with rank 5, rank_x=0, rank_y=2, rank_z=1.
	// If nodes_x>1, the node with rank 6 resides 1 cell in front of node 0.

	// *********************************
	// Determine rank_x, rank_y, rank_z:
	// *********************************
#ifndef MPI_DISABLE
	int cart_coords[3];
	MPI_Cart_coords(MPI_CartSubComm,rank,3,cart_coords);
	rank_x = cart_coords[0];
	rank_y = cart_coords[1];
	rank_z = cart_coords[2];
#else
	// if there is no MPI, then the ranks are all 0
	rank_x = 0;
	rank_y = 0;
	rank_z = 0;
#endif

	// In the below, the nodes that include the PML regions are made thinner than
	// those that include non-PML regions. This is because a PML cell requires approximately three times
	// the number of computations compared to a regular grid cell.

	double pml_versus_normal = 3.0;

	// *************************
	// Determine iback & ifront:
	// *************************

	//if there is only one node in the x direction:
	if (nodes_x==1)
	{
		iback=1;
		ifront=NCELLS_X+2*NPML;
	}
	else
	{
		//There is more than one node in the x direction:
		if (NCELLS_X<(nodes_x-2)*2*NPML)
		{
			if (rank==0)
			{
				cout << "Back and front PML regions must each be contained in one node." << endl
					<< "Too many nodes or PML too thick!" << endl;
			}
			exit(-1);
		}

		//Two nodes are reserved for the back and front PML blocks
		//Width of these nodes is NPML + a, where a is given below
		int a = (int)round((NCELLS_X-(pml_versus_normal*(nodes_x-2)+2)*NPML)/nodes_x);

		//Number of cells in the nodes that cover the rest of the main grid
		gridwidth_x = 0;
		if (nodes_x!=2)
		{
			gridwidth_x = (NCELLS_X-2*a)/(nodes_x-2);	//node width = (total width- total PML node width)/(number of nodes-2)
		}

		if (rank_x==0) //First node includes back PML block
		{
			iback=1;								//x-index of rearmost cell
			if (nodes_x==2)
			{//different treatment is needed if there are only 2 nodes in this direction
				ifront=NCELLS_X+NPML-a;	//x-index of foremost cell
			}
			else
			{
				ifront=NPML+a;	//x-index of foremost cell
			}
		}
		else if (rank_x==(nodes_x-1)) //Last node includes front PML block
		{
			iback=NCELLS_X+NPML-a+1;				//x-index of rearmost cell
			ifront=NCELLS_X+2*NPML;					//x-index of foremost cell
		}
		else if (rank_x==(nodes_x-2)) //Width might be irregular for node #(nodes-2)
		{
			iback = NPML+a+(rank_x-1)*gridwidth_x+1;	//x-index of rearmost cell
			ifront=NCELLS_X+NPML-a;						//x-index of foremost cell
		}
		else	//The rest of the nodes include parts of main grid
		{
			iback = NPML+a+(rank_x-1)*gridwidth_x+1;	//x-index of rearmost cell
			ifront = iback+gridwidth_x-1;				//x-index of foremost cell
		}
	}

	// *************************
	// Determine jleft & jright:
	// *************************

	//if there is only one node:
	if (nodes_y==1)
	{
		jleft=1;
		jright=NCELLS_Y+2*NPML;
	}
	else
	{
		//There is more than one node:
		if (NCELLS_Y<(nodes_y-2)*2*NPML)
		{
			if (rank==0)
			{
				cout << "Left and right PML regions must each be contained in one node." << endl
					<< "Too many nodes or PML too thick!" << endl;
			}
			exit(-1);
		}

		//Two nodes are reserved for the left and right PML blocks
		//Width of these nodes is NPML + a, where a is given below
		int a = (int)round((NCELLS_Y-(pml_versus_normal*(nodes_y-2)+2)*NPML)/nodes_y);

		//Number of cells in the nodes that cover the rest of the main grid
		gridwidth_y = 0;
		if (nodes_y!=2)
		{
			gridwidth_y = (NCELLS_Y-2*a)/(nodes_y-2);	//node width = (total width- total PML node width)/(number of nodes-2)
		}

		if (rank_y==0) //First node includes left PML block
		{
			jleft=1;								//y-index of leftmost cell
			if (nodes_y==2)
			{//different treatment is needed if there are only 2 nodes in this direction
				jright=NCELLS_Y+NPML-a;	//y-index of rightmost cell
			}
			else
			{
				jright=NPML+a;	//y-index of rightmost cell
			}
		}
		else if (rank_y==(nodes_y-1)) //Last node includes right PML block
		{
			jleft=NCELLS_Y+NPML-a+1;				//y-index of leftmost cell
			jright=NCELLS_Y+2*NPML;					//y-index of rightmost cell
		}
		else if (rank_y==(nodes_y-2)) //Width might be irregular for node #(nodes-2)
		{
			jleft = NPML+a+(rank_y-1)*gridwidth_y+1;	//y-index of leftmost cell
			jright=NCELLS_Y+NPML-a;						//y-index of rightmost cell
		}
		else	//The rest of the nodes include parts of main grid
		{
			jleft = NPML+a+(rank_y-1)*gridwidth_y+1;	//y-index of leftmost cell
			jright = jleft+gridwidth_y-1;				//y-index of rightmost cell
		}
	}

	// *************************
	// Determine klower & kupper:
	// *************************

	//if there is only one node:
	if (nodes_z==1)
	{
		klower=1;
		kupper=NCELLS_Z+2*NPML;
	}
	else
	{
		//There is more than one node:
		if (NCELLS_Z<(nodes_z-2)*2*NPML)
		{
			if (rank==0)
			{
				cout << "Lower and upper PML regions must each be contained in one node." << endl
					<< "Too many nodes or PML too thick!" << endl;
			}
			exit(-1);
		}

		//Two nodes are reserved for the lower and upper PML blocks
		//Width of these nodes is NPML + a, where a is given below
		int a = (int)round((NCELLS_Z-(pml_versus_normal*(nodes_z-2)+2)*NPML)/nodes_z);

		//Number of cells in the nodes that cover the rest of the main grid
		gridwidth_z = 0;
		if (nodes_z!=2)
		{
			gridwidth_z = (NCELLS_Z-2*a)/(nodes_z-2);	//node width = (total width- total PML node width)/(number of nodes-2)
		}

		if (rank_z==0) //First node includes lower PML block
		{
			klower=1;								//z-index of lowermost cell
			if (nodes_z==2)
			{//different treatment is needed if there are only 2 nodes in this direction
				kupper=NCELLS_Z+NPML-a;	//z-index of uppermost cell
			}
			else
			{
				kupper=NPML+a;	//z-index of uppermost cell
			}
		}
		else if (rank_z==(nodes_z-1)) //Last node includes upper PML block
		{
			klower=NCELLS_Z+NPML-a+1;				//z-index of lowermost cell
			kupper=NCELLS_Z+2*NPML;					//z-index of uppermost cell
		}
		else if (rank_z==(nodes_z-2)) //Width might be irregular for node #(nodes-2)
		{
			klower= NPML+a+(rank_z-1)*gridwidth_z+1;	//z-index of lowermost cell
			kupper=NCELLS_Z+NPML-a;						//z-index of uppermost cell
		}
		else	//The rest of the nodes include parts of main grid
		{
			klower=NPML+a+(rank_z-1)*gridwidth_z+1;		//z-index of lowermost cell
			kupper=klower+gridwidth_z-1;				//z-index of uppermost cell
		}
	}
	/** Sectioning info initialized**/

///** REMOVE LATER **/
//if (rank==0)
//{
//	iback=1;
//	ifront=3;
//}
//if (rank==1)
//{
//	iback=4;
//	ifront=20;
//}
//if (rank==2)
//{
//	iback=21;
//	ifront=47;
//}
//if (rank==3)
//{
//	iback=48;
//	ifront=57;
//}
//if (rank==4)
//{
//	iback=58;
//	ifront=60;
//}
//jleft=1;
//jright=NCELLS_Y+2*NPML;
//klower=1;
//kupper=NCELLS_Z+2*NPML;
///** REMOVE LATER **/

	/** If not in check mode, initialize the grid arrays **/
	if (!check_mode)
	{
		//Allocate field arrays and material pointer arrays using the sectioning info determined in init_sectioning()
		//Allocate field arrays
		Ex.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		Ey.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		Ez.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		Hx.resize(Range(iback,ifront+1),Range(jleft-1,jright+1),Range(klower-1,kupper+1));
		Hy.resize(Range(iback-1,ifront+1),Range(jleft,jright+1),Range(klower-1,kupper+1));
		Hz.resize(Range(iback-1,ifront+1),Range(jleft-1,jright+1),Range(klower,kupper+1));
		//Allocate update coefficient arrays
		Ca_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		Cb_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		Ca_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		Cb_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		Ca_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		Cb_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		Da_X.resize(Range(iback,ifront+1),Range(jleft-1,jright+1),Range(klower-1,kupper+1));
		Db_X.resize(Range(iback,ifront+1),Range(jleft-1,jright+1),Range(klower-1,kupper+1));
		Da_Y.resize(Range(iback-1,ifront+1),Range(jleft,jright+1),Range(klower-1,kupper+1));
		Db_Y.resize(Range(iback-1,ifront+1),Range(jleft,jright+1),Range(klower-1,kupper+1));
		Da_Z.resize(Range(iback-1,ifront+1),Range(jleft-1,jright+1),Range(klower,kupper+1));
		Db_Z.resize(Range(iback-1,ifront+1),Range(jleft-1,jright+1),Range(klower,kupper+1));
		//Allocate constitutive parameter index arrays
		eps_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		eps_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		eps_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		mu_x_indices.resize(Range(iback,ifront+1),Range(jleft-1,jright+1),Range(klower-1,kupper+1));
		mu_y_indices.resize(Range(iback-1,ifront+1),Range(jleft,jright+1),Range(klower-1,kupper+1));
		mu_z_indices.resize(Range(iback-1,ifront+1),Range(jleft-1,jright+1),Range(klower,kupper+1));
		cond_e_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		cond_e_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		cond_e_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		cond_h_x_indices.resize(Range(iback,ifront+1),Range(jleft-1,jright+1),Range(klower-1,kupper+1));
		cond_h_y_indices.resize(Range(iback-1,ifront+1),Range(jleft,jright+1),Range(klower-1,kupper+1));
		cond_h_z_indices.resize(Range(iback-1,ifront+1),Range(jleft-1,jright+1),Range(klower,kupper+1));

		//initialize MPI send and receive buffer arrays:
		//x direction (send-receive Hy and Hz)
		SendBuf_Hy_x.resize(Range(jleft,jright+1),Range(klower,kupper));
		SendBuf_Hz_x.resize(Range(jleft,jright),Range(klower,kupper+1));
		RecvBuf_Hy_x.resize(Range(jleft,jright+1),Range(klower,kupper));
		RecvBuf_Hz_x.resize(Range(jleft,jright),Range(klower,kupper+1));
		//y direction (send-receive Hx and Hz)
		SendBuf_Hx_y.resize(Range(iback,ifront+1),Range(klower,kupper));
		SendBuf_Hz_y.resize(Range(iback,ifront),Range(klower,kupper+1));
		RecvBuf_Hx_y.resize(Range(iback,ifront+1),Range(klower,kupper));
		RecvBuf_Hz_y.resize(Range(iback,ifront),Range(klower,kupper+1));
		//z direction (send-receive Hx and Hy)
		SendBuf_Hx_z.resize(Range(iback,ifront+1),Range(jleft,jright));
		SendBuf_Hy_z.resize(Range(iback,ifront),Range(jleft,jright+1));
		RecvBuf_Hx_z.resize(Range(iback,ifront+1),Range(jleft,jright));
		RecvBuf_Hy_z.resize(Range(iback,ifront),Range(jleft,jright+1));

		//initialize the dispersion indicator arrays
		dispersion_exists_at_Ex_position.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		dispersion_exists_at_Ey_position.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		dispersion_exists_at_Ez_position.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));

		/** TODO: When the resizeAndPreserve() overload with Range arguments is developed,
		do this resizing when a new object is placed in the grid.**/
		J_p_x.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		J_p_y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		J_p_z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		Pa_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		Pb_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		Pa_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		Pb_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		Pa_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		Pb_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		omega_p_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		omega_p_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		omega_p_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		tau_r_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
		tau_r_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
		tau_r_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));
		/************/
	}
	/** Grid arrays initialized **/

	//min. and max. indices of field components updated at the node
	Ex_min_index_in_x = iback;
	Ex_max_index_in_x = ifront;
	Ex_min_index_in_y = max(2,jleft);
	Ex_max_index_in_y = min(NCELLS_Y+2*NPML,jright+1);
	Ex_min_index_in_z = max(2,klower);
	Ex_max_index_in_z = min(NCELLS_Z+2*NPML,kupper+1);

	Ey_min_index_in_x = max(2,iback);
	Ey_max_index_in_x = min(NCELLS_X+2*NPML,ifront+1);
	Ey_min_index_in_y = jleft;
	Ey_max_index_in_y = jright;
	Ey_min_index_in_z = max(2,klower);
	Ey_max_index_in_z = min(NCELLS_Z+2*NPML,kupper+1);

	Ez_min_index_in_x = max(2,iback);
	Ez_max_index_in_x = min(NCELLS_X+2*NPML,ifront+1);
	Ez_min_index_in_y = max(2,jleft);
	Ez_max_index_in_y = min(NCELLS_Y+2*NPML,jright+1);
	Ez_min_index_in_z = klower;
	Ez_max_index_in_z = kupper;

	Hx_min_index_in_x = max(2,iback);
	Hx_max_index_in_x = min(NCELLS_X+2*NPML,ifront+1);
	Hx_min_index_in_y = max(1,jleft);
	Hx_max_index_in_y = min(NCELLS_Y+2*NPML,jright);
	Hx_min_index_in_z = max(1,klower);
	Hx_max_index_in_z = min(NCELLS_Z+2*NPML,kupper);

	Hy_min_index_in_x = max(1,iback);
	Hy_max_index_in_x = min(NCELLS_X+2*NPML,ifront);
	Hy_min_index_in_y = max(2,jleft);
	Hy_max_index_in_y = min(NCELLS_Y+2*NPML,jright+1);
	Hy_min_index_in_z = max(1,klower);
	Hy_max_index_in_z = min(NCELLS_Z+2*NPML,kupper);

	Hz_min_index_in_x = max(1,iback);
	Hz_max_index_in_x = min(NCELLS_X+2*NPML,ifront);
	Hz_min_index_in_y = max(1,jleft);
	Hz_max_index_in_y = min(NCELLS_Y+2*NPML,jright);
	Hz_min_index_in_z = max(2,klower);
	Hz_max_index_in_z = min(NCELLS_Z+2*NPML,kupper+1);

	//initialize the grid geometry
	init_geom();
}
