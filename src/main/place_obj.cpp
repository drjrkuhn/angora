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

#include "place_obj.h"

#include "material/Cmat.h"
#include "shape/Cshape.h"

extern double dt;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

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
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;

extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<update_coeff_type,3> Pa_X,Pb_X,Pa_Y,Pb_Y,Pa_Z,Pb_Z;
extern Array<omega_p_x_type,3> omega_p_x_indices;
extern Array<omega_p_y_type,3> omega_p_y_indices;
extern Array<omega_p_z_type,3> omega_p_z_indices;
extern Array<tau_r_x_type,3> tau_r_x_indices;
extern Array<tau_r_y_type,3> tau_r_y_indices;
extern Array<tau_r_z_type,3> tau_r_z_indices;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,tau_r_x,tau_r_y,tau_r_z;


void place_obj(const Cmat& material, const_Cshape_shared_ptr shapeptr)
{
	int i,j,k;

//	double BlockBack,BlockFront,BlockLeft,BlockRight,BlockLower,BlockUpper;

//	double BlockBack = shapeptr->bounding_box.back_limit;
//	double BlockFront = shapeptr->bounding_box.front_limit;
//	double BlockLeft = shapeptr->bounding_box.left_limit;
//	double BlockRight = shapeptr->bounding_box.right_limit;
//	double BlockLower = shapeptr->bounding_box.lower_limit;
//	double BlockUpper = shapeptr->bounding_box.upper_limit;

	//get the bounding box limits for quicker placement
	int BlockBack = shapeptr->bounding_box_back_cell();
	int BlockFront = shapeptr->bounding_box_front_cell();
	int BlockLeft = shapeptr->bounding_box_left_cell();
	int BlockRight = shapeptr->bounding_box_right_cell();
	int BlockLower = shapeptr->bounding_box_lower_cell();
	int BlockUpper = shapeptr->bounding_box_upper_cell();

	/** TODO: When the resizeAndPreserve() overload with Range arguments is developed,
	do the resizing here, instead of in init.cpp or initgeom.cpp **/
//	//resize the polarization current arrays as necessary
//	if ((material.omega_p_exists())&&(material.omega_p_value()!=0))
//	{
//		//x components
//		if (BlockBack<J_p_x.lbound(firstDim))
//		{
//			J_p_x.resize(Range(BlockBack,J_p_x.ubound(firstDim)),
//									Range(J_p_x.lbound(secondDim),J_p_x.ubound(secondDim)),
//									Range(J_p_x.lbound(thirdDim),J_p_x.ubound(thirdDim)));
//		}
//	}

	//place the bulk of the object
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-1))
				{
					if (material.eps_x_exists())
					{
						eps_x_indices(i,j,k)=material.eps_x_index();
					}
					if (material.cond_e_x_exists())
					{
						cond_e_x_indices(i,j,k)=material.cond_e_x_index();
					}
					if (material.omega_p_x_exists())
					{
						omega_p_x_indices(i,j,k)=material.omega_p_x_index();
					}
					if (material.tau_r_x_exists())
					{
						tau_r_x_indices(i,j,k)=material.tau_r_x_index();
					}
					//check if the pole plasma frequency is nonzero
					if (omega_p_x(omega_p_x_indices(i,j,k))!=0)
					{//there is dispersion
						dispersion_exists_at_Ex_position(i,j,k) = true;
						Pa_X(i,j,k)=(1-dt/(2.0*tau_r_x(tau_r_x_indices(i,j,k))))
								   /(1+dt/(2.0*tau_r_x(tau_r_x_indices(i,j,k))));
						Pb_X(i,j,k)=0.5*dx*(1+Pa_X(i,j,k))*(pow2(omega_p_x(omega_p_x_indices(i,j,k)))*epsilon_0*dt/2.0)
								   /(1+dt/(2.0*tau_r_x(tau_r_x_indices(i,j,k))));
						Ca_X(i,j,k)=(1-dt*(cond_e_x(cond_e_x_indices(i,j,k))+2/(dx*(1+Pa_X(i,j,k)))*Pb_X(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))
								   /(1+dt*(cond_e_x(cond_e_x_indices(i,j,k))+2/(dx*(1+Pa_X(i,j,k)))*Pb_X(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
						Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx
								   /(1+dt*(cond_e_x(cond_e_x_indices(i,j,k))+2/(dx*(1+Pa_X(i,j,k)))*Pb_X(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
					}
					else
					{//there is NO dispersion
						dispersion_exists_at_Ex_position(i,j,k) = false;
						Pa_X(i,j,k)=1.0;
						Pb_X(i,j,k)=0.0;
						Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
						Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
					}
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-1))
				{
					if (material.eps_y_exists())
					{
						eps_y_indices(i,j,k)=material.eps_y_index();
					}
					if (material.cond_e_y_exists())
					{
						cond_e_y_indices(i,j,k)=material.cond_e_y_index();
					}
					if (material.omega_p_y_exists())
					{
						omega_p_y_indices(i,j,k)=material.omega_p_y_index();
					}
					if (material.tau_r_y_exists())
					{
						tau_r_y_indices(i,j,k)=material.tau_r_y_index();
					}
					//check if the pole plasma frequency is nonzero
					if (omega_p_y(omega_p_y_indices(i,j,k))!=0)
					{//there is dispersion
						dispersion_exists_at_Ey_position(i,j,k) = true;
						Pa_Y(i,j,k)=(1-dt/(2.0*tau_r_y(tau_r_y_indices(i,j,k))))
								   /(1+dt/(2.0*tau_r_y(tau_r_y_indices(i,j,k))));
						Pb_Y(i,j,k)=0.5*dx*(1+Pa_Y(i,j,k))*(pow2(omega_p_y(omega_p_y_indices(i,j,k)))*epsilon_0*dt/2.0)
								   /(1+dt/(2.0*tau_r_y(tau_r_y_indices(i,j,k))));
						Ca_Y(i,j,k)=(1-dt*(cond_e_y(cond_e_y_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Y(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))
								   /(1+dt*(cond_e_y(cond_e_y_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Y(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
						Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx
								   /(1+dt*(cond_e_y(cond_e_y_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Y(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
					}
					else
					{//there is NO dispersion
						dispersion_exists_at_Ey_position(i,j,k) = false;
						Pa_Y(i,j,k)=1.0;
						Pb_Y(i,j,k)=0.0;
						Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
						Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
					}
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-1,k-0.5))
				{
					if (material.eps_z_exists())
					{
						eps_z_indices(i,j,k)=material.eps_z_index();
					}
					if (material.cond_e_z_exists())
					{
						cond_e_z_indices(i,j,k)=material.cond_e_z_index();
					}
					if (material.omega_p_z_exists())
					{
						omega_p_z_indices(i,j,k)=material.omega_p_z_index();
					}
					if (material.tau_r_z_exists())
					{
						tau_r_z_indices(i,j,k)=material.tau_r_z_index();
					}
					//check if the pole plasma frequency is nonzero
					if (omega_p_z(omega_p_z_indices(i,j,k))!=0)
					{//there is dispersion
						dispersion_exists_at_Ez_position(i,j,k) = true;
						Pa_Z(i,j,k)=(1-dt/(2.0*tau_r_z(tau_r_z_indices(i,j,k))))
								   /(1+dt/(2.0*tau_r_z(tau_r_z_indices(i,j,k))));
						Pb_Z(i,j,k)=(0.5*dx*(1+Pa_Z(i,j,k))*pow2(omega_p_z(omega_p_z_indices(i,j,k)))*epsilon_0*dt/2.0)
								   /(1+dt/(2.0*tau_r_z(tau_r_z_indices(i,j,k))));
						Ca_Z(i,j,k)=(1-dt*(cond_e_z(cond_e_z_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Z(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))
								   /(1+dt*(cond_e_z(cond_e_z_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Z(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
						Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx
								   /(1+dt*(cond_e_z(cond_e_z_indices(i,j,k))+2/(dx*(1+Pa_Y(i,j,k)))*Pb_Z(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
					}
					else
					{//there is NO dispersion
						dispersion_exists_at_Ez_position(i,j,k) = false;
						Pa_Z(i,j,k)=1.0;
						Pb_Z(i,j,k)=0.0;
						Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
						Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
					}
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-0.5))
				{
					if (material.mu_x_exists())
					{
						mu_x_indices(i,j,k)=material.mu_x_index();
					}
					if (material.cond_h_x_exists())
					{
						cond_h_x_indices(i,j,k)=material.cond_h_x_index();
					}
					Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
					Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-0.5))
				{
					if (material.mu_y_exists())
					{
						mu_y_indices(i,j,k)=material.mu_y_index();
					}
					if (material.cond_h_y_exists())
					{
						cond_h_y_indices(i,j,k)=material.cond_h_y_index();
					}
					Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
					Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-0.5,k-1))
				{
					if (material.mu_z_exists())
					{
						mu_z_indices(i,j,k)=material.mu_z_index();
					}
					if (material.cond_h_z_exists())
					{
						cond_h_z_indices(i,j,k)=material.cond_h_z_index();
					}
					Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
					Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
				}
			}
		}
	}

	/***************************************/
	/** TODO: Interpolation at boundaries **/
	/***************************************/
}
