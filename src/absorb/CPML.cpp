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

//Defines the PML object "CPML"

#include "headers.h"

#include "material/Cmat_types.h"

#include "CPML.h"

extern double dx,dt;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z;

extern Array<double,3> Ex,Ey,Ez;
extern Array<double,3> Hx,Hy,Hz;

extern Array<update_coeff_type,3> Ca_X,Cb_X,Ca_Y,Cb_Y,Ca_Z,Cb_Z;
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

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

extern Array<double,1> inv_kappa_e_x,inv_kappa_e_y,inv_kappa_e_z,inv_kappa_h_x,inv_kappa_h_y,inv_kappa_h_z;
extern Array<float,1> eps_x,eps_y,eps_z;

extern const int PEC;

extern int rank;
extern int rank_x, rank_y, rank_z;
extern int nodes_x, nodes_y, nodes_z;
extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern int Ex_min_index_in_x,Ex_max_index_in_x,Ex_min_index_in_y,Ex_max_index_in_y,Ex_min_index_in_z,Ex_max_index_in_z;
extern int Ey_min_index_in_x,Ey_max_index_in_x,Ey_min_index_in_y,Ey_max_index_in_y,Ey_min_index_in_z,Ey_max_index_in_z;
extern int Ez_min_index_in_x,Ez_max_index_in_x,Ez_min_index_in_y,Ez_max_index_in_y,Ez_min_index_in_z,Ez_max_index_in_z;
extern int Hx_min_index_in_x,Hx_max_index_in_x,Hx_min_index_in_y,Hx_max_index_in_y,Hx_min_index_in_z,Hx_max_index_in_z;
extern int Hy_min_index_in_x,Hy_max_index_in_x,Hy_min_index_in_y,Hy_max_index_in_y,Hy_min_index_in_z,Hy_max_index_in_z;
extern int Hz_min_index_in_x,Hz_max_index_in_x,Hz_min_index_in_y,Hz_max_index_in_y,Hz_min_index_in_z,Hz_max_index_in_z;

extern int i,j,k;


Cpml::Cpml(const int& my_pml_thickness, const double& SizeOfScatterer)
		: pml_thickness(my_pml_thickness)
{
	if (pml_thickness>0)
	{
		//PML Parameters
		m=4;	//Order for polynomial grading
		epsilon_eff_x=0;	//Effective permittivity in the yz interface plane will be found by arithmetic weighting
		int usable_thickness = 0;	//usable thickness of the layered structure
		k = NCELLS_Z+2*pml_thickness;
		do
		{
			usable_thickness++;
			epsilon_eff_x += eps_z(eps_z_indices_on_z_axis(k));  /** isotropy assumed **/
			k--;
		} while (k>=1); //terminate at k=1
		epsilon_eff_x = epsilon_eff_x/usable_thickness;

		epsilon_eff_y=epsilon_eff_x; /** isotropy assumed **/	//Same as epsilon_eff_x, from lateral symmetry
		epsilon_eff_lowerz=eps_x(eps_x_indices_on_z_axis(pml_thickness+1)); /** isotropy assumed **/  //Effective permittivity in the lower xy interface plane
		epsilon_eff_upperz=eps_x(eps_x_indices_on_z_axis(NCELLS_Z+pml_thickness+1)); /** isotropy assumed **/ //Effective permittivity in the upper xy interface plane
		sigma_coeff=1;		//sigma_max/sigma_optimum (sigma_opt as in Taflove 7.61)

		//Maximum sigma values in (x,y,z)-normal PML permittivity tensor (max. in back plane)
		sigma_x_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_x,0.5)*dx)*sigma_coeff;
		sigma_y_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_y,0.5)*dx)*sigma_coeff;
		sigma_lowerz_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_lowerz,0.5)*dx)*sigma_coeff;
		sigma_upperz_max=0.8*(m+1)/(eta_0*pow(epsilon_eff_upperz,0.5)*dx)*sigma_coeff;

		//Maximum kappa values in (x,y,z)-normal PML permittivity tensor (max. in back plane)
		kappa_x_max=1;
		kappa_y_max=1;
		kappa_z_max=1;

		//Maximum alpha in (x,y,z)-normal PML permittivity tensor (max. in front plane)
		//alpha_opt = c*eps/w , where w is the size of the scatterer
		// (see Berenger, "Numerical reflection from FDTD-PMLs - a comparison of the split PML with the unsplit and CFS PMLs", IEEE Ant.Prop. vol 50, Mar. 2002)
		alpha_x_max=c*epsilon_0/(dx*SizeOfScatterer);
		alpha_y_max=c*epsilon_0/(dx*SizeOfScatterer);
		alpha_z_max=c*epsilon_0/(dx*SizeOfScatterer);
	// 	alpha_x_max=1e-7;
	// 	alpha_y_max=1e-7;
	// 	alpha_z_max=1e-7;

		//Minimum alpha in (x,y,z)-normal PML permittivity tensor (min. in back plane)
		alpha_x_min=0;
		alpha_y_min=0;
		alpha_z_min=0;

		//Modify the coordinate-stretching kappa parameters in the PML region
		//(note the inverse definition: multiplication is faster than division)
		//BackX and FrontX
		for (int i=1;i<=pml_thickness+1;i++)
		{
			//BackX
			inv_kappa_e_x(i)=1.0/(1+(kappa_x_max-1)*pow((pml_thickness-i+1.0)/pml_thickness,m));
			//FrontX
			inv_kappa_e_x(NCELLS_X+2*pml_thickness+2-i)=inv_kappa_e_x(i);
		}

		for (int i=1;i<=pml_thickness;i++)
		{
			//BackX
			inv_kappa_h_x(i)=1.0/(1+(kappa_x_max-1)*pow((pml_thickness-i+0.5)/pml_thickness,m));
			//FrontX
			inv_kappa_h_x(NCELLS_X+2*pml_thickness+1-i)=inv_kappa_h_x(i);
		}

		//LeftY and RightY
		for (int j=1;j<=pml_thickness+1;j++)
		{
			//LeftY
			inv_kappa_e_y(j)=1.0/(1+(kappa_y_max-1)*pow((pml_thickness-j+1.0)/pml_thickness,m));
			//RightY
			inv_kappa_e_y(NCELLS_Y+2*pml_thickness+2-j)=inv_kappa_e_y(j);
		}
		for (int j=1;j<=pml_thickness;j++)
		{
			//LeftY
			inv_kappa_h_y(j)=1.0/(1+(kappa_y_max-1)*pow((pml_thickness-j+0.5)/pml_thickness,m));
			//RightY
			inv_kappa_h_y(NCELLS_Y+2*pml_thickness+1-j)=inv_kappa_h_y(j);
		}

		//UpperZ and LowerZ
		for (int k=1;k<=pml_thickness+1;k++)
		{
			//LowerZ
			inv_kappa_e_z(k)=1.0/(1+(kappa_z_max-1)*pow((pml_thickness-k+1.0)/pml_thickness,m));
			//UpperZ
			inv_kappa_e_z(NCELLS_Z+2*pml_thickness+2-k)=inv_kappa_e_z(k);
		}
		for (int k=1;k<=pml_thickness;k++)
		{
			//LowerZ
			inv_kappa_h_z(k)=1.0/(1+(kappa_z_max-1)*pow((pml_thickness-k+0.5)/pml_thickness,m));
			//UpperZ
			inv_kappa_h_z(NCELLS_Z+2*pml_thickness+1-k)=inv_kappa_h_z(k);
		}

		//Allocate and initialize PML current update parameters
		//	for updating the E-fields perpendicular to stretching axis:
		//Evaluated at the edge centers
		Apml_e_backx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the back PML slab
		Apml_e_frontx.resize(Range(NCELLS_X+pml_thickness+1,NCELLS_X+2*pml_thickness+1));	//for Ey and Ez in the front PML slab
		Apml_e_lefty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the left PML slab
		Apml_e_righty.resize(Range(NCELLS_Y+pml_thickness+1,NCELLS_Y+2*pml_thickness+1));	//for Ex and Ez in the right PML slab
		Apml_e_lowerz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the lower PML slab
		Apml_e_upperz.resize(Range(NCELLS_Z+pml_thickness+1,NCELLS_Z+2*pml_thickness+1));	//for Ex and Ey in the upper PML slab
		Bpml_e_backx.resize(Range(1,pml_thickness+1));	//for Ey and Ez in the back PML slab
		Bpml_e_frontx.resize(Range(NCELLS_X+pml_thickness+1,NCELLS_X+2*pml_thickness+1));	//for Ey and Ez in the front PML slab
		Bpml_e_lefty.resize(Range(1,pml_thickness+1));	//for Ex and Ez in the left PML slab
		Bpml_e_righty.resize(Range(NCELLS_Y+pml_thickness+1,NCELLS_Y+2*pml_thickness+1));	//for Ex and Ez in the right PML slab
		Bpml_e_lowerz.resize(Range(1,pml_thickness+1));	//for Ex and Ey in the lower PML slab
		Bpml_e_upperz.resize(Range(NCELLS_Z+pml_thickness+1,NCELLS_Z+2*pml_thickness+1));	//for Ex and Ey in the upper PML slab
		//	for updating the H-fields perpendicular to stretching axis:
		//Evaluated at the face centers
		Apml_h_backx.resize(Range(1,pml_thickness));		//for Hy and Hz in the back PML slab
		Apml_h_frontx.resize(Range(NCELLS_X+pml_thickness+1,NCELLS_X+2*pml_thickness));	//for Hy and Hz in the front PML slab
		Apml_h_lefty.resize(Range(1,pml_thickness));		//for Hx and Hz in the left PML slab
		Apml_h_righty.resize(Range(NCELLS_Y+pml_thickness+1,NCELLS_Y+2*pml_thickness));	//for Hx and Hz in the right PML slab
		Apml_h_lowerz.resize(Range(1,pml_thickness));	//for Hx and Hy in the lower PML slab
		Apml_h_upperz.resize(Range(NCELLS_Z+pml_thickness+1,NCELLS_Z+2*pml_thickness));	//for Hx and Hy in the upper PML slab
		Bpml_h_backx.resize(Range(1,pml_thickness));		//for Hy and Hz in the back PML slab
		Bpml_h_frontx.resize(Range(NCELLS_X+pml_thickness+1,NCELLS_X+2*pml_thickness));	//for Hy and Hz in the front PML slab
		Bpml_h_lefty.resize(Range(1,pml_thickness));		//for Hx and Hz in the left PML slab
		Bpml_h_righty.resize(Range(NCELLS_Y+pml_thickness+1,NCELLS_Y+2*pml_thickness));	//for Hx and Hz in the right PML slab
		Bpml_h_lowerz.resize(Range(1,pml_thickness));	//for Hx and Hy in the lower PML slab
		Bpml_h_upperz.resize(Range(NCELLS_Z+pml_thickness+1,NCELLS_Z+2*pml_thickness));	//for Hx and Hy in the upper PML slab
		//Initialize PML current update parameters (1-D)
		//	for the E-field:
		//Back and Front PML slab
		for (int i=1;i<=pml_thickness+1;i++)
		{
			double sigma_x,kappa_x,alpha_x;
			double r = pow((i-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_x=r*sigma_x_max;
			kappa_x=1+(kappa_x_max-1)*r;
			alpha_x=alpha_x_max+(alpha_x_min-alpha_x_max)*r;
			//BackX
			Apml_e_backx(pml_thickness+2-i)=sigma_x/(sigma_x*kappa_x+kappa_x*kappa_x*alpha_x)*
				(exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0)-1.0);
			Bpml_e_backx(pml_thickness+2-i)=exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0);
			//FrontX
			Apml_e_frontx(i+NCELLS_X+pml_thickness)=Apml_e_backx(pml_thickness+2-i);
			Bpml_e_frontx(i+NCELLS_X+pml_thickness)=Bpml_e_backx(pml_thickness+2-i);
		}
		//Left and Right PML slab
		for (int j=1;j<=pml_thickness+1;j++)
		{
			double sigma_y,kappa_y,alpha_y;
			double r = pow((j-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_y=r*sigma_y_max;
			kappa_y=1+(kappa_y_max-1)*r;
			alpha_y=alpha_y_max+(alpha_y_min-alpha_y_max)*r;
			//Left
			Apml_e_lefty(pml_thickness+2-j)=sigma_y/(sigma_y*kappa_y+kappa_y*kappa_y*alpha_y)*
				(exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0)-1.0);
			Bpml_e_lefty(pml_thickness+2-j)=exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0);
			//Right
			Apml_e_righty(j+NCELLS_Y+pml_thickness)=Apml_e_lefty(pml_thickness+2-j);
			Bpml_e_righty(j+NCELLS_Y+pml_thickness)=Bpml_e_lefty(pml_thickness+2-j);
		}
		//Lower and Upper PML slab
		for (int k=1;k<=pml_thickness+1;k++)
		{
			double sigma_z,kappa_z,alpha_z;
			double r = pow((k-1.0)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_z=r*sigma_lowerz_max;
			kappa_z=1+(kappa_z_max-1)*r;
			alpha_z=alpha_z_max+(alpha_z_min-alpha_z_max)*r;
			//Lower
			Apml_e_lowerz(pml_thickness+2-k)=sigma_z/(sigma_z*kappa_z+kappa_z*kappa_z*alpha_z)*
				(exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0)-1.0);
			Bpml_e_lowerz(pml_thickness+2-k)=exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0);
			//Upper
			Apml_e_upperz(k+NCELLS_Z+pml_thickness)=Apml_e_lowerz(pml_thickness+2-k);
			Bpml_e_upperz(k+NCELLS_Z+pml_thickness)=Bpml_e_lowerz(pml_thickness+2-k);
		}

		//	for the H-field:
		//Back and Front PML slab
		for (int i=1;i<=pml_thickness;i++)
		{
			double sigma_x,kappa_x,alpha_x;
			double r = pow((i-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_x=r*sigma_x_max;
			kappa_x=1+(kappa_x_max-1)*r;
			alpha_x=alpha_x_max*(1-r);
			//BackX
			Apml_h_backx(pml_thickness+1-i)=sigma_x/(sigma_x*kappa_x+kappa_x*kappa_x*alpha_x)*
				(exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0)-1.0);
			Bpml_h_backx(pml_thickness+1-i)=exp(-(sigma_x/kappa_x+alpha_x)*dt/epsilon_0);
			//FrontX
			Apml_h_frontx(i+NCELLS_X+pml_thickness)=Apml_h_backx(pml_thickness+1-i);
			Bpml_h_frontx(i+NCELLS_X+pml_thickness)=Bpml_h_backx(pml_thickness+1-i);
		}
		//Left and Right PML slab
		for (int j=1;j<=pml_thickness;j++)
		{
			double sigma_y,kappa_y,alpha_y;
			double r = pow((j-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_y=r*sigma_y_max;
			kappa_y=1+(kappa_y_max-1)*r;
			alpha_y=alpha_y_max*(1-r);
			//Left
			Apml_h_lefty(pml_thickness+1-j)=sigma_y/(sigma_y*kappa_y+kappa_y*kappa_y*alpha_y)*
				(exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0)-1.0);
			Bpml_h_lefty(pml_thickness+1-j)=exp(-(sigma_y/kappa_y+alpha_y)*dt/epsilon_0);
			//Right
			Apml_h_righty(j+NCELLS_Y+pml_thickness)=Apml_h_lefty(pml_thickness+1-j);
			Bpml_h_righty(j+NCELLS_Y+pml_thickness)=Bpml_h_lefty(pml_thickness+1-j);
		}
		//Lower and Upper PML slab
		for (int k=1;k<=pml_thickness;k++)
		{
			double sigma_z,kappa_z,alpha_z;
			double r = pow((k-0.5)/pml_thickness,m);	//r=0 at the interface, 1 at the backing plane
			sigma_z=r*sigma_lowerz_max;
			kappa_z=1+(kappa_z_max-1)*r;
			alpha_z=alpha_z_max*(1-r);
			//Lower
			Apml_h_lowerz(pml_thickness+1-k)=sigma_z/(sigma_z*kappa_z+kappa_z*kappa_z*alpha_z)*
				(exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0)-1.0);
			Bpml_h_lowerz(pml_thickness+1-k)=exp(-(sigma_z/kappa_z+alpha_z)*dt/epsilon_0);
			//Upper
			Apml_h_upperz(k+NCELLS_Z+pml_thickness)=Apml_h_lowerz(pml_thickness+1-k);
			Bpml_h_upperz(k+NCELLS_Z+pml_thickness)=Bpml_h_lowerz(pml_thickness+1-k);
		}

		Ex_min_index_in_y_left = max(2,jleft);
		Ex_max_index_in_y_left = min(pml_thickness+1,jright+1);
		Ex_min_index_in_y_right = max(NCELLS_Y+pml_thickness+1,jleft);
		Ex_max_index_in_y_right = min(NCELLS_Y+2*pml_thickness,jright+1);
		Ex_min_index_in_z_lower = max(2,klower);
		Ex_max_index_in_z_lower = min(pml_thickness+1,kupper+1);
		Ex_min_index_in_z_upper = max(NCELLS_Z+pml_thickness+1,klower);
		Ex_max_index_in_z_upper = min(NCELLS_Z+2*pml_thickness,kupper+1);
		Ey_min_index_in_x_back = max(2,iback);
		Ey_max_index_in_x_back = min(pml_thickness+1,ifront+1);
		Ey_min_index_in_x_front = max(NCELLS_X+pml_thickness+1,iback);
		Ey_max_index_in_x_front = min(NCELLS_X+2*pml_thickness,ifront+1);
		Ey_min_index_in_z_lower = max(2,klower);
		Ey_max_index_in_z_lower = min(pml_thickness+1,kupper+1);
		Ey_min_index_in_z_upper = max(NCELLS_Z+pml_thickness+1,klower);
		Ey_max_index_in_z_upper = min(NCELLS_Z+2*pml_thickness,kupper+1);
		Ez_min_index_in_x_back = max(2,iback);
		Ez_max_index_in_x_back = min(pml_thickness+1,ifront+1);
		Ez_min_index_in_x_front = max(NCELLS_X+pml_thickness+1,iback);
		Ez_max_index_in_x_front = min(NCELLS_X+2*pml_thickness,ifront+1);
		Ez_min_index_in_y_left = max(2,jleft);
		Ez_max_index_in_y_left = min(pml_thickness+1,jright+1);
		Ez_min_index_in_y_right = max(NCELLS_Y+pml_thickness+1,jleft);
		Ez_max_index_in_y_right = min(NCELLS_Y+2*pml_thickness,jright+1);

		Hx_min_index_in_y_left = max(1,jleft);
		Hx_max_index_in_y_left = min(pml_thickness,jright);
		Hx_min_index_in_y_right = max(NCELLS_Y+pml_thickness+1,jleft);
		Hx_max_index_in_y_right = min(NCELLS_Y+2*pml_thickness,jright);
		Hx_min_index_in_z_lower = max(1,klower);
		Hx_max_index_in_z_lower = min(pml_thickness,kupper);
		Hx_min_index_in_z_upper = max(NCELLS_Z+pml_thickness+1,klower);
		Hx_max_index_in_z_upper = min(NCELLS_Z+2*pml_thickness,kupper);
		Hy_min_index_in_x_back = max(1,iback);
		Hy_max_index_in_x_back = min(pml_thickness,ifront);
		Hy_min_index_in_x_front = max(NCELLS_X+pml_thickness+1,iback);
		Hy_max_index_in_x_front = min(NCELLS_X+2*pml_thickness,ifront);
		Hy_min_index_in_z_lower = max(1,klower);
		Hy_max_index_in_z_lower = min(pml_thickness,kupper);
		Hy_min_index_in_z_upper = max(NCELLS_Z+pml_thickness+1,klower);
		Hy_max_index_in_z_upper = min(NCELLS_Z+2*pml_thickness,kupper);
		Hz_min_index_in_x_back = max(1,iback);
		Hz_max_index_in_x_back = min(pml_thickness,ifront);
		Hz_min_index_in_x_front = max(NCELLS_X+pml_thickness+1,iback);
		Hz_max_index_in_x_front = min(NCELLS_X+2*pml_thickness,ifront);
		Hz_min_index_in_y_left = max(1,jleft);
		Hz_max_index_in_y_left = min(pml_thickness,jright);
		Hz_min_index_in_y_right = max(NCELLS_Y+pml_thickness+1,jleft);
		Hz_max_index_in_y_right = min(NCELLS_Y+2*pml_thickness,jright);

		//Allocate and initialize PML current arrays
		if (Ey_max_index_in_x_back>=Ey_min_index_in_x_back)
		{
			dx_times_Psi_Eyx_BackX.resize(Range(Ey_min_index_in_x_back,Ey_max_index_in_x_back),Range(jleft,jright),Range(klower,kupper+1));
			dx_times_Psi_Ezx_BackX.resize(Range(Ez_min_index_in_x_back,Ez_max_index_in_x_back),Range(jleft,jright+1),Range(klower,kupper));
			dx_times_Psi_Eyx_BackX=0;
			dx_times_Psi_Ezx_BackX=0;
		}
		if (Hy_max_index_in_x_back>=Hy_min_index_in_x_back)
		{
			dx_times_Psi_Hyx_BackX.resize(Range(Hy_min_index_in_x_back,Hy_max_index_in_x_back),Range(jleft,jright+1),Range(klower,kupper));
			dx_times_Psi_Hzx_BackX.resize(Range(Hz_min_index_in_x_back,Hz_max_index_in_x_back),Range(jleft,jright),Range(klower,kupper+1));
			dx_times_Psi_Hyx_BackX=0;
			dx_times_Psi_Hzx_BackX=0;
		}

		if (Ey_max_index_in_x_front>=Ey_min_index_in_x_front)
		{
			dx_times_Psi_Eyx_FrontX.resize(Range(Ey_min_index_in_x_front,Ey_max_index_in_x_front),Range(jleft,jright),Range(klower,kupper+1));
			dx_times_Psi_Ezx_FrontX.resize(Range(Ez_min_index_in_x_front,Ez_max_index_in_x_front),Range(jleft,jright+1),Range(klower,kupper));
			dx_times_Psi_Eyx_FrontX=0;
			dx_times_Psi_Ezx_FrontX=0;
		}
		if (Hy_max_index_in_x_front>=Hy_min_index_in_x_front)
		{
			dx_times_Psi_Hyx_FrontX.resize(Range(Hy_min_index_in_x_front,Hy_max_index_in_x_front),Range(jleft,jright+1),Range(klower,kupper));
			dx_times_Psi_Hzx_FrontX.resize(Range(Hz_min_index_in_x_front,Hz_max_index_in_x_front),Range(jleft,jright),Range(klower,kupper+1));
			dx_times_Psi_Hyx_FrontX=0;
			dx_times_Psi_Hzx_FrontX=0;
		}

		if (Ex_max_index_in_y_left>=Ex_min_index_in_y_left)
		{
			dx_times_Psi_Exy_LeftY.resize(Range(iback,ifront),Range(Ex_min_index_in_y_left,Ex_max_index_in_y_left),Range(klower,kupper+1));
			dx_times_Psi_Ezy_LeftY.resize(Range(iback,ifront+1),Range(Ez_min_index_in_y_left,Ez_max_index_in_y_left),Range(klower,kupper));
			dx_times_Psi_Exy_LeftY=0;
			dx_times_Psi_Ezy_LeftY=0;
		}
		if (Hx_max_index_in_y_left>=Hx_min_index_in_y_left)
		{
			dx_times_Psi_Hxy_LeftY.resize(Range(iback,ifront+1),Range(Hx_min_index_in_y_left,Hx_max_index_in_y_left),Range(klower,kupper));
			dx_times_Psi_Hzy_LeftY.resize(Range(iback,ifront),Range(Hz_min_index_in_y_left,Hz_max_index_in_y_left),Range(klower,kupper+1));
			dx_times_Psi_Hxy_LeftY=0;
			dx_times_Psi_Hzy_LeftY=0;
		}

		if (Ex_max_index_in_y_right>=Ex_min_index_in_y_right)
		{
			dx_times_Psi_Exy_RightY.resize(Range(iback,ifront),Range(Ex_min_index_in_y_right,Ex_max_index_in_y_right),Range(klower,kupper+1));
			dx_times_Psi_Ezy_RightY.resize(Range(iback,ifront+1),Range(Ez_min_index_in_y_right,Ez_max_index_in_y_right),Range(klower,kupper));
			dx_times_Psi_Exy_RightY=0;
			dx_times_Psi_Ezy_RightY=0;
		}
		if (Hx_max_index_in_y_right>=Hx_min_index_in_y_right)
		{
			dx_times_Psi_Hxy_RightY.resize(Range(iback,ifront+1),Range(Hx_min_index_in_y_right,Hx_max_index_in_y_right),Range(klower,kupper));
			dx_times_Psi_Hzy_RightY.resize(Range(iback,ifront),Range(Hz_min_index_in_y_right,Hz_max_index_in_y_right),Range(klower,kupper+1));
			dx_times_Psi_Hxy_RightY=0;
			dx_times_Psi_Hzy_RightY=0;
		}

		if (Ex_max_index_in_z_lower>=Ex_min_index_in_z_lower)
		{
			dx_times_Psi_Exz_LowerZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(Ex_min_index_in_z_lower,Ex_max_index_in_z_lower));
			dx_times_Psi_Eyz_LowerZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(Ey_min_index_in_z_lower,Ey_max_index_in_z_lower));
			dx_times_Psi_Exz_LowerZ=0;
			dx_times_Psi_Eyz_LowerZ=0;
		}
		if (Hx_max_index_in_z_lower>=Hx_min_index_in_z_lower)
		{
			dx_times_Psi_Hxz_LowerZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(Hx_min_index_in_z_lower,Hx_max_index_in_z_lower));
			dx_times_Psi_Hyz_LowerZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(Hy_min_index_in_z_lower,Hy_max_index_in_z_lower));
			dx_times_Psi_Hxz_LowerZ=0;
			dx_times_Psi_Hyz_LowerZ=0;
		}

		if (Ex_max_index_in_z_upper>=Ex_min_index_in_z_upper)
		{
			dx_times_Psi_Exz_UpperZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(Ex_min_index_in_z_upper,Ex_max_index_in_z_upper));
			dx_times_Psi_Eyz_UpperZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(Ey_min_index_in_z_upper,Ey_max_index_in_z_upper));
			dx_times_Psi_Exz_UpperZ=0;
			dx_times_Psi_Eyz_UpperZ=0;
		}
		if (Hx_max_index_in_z_upper>=Hx_min_index_in_z_upper)
		{
			dx_times_Psi_Hxz_UpperZ.resize(Range(iback,ifront+1),Range(jleft,jright),Range(Hx_min_index_in_z_upper,Hx_max_index_in_z_upper));
			dx_times_Psi_Hyz_UpperZ.resize(Range(iback,ifront),Range(jleft,jright+1),Range(Hy_min_index_in_z_upper,Hy_max_index_in_z_upper));
			dx_times_Psi_Hxz_UpperZ=0;
			dx_times_Psi_Hyz_UpperZ=0;
		}
	}
}

void Cpml::UpdateE()
{
	if (pml_thickness>0)
	{
		UpdateE_X();
		UpdateE_Y();
		UpdateE_Z();
	}
}

void Cpml::UpdateH()
{
	if (pml_thickness>0)
	{
		UpdateH_X();
		UpdateH_Y();
		UpdateH_Z();
	}
}

void Cpml::UpdateE_X()
{
		//dx_times_Psi_Eyx_BackX
		for (i=Ey_min_index_in_x_back; i<=Ey_max_index_in_x_back; i++){
			for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
				for (k=Ey_min_index_in_z; k<=Ey_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Eyx_BackX(i,j,k);
					//BackX
					temp=Bpml_e_backx(i)*temp + Apml_e_backx(i)*(Hz(i,j,k)-Hz(i-1,j,k));
					Ey(i,j,k)+=Cb_Y(i,j,k)*(- temp);
				}
			}
		}
		//dx_times_Psi_Ezx_BackX
		for (i=Ez_min_index_in_x_back; i<=Ez_max_index_in_x_back; i++){
			for (j=Ez_min_index_in_y; j<=Ez_max_index_in_y; j++){
				for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Ezx_BackX(i,j,k);
					//BackX
					temp=Bpml_e_backx(i)*temp + Apml_e_backx(i)*(Hy(i,j,k)-Hy(i-1,j,k));
					Ez(i,j,k)+=Cb_Z(i,j,k)*temp;
				}
			}
		}

		//dx_times_Psi_Eyx_FrontX
		for (i=Ey_min_index_in_x_front; i<=Ey_max_index_in_x_front; i++){
			for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
				for (k=Ey_min_index_in_z; k<=Ey_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Eyx_FrontX(i,j,k);
					//FrontX
					temp=Bpml_e_frontx(i)*temp + Apml_e_frontx(i)*(Hz(i,j,k)-Hz(i-1,j,k));
					Ey(i,j,k)+=Cb_Y(i,j,k)*(- temp);
				}
			}
		}
		//dx_times_Psi_Ezx_FrontX
		for (i=Ez_min_index_in_x_front; i<=Ez_max_index_in_x_front; i++){
			for (j=Ez_min_index_in_y; j<=Ez_max_index_in_y; j++){
				for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Ezx_FrontX(i,j,k);
					//FrontX
					temp=Bpml_e_frontx(i)*temp + Apml_e_frontx(i)*(Hy(i,j,k)-Hy(i-1,j,k));
					Ez(i,j,k)+=Cb_Z(i,j,k)*temp;
				}
			}
		}
}

void Cpml::UpdateH_X()
{
		//dx_times_Psi_Hyx_BackX
		for (i=Hy_min_index_in_x_back; i<=Hy_max_index_in_x_back; i++){
			for (j=Hy_min_index_in_y; j<=Hy_max_index_in_y; j++){
				for (k=Hy_min_index_in_z; k<=Hy_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hyx_BackX(i,j,k);
					//BackX
					temp=Bpml_h_backx(i)*temp + Apml_h_backx(i)*(Ez(i+1,j,k)-Ez(i,j,k));
					Hy(i,j,k)+=Db_Y(i,j,k)*temp;
				}
			}
		}
		//dx_times_Psi_Hzx_BackX
		for (i=Hz_min_index_in_x_back; i<=Hz_max_index_in_x_back; i++){
			for (j=Hz_min_index_in_y; j<=Hz_max_index_in_y; j++){
				for (k=Hz_min_index_in_z; k<=Hz_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hzx_BackX(i,j,k);
					//BackX
					temp=Bpml_h_backx(i)*temp + Apml_h_backx(i)*(Ey(i+1,j,k)-Ey(i,j,k));
					Hz(i,j,k)+=Db_Z(i,j,k)*(- temp);
				}
			}
		}

		//dx_times_Psi_Hyx_FrontX
		for (i=Hy_min_index_in_x_front; i<=Hy_max_index_in_x_front; i++){
			for (j=Hy_min_index_in_y; j<=Hy_max_index_in_y; j++){
				for (k=Hy_min_index_in_z; k<=Hy_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hyx_FrontX(i,j,k);
					//FrontX
					temp=Bpml_h_frontx(i)*temp + Apml_h_frontx(i)*(Ez(i+1,j,k)-Ez(i,j,k));
					Hy(i,j,k)+=Db_Y(i,j,k)*temp;
				}
			}
		}
		//dx_times_Psi_Hzx_FrontX
		for (i=Hz_min_index_in_x_front; i<=Hz_max_index_in_x_front; i++){
			for (j=Hz_min_index_in_y; j<=Hz_max_index_in_y; j++){
				for (k=Hz_min_index_in_z; k<=Hz_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hzx_FrontX(i,j,k);
					//FrontX
					temp=Bpml_h_frontx(i)*temp + Apml_h_frontx(i)*(Ey(i+1,j,k)-Ey(i,j,k));
					Hz(i,j,k)+=Db_Z(i,j,k)*(- temp);
				}
			}
		}
}

void Cpml::UpdateE_Y()
{
		//dx_times_Psi_Exy_LeftY
		for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
			for (j=Ex_min_index_in_y_left; j<=Ex_max_index_in_y_left; j++){
				for (k=Ex_min_index_in_z; k<=Ex_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Exy_LeftY(i,j,k);
					//LeftY
					temp=Bpml_e_lefty(j)*temp + Apml_e_lefty(j)*(Hz(i,j,k)-Hz(i,j-1,k));
					Ex(i,j,k)+=Cb_X(i,j,k)*temp;
				}
			}
		}
		//dx_times_Psi_Ezy_LeftY
		for (i=Ez_min_index_in_x; i<=Ez_max_index_in_x; i++){
			for (j=Ez_min_index_in_y_left; j<=Ez_max_index_in_y_left; j++){
				for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Ezy_LeftY(i,j,k);
					//LeftY
					temp=Bpml_e_lefty(j)*temp + Apml_e_lefty(j)*(Hx(i,j,k)-Hx(i,j-1,k));
					Ez(i,j,k)+=Cb_Z(i,j,k)*(- temp);
				}
			}
		}

		//dx_times_Psi_Exy_RightY
		for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
			for (j=Ex_min_index_in_y_right; j<=Ex_max_index_in_y_right; j++){
				for (k=Ex_min_index_in_z; k<=Ex_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Exy_RightY(i,j,k);
					//RightY
					temp=Bpml_e_righty(j)*temp + Apml_e_righty(j)*(Hz(i,j,k)-Hz(i,j-1,k));
					Ex(i,j,k)+=Cb_X(i,j,k)*temp;
				}
			}
		}
		//dx_times_Psi_Ezy_RightY
		for (i=Ez_min_index_in_x; i<=Ez_max_index_in_x; i++){
			for (j=Ez_min_index_in_y_right; j<=Ez_max_index_in_y_right; j++){
				for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Ezy_RightY(i,j,k);
					//RightY
					temp=Bpml_e_righty(j)*temp + Apml_e_righty(j)*(Hx(i,j,k)-Hx(i,j-1,k));
					Ez(i,j,k)+=Cb_Z(i,j,k)*(- temp);
				}
			}
		}
}

void Cpml::UpdateH_Y()
{
		//dx_times_Psi_Hxy_LeftY
		for (i=Hx_min_index_in_x; i<=Hx_max_index_in_x; i++){
			for (j=Hx_min_index_in_y_left; j<=Hx_max_index_in_y_left; j++){
				for (k=Hx_min_index_in_z; k<=Hx_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hxy_LeftY(i,j,k);
					//LeftY
					temp=Bpml_h_lefty(j)*temp + Apml_h_lefty(j)*(Ez(i,j+1,k)-Ez(i,j,k));
					Hx(i,j,k)+=Db_X(i,j,k)*(- temp);
				}
			}
		}
		//dx_times_Psi_Hzy_LeftY
		for (i=Hz_min_index_in_x; i<=Hz_max_index_in_x; i++){
			for (j=Hz_min_index_in_y_left; j<=Hz_max_index_in_y_left; j++){
				for (k=Hz_min_index_in_z; k<=Hz_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hzy_LeftY(i,j,k);
					//LeftY
					temp=Bpml_h_lefty(j)*temp + Apml_h_lefty(j)*(Ex(i,j+1,k)-Ex(i,j,k));
					Hz(i,j,k)+=Db_Z(i,j,k)*temp;
				}
			}
		}

		//dx_times_Psi_Hxy_RightY
		for (i=Hx_min_index_in_x; i<=Hx_max_index_in_x; i++){
			for (j=Hx_min_index_in_y_right; j<=Hx_max_index_in_y_right; j++){
				for (k=Hx_min_index_in_z; k<=Hx_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hxy_RightY(i,j,k);
					//RightY
					temp=Bpml_h_righty(j)*temp + Apml_h_righty(j)*(Ez(i,j+1,k)-Ez(i,j,k));
					Hx(i,j,k)+=Db_X(i,j,k)*(- temp);
				}
			}
		}
		//dx_times_Psi_Hzy_RightY
		for (i=Hz_min_index_in_x; i<=Hz_max_index_in_x; i++){
			for (j=Hz_min_index_in_y_right; j<=Hz_max_index_in_y_right; j++){
				for (k=Hz_min_index_in_z; k<=Hz_max_index_in_z; k++)
				{
					double& temp = dx_times_Psi_Hzy_RightY(i,j,k);
					//RightY
					temp=Bpml_h_righty(j)*temp + Apml_h_righty(j)*(Ex(i,j+1,k)-Ex(i,j,k));
					Hz(i,j,k)+=Db_Z(i,j,k)*temp;
				}
			}
		}
}

void Cpml::UpdateE_Z()
{
		//dx_times_Psi_Exz_LowerZ
		for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
			for (j=Ex_min_index_in_y; j<=Ex_max_index_in_y; j++){
				for (k=Ex_min_index_in_z_lower; k<=Ex_max_index_in_z_lower; k++)
				{
					double& temp = dx_times_Psi_Exz_LowerZ(i,j,k);
					//LowerZ
					temp=Bpml_e_lowerz(k)*temp + Apml_e_lowerz(k)*(Hy(i,j,k)-Hy(i,j,k-1));
					Ex(i,j,k)+=Cb_X(i,j,k)*(- temp);
				}
			}
		}
		//dx_times_Psi_Eyz_LowerZ
		for (i=Ey_min_index_in_x; i<=Ey_max_index_in_x; i++){
			for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
				for (k=Ey_min_index_in_z_lower; k<=Ey_max_index_in_z_lower; k++)
				{
					double& temp = dx_times_Psi_Eyz_LowerZ(i,j,k);
					//LowerZ
					temp=Bpml_e_lowerz(k)*temp + Apml_e_lowerz(k)*(Hx(i,j,k)-Hx(i,j,k-1));
					Ey(i,j,k)+=Cb_Y(i,j,k)*temp;
				}
			}
		}

		//dx_times_Psi_Exz_UpperZ
		for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
			for (j=Ex_min_index_in_y; j<=Ex_max_index_in_y; j++){
				for (k=Ex_min_index_in_z_upper; k<=Ex_max_index_in_z_upper; k++)
				{
					double& temp = dx_times_Psi_Exz_UpperZ(i,j,k);
					//UpperZ
					temp=Bpml_e_upperz(k)*temp + Apml_e_upperz(k)*(Hy(i,j,k)-Hy(i,j,k-1));
					Ex(i,j,k)+=Cb_X(i,j,k)*(-temp);
				}
			}
		}
		//dx_times_Psi_Eyz_UpperZ
		for (i=Ey_min_index_in_x; i<=Ey_max_index_in_x; i++){
			for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
				for (k=Ey_min_index_in_z_upper; k<=Ey_max_index_in_z_upper; k++)
				{
					double& temp = dx_times_Psi_Eyz_UpperZ(i,j,k);
					//UpperZ
					temp=Bpml_e_upperz(k)*temp + Apml_e_upperz(k)*(Hx(i,j,k)-Hx(i,j,k-1));
					Ey(i,j,k)+=Cb_Y(i,j,k)*temp;
				}
			}
		}
}

void Cpml::UpdateH_Z()
{
		//dx_times_Psi_Hxz_LowerZ
		for (i=Hx_min_index_in_x; i<=Hx_max_index_in_x; i++){
			for (j=Hx_min_index_in_y; j<=Hx_max_index_in_y; j++){
				for (k=Hx_min_index_in_z_lower; k<=Hx_max_index_in_z_lower; k++)
				{
					double& temp = dx_times_Psi_Hxz_LowerZ(i,j,k);
					//LowerZ
					temp=Bpml_h_lowerz(k)*temp + Apml_h_lowerz(k)*(Ey(i,j,k+1)-Ey(i,j,k));
					Hx(i,j,k)+=Db_X(i,j,k)*temp;
				}
			}
		}
		//dx_times_Psi_Hyz_LowerZ
		for (i=Hy_min_index_in_x; i<=Hy_max_index_in_x; i++){
			for (j=Hy_min_index_in_y; j<=Hy_max_index_in_y; j++){
				for (k=Hy_min_index_in_z_lower; k<=Hy_max_index_in_z_lower; k++)
				{
					double& temp = dx_times_Psi_Hyz_LowerZ(i,j,k);
					//LowerZ
					temp=Bpml_h_lowerz(k)*temp + Apml_h_lowerz(k)*(Ex(i,j,k+1)-Ex(i,j,k));
					Hy(i,j,k)+=Db_Y(i,j,k)*(- temp);
				}
			}
		}

		//dx_times_Psi_Hxz_UpperZ
		for (i=Hx_min_index_in_x; i<=Hx_max_index_in_x; i++){
			for (j=Hx_min_index_in_y; j<=Hx_max_index_in_y; j++){
				for (k=Hx_min_index_in_z_upper; k<=Hx_max_index_in_z_upper; k++)
				{
					double& temp = dx_times_Psi_Hxz_UpperZ(i,j,k);
					//UpperZ
					temp=Bpml_h_upperz(k)*temp + Apml_h_upperz(k)*(Ey(i,j,k+1)-Ey(i,j,k));
					Hx(i,j,k)+=Db_X(i,j,k)*(temp);
				}
			}
		}
		//dx_times_Psi_Hyz_UpperZ
		for (i=Hy_min_index_in_x; i<=Hy_max_index_in_x; i++){
			for (j=Hy_min_index_in_y; j<=Hy_max_index_in_y; j++){
				for (k=Hy_min_index_in_z_upper; k<=Hy_max_index_in_z_upper; k++)
				{
					double& temp = dx_times_Psi_Hyz_UpperZ(i,j,k);
					//UpperZ
					temp=Bpml_h_upperz(k)*temp + Apml_h_upperz(k)*(Ex(i,j,k+1)-Ex(i,j,k));
					Hy(i,j,k)+=Db_Y(i,j,k)*(- temp);
				}
			}
		}
}
