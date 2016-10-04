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

#ifndef CPML_H
#define CPML_H

//Declaration of the PML object "CPML"


class Cpml
{
 public:
	 Cpml(const int& my_pml_thickness, const double& SizeOfScatterer);
	 void UpdateE();		//Do the E-field updates in the PML region
	 void UpdateH();		//Do the H-field updates in the PML region
 private:
	 //PML Current Arrays
	 //FrontX : (NPML)x(NCELLS_Y+NPML)x(NCELLS_Z+NPML)	cells
	 //RightY : (NCELLS_X+NPML)x(NPML)x(NCELLS_Z+NPML)	cells
	 //UpperZ : (NCELLS_X+NPML)x(NCELLS_Y+NPML)xNPML	cells
	 Array<double,3>
		 dx_times_Psi_Eyx_BackX,dx_times_Psi_Ezx_BackX,dx_times_Psi_Eyx_FrontX,dx_times_Psi_Ezx_FrontX,dx_times_Psi_Exy_LeftY,dx_times_Psi_Ezy_LeftY,dx_times_Psi_Exy_RightY,dx_times_Psi_Ezy_RightY,dx_times_Psi_Exz_UpperZ,dx_times_Psi_Eyz_UpperZ,dx_times_Psi_Exz_LowerZ,dx_times_Psi_Eyz_LowerZ,
		 dx_times_Psi_Hyx_BackX,dx_times_Psi_Hzx_BackX,dx_times_Psi_Hyx_FrontX,dx_times_Psi_Hzx_FrontX,dx_times_Psi_Hxy_LeftY,dx_times_Psi_Hzy_LeftY,dx_times_Psi_Hxy_RightY,dx_times_Psi_Hzy_RightY,dx_times_Psi_Hxz_UpperZ,dx_times_Psi_Hyz_UpperZ,dx_times_Psi_Hxz_LowerZ,dx_times_Psi_Hyz_LowerZ;
	 //PML Parameters
	 int m;	//Order for polynomial grading
	 double epsilon_eff_x,epsilon_eff_y,epsilon_eff_lowerz,epsilon_eff_upperz;
	 double sigma_coeff,sigma_x_max,sigma_y_max,sigma_lowerz_max,sigma_upperz_max;
	 double kappa_x_max,kappa_y_max,kappa_z_max;
	 double alpha_x_max,alpha_y_max,alpha_z_max;
	 double alpha_x_min,alpha_y_min,alpha_z_min;
	 //PML Update Arrays (1-D)
	 Array<double,1> Apml_e_backx,Apml_e_frontx,Apml_e_lefty,Apml_e_righty,Apml_e_lowerz,Apml_e_upperz,
		 Bpml_e_backx,Bpml_e_frontx,Bpml_e_lefty,Bpml_e_righty,Bpml_e_lowerz,Bpml_e_upperz;
	 Array<double,1> Apml_h_backx,Apml_h_frontx,Apml_h_lefty,Apml_h_righty,Apml_h_lowerz,Apml_h_upperz,
		 Bpml_h_backx,Bpml_h_frontx,Bpml_h_lefty,Bpml_h_righty,Bpml_h_lowerz,Bpml_h_upperz;

	 const int pml_thickness;

	 int Ex_min_index_in_y_left,Ex_max_index_in_y_left,Ex_min_index_in_y_right,Ex_max_index_in_y_right;
	 int Ex_min_index_in_z_lower,Ex_max_index_in_z_lower,Ex_min_index_in_z_upper,Ex_max_index_in_z_upper;
	 int Ey_min_index_in_x_back,Ey_max_index_in_x_back,Ey_min_index_in_x_front,Ey_max_index_in_x_front;
	 int Ey_min_index_in_z_lower,Ey_max_index_in_z_lower,Ey_min_index_in_z_upper,Ey_max_index_in_z_upper;
	 int Ez_min_index_in_x_back,Ez_max_index_in_x_back,Ez_min_index_in_x_front,Ez_max_index_in_x_front;
	 int Ez_min_index_in_y_left,Ez_max_index_in_y_left,Ez_min_index_in_y_right,Ez_max_index_in_y_right;

	 int Hx_min_index_in_y_left,Hx_max_index_in_y_left,Hx_min_index_in_y_right,Hx_max_index_in_y_right;
	 int Hx_min_index_in_z_lower,Hx_max_index_in_z_lower,Hx_min_index_in_z_upper,Hx_max_index_in_z_upper;
	 int Hy_min_index_in_x_back,Hy_max_index_in_x_back,Hy_min_index_in_x_front,Hy_max_index_in_x_front;
	 int Hy_min_index_in_z_lower,Hy_max_index_in_z_lower,Hy_min_index_in_z_upper,Hy_max_index_in_z_upper;
	 int Hz_min_index_in_x_back,Hz_max_index_in_x_back,Hz_min_index_in_x_front,Hz_max_index_in_x_front;
	 int Hz_min_index_in_y_left,Hz_max_index_in_y_left,Hz_min_index_in_y_right,Hz_max_index_in_y_right;

	 void UpdateE_X();
	 void UpdateH_X();
	 void UpdateE_Y();
	 void UpdateH_Y();
	 void UpdateE_Z();
	 void UpdateH_Z();

};
#endif
