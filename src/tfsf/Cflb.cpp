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

//Definition of the class "Cflb" for a TF/SF focused Hermite-Gaussian laser beam

#include "headers.h"

#include "Cflb.h"

#include "Cpw_fs.h"
#include "Cpw_2l.h"
#include "Cpw_ml.h"

//Uses TinyVector operations
#include <blitz/tinyvec-et.h>

//definition of Cwf needed
#include "waveforms/Cwf.h"

//#define ANGORA_CFLB_USE_CUBATURE
//#define ANGORA_CFLB_USE_HERMITIAN_CUBATURE
//#define ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
//#define ANGORA_CFLB_GAUSSIAN_QUAD_IN_PHI
//#define ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA

#ifdef ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
//gauss-legendre quadrature rule generator
extern void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w);
#endif

extern double dx;
extern int NCELLS_X,NCELLS_Y,NCELLS_Z;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;

extern int number_of_layers;

extern double Hermite(const int& N, const double& x);


Cflb::Cflb(const FLBDataType& MyData, const int& Index)
		:Data(MyData), FLBIndex(Index)
{
	/** incidence, rotation and polarization angles of the beam **/
	theta = Data.THETA;
	phi = Data.PHI;
	alpha = Data.alpha;
	psi = Data.PSI;

	if ((cos(Data.ap_half_angle)<0)&&(sin(Data.ap_half_angle)<0))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<double>(func_name,Data.ap_half_angle,
			"(should be between 0 and 90deg)");
	}

	if (cos(theta)>0)
	{//incident from the upper half-space
		epsilon_r_i = epsilon_r_upper;
		mu_r_i = mu_r_upper;
	}
	else
	{//incident from the lower half-space
		epsilon_r_i = epsilon_r_lower;
		mu_r_i = mu_r_lower;
	}

	n_i = sqrt(epsilon_r_i*mu_r_i);

	//maximum direction cosine in the incident ray bundle
	double s_max = sin(Data.ap_half_angle);

	// for lambda_min and lambda_max, -40 dB might be too tight, so we use -20dB
	double lambda_min = c/sqrt(epsilon_r_i*mu_r_i)/(Data.waveform->w_max_20()/2/M_PI);
	double lambda_max = c/sqrt(epsilon_r_i*mu_r_i)/(Data.waveform->w_min_20()/2/M_PI);

#ifndef ANGORA_CFLB_USE_CUBATURE

#ifdef ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
//	N_X1 = 16; //number of s values
//	N_X2 = 28; //number of phi values
//	N_X1 = 19; //number of s values  // ODD NUMBER CREATES USELESS POINTS AT THE ORIGIN!
//	N_X2 = 15; //number of phi values
//	N_X1 = 19; //number of s values  // ODD NUMBER CREATES USELESS POINTS AT THE ORIGIN!
//	N_X2 = 12; //number of phi values
	N_X1 = 20; //number of s values
	N_X2 = 8; //number of phi values
//	N_X1 = 28; //number of s values
//	N_X2 = 18; //number of phi values
//	N_X1 = 34; //number of s values
//	N_X2 = 24; //number of phi values

	double d_phi;
	Array<double,1> GLx,GLw;	//Gauss-Legendre quadrature points and weights between [-1,1]
#ifndef ANGORA_CFLB_GAUSSIAN_QUAD_IN_PHI
	//spacing of phi_pw (between [0,pi])
	d_phi = M_PI/N_X2;
#else
	Array<double,1> GL_phi_x,GL_phi_w;	//Gauss-Legendre quadrature points and weights between [-1,1]
	//Gauss-Legendre parameters for phi quadrature
	GL_phi_x.resize(N_X2); //Gauss-Legendre coordinates
	GL_phi_w.resize(N_X2); //Gauss-Legendre weights

	//calculate GL quadrature rule
	gaussquadrule(N_X2,GL_phi_x,GL_phi_w);

	GL_phi_x = M_PI/2*(1.0-GL_phi_x);	//-GLx because it is calculated in reverse direction in gaussquadrule
	GL_phi_w = M_PI/2*GL_phi_w;
#endif
	//Gauss-Legendre parameters for rho quadrature [ rho=sin(theta_pw), -sin(ap_half_angle)<rho<sin(ap_half_angle) ]
	GLx.resize(N_X1); //Gauss-Legendre coordinates
	GLw.resize(N_X1); //Gauss-Legendre weights

	//calculate GL quadrature rule
	gaussquadrule(N_X1,GLx,GLw);
	//rescale the positions and weights to the correct ranges of sin(theta_pw) (between [-sin(ap_half_angle),sin(ap_half_angle)])
	//see http://en.wikipedia.org/wiki/Gaussian_quadrature for the following formulas
	//GLx and GLw are for sin(theta_pw)
#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA
	GLx = abs(Data.ap_half_angle)*(-GLx);	//-GLx because it is calculated in reverse direction in gaussquadrule
	GLw = abs(Data.ap_half_angle)*GLw;
#else
	GLx = abs(s_max)*(-GLx);	//-GLx because it is calculated in reverse direction in gaussquadrule
	GLw = abs(s_max)*GLw;
#endif
#else //ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
	//approximate half-width (-60dB) of the focused beam for the maximum wavelength
	double max_beam_half_width = lambda_max/(M_PI*Data.filling_factor*s_max);

	//number of half-widths in the lateral extent of the beam
	double half_widths_in_beam = 5.2; // = 2x2.6  (beam is narrower without the factor of 2 in the denominator)
									  // (with the 2 factor, this would have been 7.4=2x3.7)

	//maximum diagonal extent of the TF/SF box (in m)
//	double TFSF_max_extent = sqrt(pow2((NCELLS_X-(Data.FLBMarginFrontX+Data.FLBMarginBackX))*dx)
//							 +pow2((NCELLS_Y-(Data.FLBMarginRightY+Data.FLBMarginLeftY))*dx)
//							 +pow2((NCELLS_Z-(Data.FLBMarginUpperZ+Data.FLBMarginLowerZ))*dx));
	/** This is specific to normal incidence. A more flexible way should be found later. **/
	double TFSF_max_extent = sqrt(pow2((NCELLS_X-(Data.FLBMarginFrontX+Data.FLBMarginBackX))*dx)
							 +pow2((NCELLS_Y-(Data.FLBMarginRightY+Data.FLBMarginLeftY))*dx));
	/** This is specific to normal incidence. A more flexible way should be found later. **/

	double TFSF_max_horizontal_extent = TFSF_max_extent;
	double TFSF_max_vertical_extent = TFSF_max_extent/2;

	//maximum lateral extent of the beam (ignoring the lateral extent of the TF/SF box)
	double max_extent_of_beam = max(
									max_beam_half_width*half_widths_in_beam,
								    4*lambda_max/s_max //the lower limit for an infinitely-overfilled back aperture
								   );
//								+TFSF_max_vertical_extent*tan(Data.ap_half_angle); //off-focus spread of the beam
	//minimum aliasing distance for negligible overlap
	double min_aliasing_displacement = max(max_extent_of_beam,TFSF_max_horizontal_extent);

	//uninitialized in config file, use sampling theorem for default spacing
	dsx = lambda_min/min_aliasing_displacement;
/////**REMOVE LATER **/
//dsx *=2;
/////**REMOVE LATER **/
	//The range 0->(dsx*N_X1) is divided into N_X1 regions and sx is placed at the midpoint of each region.
	//This provides consistency with previous publications.
	//The point (sx,sy) on the 2D plane may fall out of the sin(th)<1 circle for N_X1=N_X2=2, but neither N_X1 nor N_X2 are supposed to be that small anyway. Note that dsx<(lambda_center/beam_half_width)/(half_widths_in_beam)
	//						   <1/half_widths_in_beam
	// Therefore N_X1 is larger than (half_widths_in_beam), which is larger than 2.
	N_X1 = int(2*s_max/dsx)+1;

	//uninitialized in config file, use sampling theorem for default spacing
	dsy = lambda_min/min_aliasing_displacement;
///**REMOVE LATER **/
//dsy *=2;
///**REMOVE LATER **/
	//(see note above for sx)
	N_X2 = int(2*s_max/dsy)+1;
//cout << N_X1*N_X2 << endl;

#endif //ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
#endif //ANGORA_CFLB_USE_CUBATURE

	double sx,sy;

/////** REMOVE LATER **/
//int pwindex=0;
/////** REMOVE LATER **/
//ofstream temp;
//temp.open("tempfile",ios::out);
	double pw_factor; //extra amplitude factor applied to each plane wave (angle-dependent)
#ifdef ANGORA_CFLB_USE_CUBATURE

#include "cubature.cc"

#ifdef ANGORA_CFLB_USE_HERMITIAN_CUBATURE
		for (int cb_index=0;cb_index<N_cb;cb_index++)
		{
			sx = s_max*Data.filling_factor*cb_xi[cb_index];
			sy = s_max*Data.filling_factor*cb_yi[cb_index];
			/** determine the incidence and polarization angles of the individual plane wave **/
			//the plane wave is always incident from the upper half space
			double sz = real(sqrt((complex<double>)(1-pow2(sx)-pow2(sy)))); //cast into complex, in case there are roundoff errors for sx^2+sy^2=1

			theta_pw = acos(sz);	//sz = cos(theta_pw), result of acos is always in the [0,pi] range

			if ((abs(sx)<LIBSTD_DBL_EPSILON*100)&&(abs(sy)<LIBSTD_DBL_EPSILON*100)) //is theta_pw=0?
			{
				phi_pw = 0; //just assign an arbitrary value, since phi_pw is undefined here
			}
			else
			{
				phi_pw = atan2(sy,sx);
			}
			//atan2 is between [-pi,pi]
			//if less than 0, make phi_pw positive
			if (phi_pw<0)
			{
				phi_pw += 2*M_PI;
			}
//temp << sx << "," << sy << "," << pow2(s_max)*cb_weights[cb_index] << endl;
			/** incidence and polarization angles of the individual plane wave are determined. **/
			pw_factor = pow2(s_max*Data.filling_factor)*cb_weights[cb_index]/sz // = infinitesimal solid angle = dsx*dsy/sz
					*sqrt(Data.n_obj/n_i)*sqrt(abs(cos(theta_pw)))  //radiometric factor (Novotny "Principles of Nano-Optics" p. 57)
					*(n_i*Data.f1/(2*M_PI*c))  //other factors (from Richards&Wolf'59, Novotny "Principles of Nano-Optics" p. 60)
												  //note that k always belongs to the image space of the focusing lens
					//H_m(sqrt(2)x/w0)H_n(sqrt(2)y/w0)exp(-(x^2+y^2)/w0),
					// where w0=(filling_factor)*f1*s_max
					// x = f1*sin(theta)*cos(phi) = f1*sx
					// y = f1*sin(theta)*sin(phi) = f1*sy
					*Hermite(Data.x_order,(sqrt(2)/Data.filling_factor)*sx/s_max)  //NOTE: Exponential dependence is in the cubature!!!!
					*Hermite(Data.y_order,(sqrt(2)/Data.filling_factor)*sy/s_max); //NOTE: Exponential dependence is in the cubature!!!!
#else //ANGORA_CFLB_USE_HERMITIAN_CUBATURE
		for (int cb_index=0;cb_index<N_cb;cb_index++)
		{
			sx = s_max*cb_xi[cb_index];
			sy = s_max*cb_yi[cb_index];
			/** determine the incidence and polarization angles of the individual plane wave **/
			//the plane wave is always incident from the upper half space
			double sz = real(sqrt((complex<double>)(1-pow2(sx)-pow2(sy)))); //cast into complex, in case there are roundoff errors for sx^2+sy^2=1

			theta_pw = acos(sz);	//sz = cos(theta_pw), result of acos is always in the [0,pi] range

			if ((abs(sx)<LIBSTD_DBL_EPSILON*100)&&(abs(sy)<LIBSTD_DBL_EPSILON*100)) //is theta_pw=0?
			{
				phi_pw = 0; //just assign an arbitrary value, since phi_pw is undefined here
			}
			else
			{
				phi_pw = atan2(sy,sx);
			}
			//atan2 is between [-pi,pi]
			//if less than 0, make phi_pw positive
			if (phi_pw<0)
			{
				phi_pw += 2*M_PI;
			}
//temp << sx << "," << sy << "," << pow2(s_max)*cb_weights[cb_index] << endl;
			/** incidence and polarization angles of the individual plane wave are determined. **/
			pw_factor = pow2(s_max)*cb_weights[cb_index]/sz // = infinitesimal solid angle = dsx*dsy/sz
					*sqrt(Data.n_obj/n_i)*sqrt(abs(cos(theta_pw)))  //radiometric factor (Novotny "Principles of Nano-Optics" p. 57)
					*(n_i*Data.f1/(2*M_PI*c))  //other factors (from Richards&Wolf'59, Novotny "Principles of Nano-Optics" p. 60)
												  //note that k always belongs to the image space of the focusing lens
					//H_m(sqrt(2)x/w0)H_n(sqrt(2)y/w0)exp(-(x^2+y^2)/w0),
					// where w0=(filling_factor)*f1*s_max
					// x = f1*sin(theta)*cos(phi) = f1*sx
					// y = f1*sin(theta)*sin(phi) = f1*sy
					*Hermite(Data.x_order,(sqrt(2)/Data.filling_factor)*sx/s_max)*exp(-pow2(sx/s_max/Data.filling_factor))
					*Hermite(Data.y_order,(sqrt(2)/Data.filling_factor)*sy/s_max)*exp(-pow2(sy/s_max/Data.filling_factor));
#endif //ANGORA_CFLB_USE_HERMITIAN_CUBATURE

#else //ANGORA_CFLB_USE_CUBATURE
	//start adding the plane waves
	for (int i=0;i<N_X1;i++)
	{
		for (int j=0;j<N_X2;j++)
		{
#ifdef ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA
			theta_pw = GLx(i); //GLx and GLw are for theta_pw
#else
			theta_pw = asin(GLx(i)); //GLx and GLw are for sin(theta_pw), result is between [-pi/2,pi/2]
#endif

#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_PHI
			phi_pw = GL_phi_x(j); //GL_phi_x and GL_phi_w are for phi
#else
			phi_pw = (j+0.5)*d_phi;
#endif
			/** convert to traditional spherical angles **/
			//This is important, because the polarizations of each plane wave below are determined with respect to the "phi" angle,
			//which only works if "theta_pw" is between 0 and pi.
			if (sin(theta_pw)<0)
			{
				theta_pw = 2*M_PI-theta_pw;
				phi_pw += M_PI;
			}

			sx = sin(theta_pw)*cos(phi_pw);
			sy = sin(theta_pw)*sin(phi_pw);

			//amplitude factor for every plane wave
#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA
#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_PHI
			pw_factor = abs(sin(theta_pw))*GLw(i)*GL_phi_w(j)  //could also use GLx(i) for theta_pw
#else
			pw_factor = abs(sin(theta_pw))*GLw(i)*d_phi  //could also use GLx(i) for theta_pw
#endif
#else //ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA
#ifdef ANGORA_CFLB_GAUSSIAN_QUAD_IN_PHI
			pw_factor = abs(sin(theta_pw)/cos(theta_pw))*GLw(i)*GL_phi_w(j)  //could also use GLx(i) for sin(theta_pw)
#else
			pw_factor = abs(sin(theta_pw)/cos(theta_pw))*GLw(i)*d_phi  //could also use GLx(i) for sin(theta_pw)
#endif
#endif //ANGORA_CFLB_GAUSSIAN_QUAD_IN_THETA
										 *sqrt(Data.n_obj/n_i)*sqrt(abs(cos(theta_pw)))  //radiometric factor (Novotny "Principles of Nano-Optics" p. 57)
										 *(n_i*Data.f1/(2*M_PI*c))  //other factors (from Richards&Wolf'59, Novotny "Principles of Nano-Optics" p. 60)
												  //note that k always belongs to the image space of the focusing lens
										//H_m(sqrt(2)x/w0)H_n(sqrt(2)y/w0)exp(-(x^2+y^2)/w0),
										// where w0=(filling_factor)*f1*s_max
										// x = f1*sin(theta)*cos(phi) = f1*sx
										// y = f1*sin(theta)*sin(phi) = f1*sy
										*Hermite(Data.x_order,(sqrt(2)/Data.filling_factor)*sx/s_max)*exp(-pow2(sx/s_max/Data.filling_factor))
										*Hermite(Data.y_order,(sqrt(2)/Data.filling_factor)*sy/s_max)*exp(-pow2(sy/s_max/Data.filling_factor));
#else //ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE
			sx = dsx*(i-(N_X1-1.0)/2.0);
			sy = dsy*(j-(N_X2-1.0)/2.0);
			/** determine the incidence and polarization angles of the individual plane wave **/
			//the plane wave is always incident from the upper half space
			double sz = real(sqrt((complex<double>)(1-pow2(sx)-pow2(sy)))); //cast into complex, in case there are roundoff errors for sx^2+sy^2=1

			theta_pw = acos(sz);	//sz = cos(theta_pw), result of acos is always in the [0,pi] range

			if ((abs(sx)<LIBSTD_DBL_EPSILON*100)&&(abs(sy)<LIBSTD_DBL_EPSILON*100)) //is theta_pw=0?
			{
				phi_pw = 0; //just assign an arbitrary value, since phi_pw is undefined here
			}
			else
			{
				phi_pw = atan2(sy,sx);
			}
			//atan2 is between [-pi,pi]
			//if less than 0, make phi_pw positive
			if (phi_pw<0)
			{
				phi_pw += 2*M_PI;
			}

			/** incidence and polarization angles of the individual plane wave are determined. **/
			pw_factor = dsx*dsy/sz // = infinitesimal solid angle = sin(theta)*d_theta*d_phi
					*sqrt(Data.n_obj/n_i)*sqrt(abs(cos(theta_pw)))  //radiometric factor (Novotny "Principles of Nano-Optics" p. 57)
					*(n_i*Data.f1/(2*M_PI*c))  //other factors (from Richards&Wolf'59, Novotny "Principles of Nano-Optics" p. 60)
												  //note that k always belongs to the image space of the focusing lens
					//H_m(sqrt(2)x/w0)H_n(sqrt(2)y/w0)exp(-(x^2+y^2)/w0),
					// where w0=(filling_factor)*f1*s_max
					// x = f1*sin(theta)*cos(phi) = f1*sx
					// y = f1*sin(theta)*sin(phi) = f1*sy
					*Hermite(Data.x_order,(sqrt(2)/Data.filling_factor)*sx/s_max)*exp(-pow2(sx/s_max/Data.filling_factor))
					*Hermite(Data.y_order,(sqrt(2)/Data.filling_factor)*sy/s_max)*exp(-pow2(sy/s_max/Data.filling_factor));
#endif //ANGORA_CFLB_USE_GAUSSIAN_QUADRATURE

#endif //ANGORA_CFLB_USE_CUBATURE

			double polarization = psi;  //psi is defined like phi with respect to the local x axis
			//psi_pw is adjusted to give zero cross-pol component
			double angle_meridional = phi_pw - polarization; //angle with respect to the meridional plane
			// the angle with the meridional plane translates differently to psi_pw depending on the incidence half space
			if (cos(theta_pw)>=0) //incident from upper half space
			{
				psi_pw = 3*M_PI/2 - angle_meridional;
			}
			else //incident from lower half space
			{
				psi_pw = M_PI/2 + angle_meridional; //pw_pol is ccw wrt the +x axis, like the spherical phi_pw angle
			}

//			if ((pow2(sx)+pow2(sy))<=pow2(s_max))  //is the direction within the illumination cone?
			if (sin(theta_pw)<=s_max)  //is the direction within the illumination cone?
			{
//cout << pwindex++ << endl;
				/****************************************************************/
				/** Rotate the incidence directions and electric field vectors **/
				/****************************************************************/
				//ccw (w.r.t. z) azimuthal rotation angle from the global x axis to the local x axis of the beam
				double rot_az = ((sin(theta)>=0)?phi+M_PI/2+alpha:phi-M_PI/2+alpha);

				Array<double,2> rot_theta(3,3),rot_alpha(3,3),rot_tot(3,3);
				// Matrix representing counter-clockwise (right-hand) rotation w.r.t. the z axis by rot_az
				//(http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
				rot_alpha = cos(rot_az) , -sin(rot_az) , 0,
							sin(rot_az) ,  cos(rot_az) , 0,
									0  ,          0  , 1;

				double ux = -sin(phi);
				double uy = cos(phi);
				// Matrix representing counter-clockwise (right-hand) rotation w.r.t. the axis (ux,uy,0) by theta
				//(http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
				rot_theta = (cos(theta)+pow2(ux)*(1-cos(theta))) , ux*uy*(1-cos(theta))             , uy*sin(theta),
							ux*uy*(1-cos(theta))             , (cos(theta)+pow2(uy)*(1-cos(theta))) , -ux*sin(theta),
							-uy*sin(theta)                   , ux*sin(theta)                    , cos(theta);
				//net rotation matrix
				firstIndex iii;
				secondIndex jjj;
				thirdIndex kkk;
				//azimuthal rotation is done before the theta rotation
				rot_tot = sum(rot_theta(iii,kkk)*rot_alpha(kkk,jjj),kkk);

				TinyVector<double,3> rot_tot_row1,rot_tot_row2,rot_tot_row3;
				rot_tot_row1 = (rot_tot(0,0)), (rot_tot(0,1)), (rot_tot(0,2));
				rot_tot_row2 = (rot_tot(1,0)), (rot_tot(1,1)), (rot_tot(1,2));
				rot_tot_row3 = (rot_tot(2,0)), (rot_tot(2,1)), (rot_tot(2,2));

				TinyVector<double,3> dircos;
				dircos = (sin(theta_pw)*cos(phi_pw)) , (sin(theta_pw)*sin(phi_pw)) , (cos(theta_pw));
				TinyVector<double,3> k_inc = -dircos;
				// Note that sin(theta_pw) is always >=0, since theta_pw is always between 0 and pi
				TinyVector<double,3> k_inc_lateral;
				k_inc_lateral = (-cos(phi_pw)),(-sin(phi_pw)),(0.0);

				TinyVector<double,3> unit_z;
				unit_z = (0.0),(0.0),(1.0);
				TinyVector<double,3> k_E;
				k_E = cos(psi_pw)*cross(k_inc_lateral,unit_z)+sin(psi_pw)*cross(cross(k_inc_lateral,unit_z),k_inc);
				TinyVector<double,3>  E;
				E = (k_E(0)), (k_E(1)), (k_E(2));

				TinyVector<double,3> dircos_rot;
				dircos_rot = (dot(rot_tot_row1,dircos)), (dot(rot_tot_row2,dircos)), (dot(rot_tot_row3,dircos));
				TinyVector<double,3> E_rot;
				E_rot = (dot(rot_tot_row1,E)), (dot(rot_tot_row2,E)), (dot(rot_tot_row3,E));

				double theta_rot = acos(dircos_rot(thirdDim));
				double phi_rot = atan2(dircos_rot(secondDim),dircos_rot(firstDim));
				TinyVector<double,3> k_inc_rot,k_inc_rot_lateral;
				k_inc_rot = (-sin(theta_rot)*cos(phi_rot)) , (-sin(theta_rot)*sin(phi_rot)) , (-cos(theta_rot));
				if (sin(theta_rot)>=0)//theta_rot could be any wild direction
				{
					k_inc_rot_lateral = (-cos(phi_rot)),(-sin(phi_rot)),(0.0);
				}
				else
				{
					k_inc_rot_lateral = (cos(phi_rot)),(sin(phi_rot)),(0.0);
				}
				double psi_rot = atan2(dot(E_rot,cross(cross(k_inc_rot_lateral,unit_z),k_inc_rot)),dot(E_rot,cross(k_inc_rot_lateral,unit_z)));

				/*************************************************************/
				/** Incidence directions and electric field vectors rotated **/
				/*************************************************************/

				PWDataType PWData;
				PWData.THETA = theta_rot;
				PWData.PHI = phi_rot;
				PWData.PSI = psi_rot;

				PWData.PWMarginBackX = Data.FLBMarginBackX;
				PWData.PWMarginFrontX = Data.FLBMarginFrontX;
				PWData.PWMarginLeftY = Data.FLBMarginLeftY;
				PWData.PWMarginRightY = Data.FLBMarginRightY;
				PWData.PWMarginLowerZ = Data.FLBMarginLowerZ;
				PWData.PWMarginUpperZ = Data.FLBMarginUpperZ;
				PWData.PWOriginX = Data.FLBOriginX;
				PWData.PWOriginY = Data.FLBOriginY;
				PWData.PWOriginZ = Data.FLBOriginZ;
				PWData.L_req = Data.L_req;

				PWData.DisplayWarnings = Data.DisplayWarnings;

				PWData.E0 = Data.E0*pw_factor;

				//waveform of the constituent plane waves is the derivative of the incident plane-wave waveform
				PWData.waveform = Data.waveform->Derivative();

				//Add the plane wave source to the Hermite-Gaussian-beam object
				if (number_of_layers==1)
				{//attach a free-space plane-wave source
					boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_fs(PWData));
					PlaneWaves.push_back(new_pw_ptr);
				}
				else if (number_of_layers==2)
				{//attach a 2-layered-medium plane-wave source
					if (cos(theta)*cos(theta_rot)<ANGORA_CPW_ML_MIN_GRAZING_COS_Z)
					{//grazing angle too low, some plane waves fall within opposite half space
					 //this is not physical, so don't allow it
						throw FLBGrazingIncidenceException();
					}
					else
					{
						boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_2l(PWData));
						PlaneWaves.push_back(new_pw_ptr);
					}
				}
				else
				{//attach a multilayered-medium plane-wave source
					if (cos(theta)*cos(theta_rot)<ANGORA_CPW_ML_MIN_GRAZING_COS_Z)
					{//grazing angle too low, some plane waves fall within opposite half space
					 //this is not physical, so don't allow it
						throw FLBGrazingIncidenceException();
					}
					else
					{
						boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_ml(PWData));
						PlaneWaves.push_back(new_pw_ptr);
					}
				}
			}
		}
#ifndef ANGORA_CFLB_USE_CUBATURE
	}
#endif
//temp.close();
}

void Cflb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the focused laser beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Cflb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the focused laser beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Cflb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
}

void Cflb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Cflb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
}
