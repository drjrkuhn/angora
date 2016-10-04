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

#ifndef CFLB_H
#define CFLB_H

//Declaration of the class "Cflb" for a TF/SF focused Hermite-Gaussian laser beam

//only the declaration of Cpw needed: use forward declaration
class Cpw;

//base Angora exception class
#include "angora_excp.h"

//declaration of shared-pointers to Cwf objects
#include "waveforms/Cwf_shared_ptr.h"

//for the vector STL class
#include <vector>
//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

extern int OriginX,OriginY,OriginZ;


class FLBGrazingIncidenceException: public AngoraException
{// the exception class for beams cutting across planar layer interfaces
	public:
	virtual const string getError() const
	{
		return "Error in Cflb (focused laser-beam class): Laser-beam grazing angle too low. Try larger \"theta\" or narrower aperture.";
	}
};

struct FLBDataType
{
	FLBDataType():
		FLBMarginBackX(6),
		FLBMarginFrontX(6),
		FLBMarginLeftY(6),
		FLBMarginRightY(6),
		FLBMarginLowerZ(6),
		FLBMarginUpperZ(6),
		FLBOriginX(OriginX),
		FLBOriginY(OriginY),
		FLBOriginZ(OriginZ),
		L_req(15),

		DisplayWarnings(true),

		E0(1)
	{};

	double THETA,PHI,PSI;	//incidence angles
	double alpha;  //angle by which the beam is rotated around its propagation axis (clockwise, or left-handed)
	int x_order; //order of the Hermite polynomial in the x direction
	int y_order; //order of the Hermite polynomial in the y direction
	double ap_half_angle;	//maximum theta value of the incidence cone (in radians)
	double f1;	//back focal length (in m)
	double filling_factor;  // = (w0)/(f1*sin(theta_max)), where w0 is the beam width: exp(-(x^2+y^2)/w0^2)
	double n_obj; // object-side refractive index of the focusing lens

	//distances (in cells) of the TF/SF box from the PML boundary
	int FLBMarginBackX,FLBMarginFrontX,FLBMarginLeftY,FLBMarginRightY,FLBMarginLowerZ,FLBMarginUpperZ;
	//the coordinates of the origin of the TF/SF box //(measured from the rear-left-lower corner of the grid)
	double FLBOriginX,FLBOriginY,FLBOriginZ;
	double L_req;	//minimum number of grid cells required per minimum wavelength (checked against if warnings are enabled)

	bool DisplayWarnings;	//are warnings displayed when (grid cells/wavelength) is low?

	const_Cwf_shared_ptr waveform;	//pointer to the waveform object Cwf, representing the paraxial Hermite-Gaussian-beam waveform on the entrance pupil of the focusing lens. The Cwf object cannot be changed using this pointer.
	double E0;	//extra amplitude factor applied to the waveform on the entrance pupil of the focusing lens
};

class Cflb
{
 public:
	 Cflb(const FLBDataType& MyData, const int& Index);	//constructor

	 //Applies field corrections on the TF/SF box
	 void CorrectE(const int& n);
	 void CorrectH(const int& n);

	 int NumberOfPlaneWaves() const
	 {
		 return PlaneWaves.size();
	 };

 	 void WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const;	//cannot modify the Cflb object
	 void WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const; //cannot modify the Cflb object
	 void WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const;	//cannot modify the Cflb object

 private:
	 const int FLBIndex;	//index of the current beam in the TF/SF object Ctfsf

	 const FLBDataType Data;	//focused Hermite-Gaussian-beam parameters

	 //relative permittivity and permeability in the incidence half space (image space of the focusing lens)
	 double epsilon_r_i,mu_r_i;
	 //refractive index in the incidence half space (image space of the focusing lens)
	 double n_i;

	 //MIDPOINT QUAD. IN SX-SY: number of x and y direction cosines
	 //GAUSSIAN QUAD. IN S-PHI: number of s and phi values
	 int N_X1, N_X2;

	 double dsx,dsy;	//uniform direction-cosine spacings
	 double theta,phi,psi;	//incidence angles and the polarization angle for the beam
	 double alpha;  //angle by which the beam is rotated around its propagation axis (clockwise, or left-hand)
	 double theta_pw,phi_pw,psi_pw;	//incidence angles and the polarization angle for a single plane-wave component

	 //*********** Plane wave components ***************//
	 vector<boost::shared_ptr<Cpw> > PlaneWaves;	//array of pointers to plane wave objects
	 //*********** Plane wave components ***************//
};


#endif // CFLB_H
