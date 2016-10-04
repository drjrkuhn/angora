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

//Includes the routine that reads the total-field/scattered-field (TFSF) definitions

#include "headers.h"

#include "read_tfsf.h"

//definition of PWDataType needed
#include "Cpw.h"

//definition of FLBDataType needed
#include "Cflb.h"

////definition of FBDataType needed
//#include "Cfb.h"
//
////definition of HGBDataType needed
//#include "Chgb.h"
//
////definition of GSMBDataType needed
//#include "Cgsmb.h"
//
////definition of KBDataType needed
//#include "Ckohler.h"

//definition of Ctfsf needed
#include "Ctfsf.h"

//definition of Cwfs needed
#include "waveforms/Cwfs.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;

extern double dx;

extern int rank;

extern void add_slash_to_path(string& path);
extern int create_path(const string& path);


void read_tfsf(Ctfsf &TFSF, const Config& fdtdconfig, const Config& validsettings, const Cwfs& Waveforms)
{
	//Read TFSF settings
	string tfsf_setting_path = "TFSF";
	if (fdtdconfig.exists(tfsf_setting_path))
	{
		Setting& TFSFsettings = fdtdconfig.lookup(tfsf_setting_path);
		//check group for invalid settings
		CheckAngoraGroupSetting(TFSFsettings,validsettings);

		if (SettingEnabledForGrid(TFSFsettings))		//apply only if enabled for this grid
		{
			//Read plane-wave settings
			string pw_setting_path = "PlaneWaves";
			if (TFSFsettings.exists(pw_setting_path))
			{
				Setting& PWCollec = read_list_from_group(TFSFsettings,pw_setting_path);

				int num_of_pws = PWCollec.getLength();
				for (int pwindex=0; pwindex<num_of_pws; pwindex++)
				{
					Setting& PWsettings = PWCollec[pwindex];	//go to the pwindex'th plane-wave setting
					//check group for invalid settings
					CheckAngoraGroupSetting(PWsettings,validsettings);

					if (SettingEnabledForGrid(PWsettings))		//apply only if enabled for this grid
					{
						PWDataType PWData;	//call the constructor for each new plane wave

						read_value_from_group<double>(PWsettings,"theta",PWData.THETA);
						PWData.THETA *= M_PI/180;	//convert to radians

						read_value_from_group<double>(PWsettings,"phi",PWData.PHI);
						PWData.PHI *= M_PI/180;	//convert to radians

						read_value_from_group<double>(PWsettings,"psi",PWData.PSI);
						PWData.PSI *= M_PI/180;	//convert to radians

						double tfsf_back_margin_x,tfsf_front_margin_x,tfsf_left_margin_y,tfsf_right_margin_y,tfsf_lower_margin_z,tfsf_upper_margin_z;
						bool given_in_meters;
//						try{
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_back_margin_x",tfsf_back_margin_x,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginBackX = (int)round(tfsf_back_margin_x/dx);
							}
							else
							{
								PWData.PWMarginBackX = (int)round(tfsf_back_margin_x);
							}
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_front_margin_x",tfsf_front_margin_x,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginFrontX = (int)round(tfsf_front_margin_x/dx);
							}
							else
							{
								PWData.PWMarginFrontX = (int)round(tfsf_front_margin_x);
							}
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_left_margin_y",tfsf_left_margin_y,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginLeftY = (int)round(tfsf_left_margin_y/dx);
							}
							else
							{
								PWData.PWMarginLeftY = (int)round(tfsf_left_margin_y);
							}
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_right_margin_y",tfsf_right_margin_y,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginRightY = (int)round(tfsf_right_margin_y/dx);
							}
							else
							{
								PWData.PWMarginRightY = (int)round(tfsf_right_margin_y);
							}
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_lower_margin_z",tfsf_lower_margin_z,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginLowerZ = (int)round(tfsf_lower_margin_z/dx);
							}
							else
							{
								PWData.PWMarginLowerZ = (int)round(tfsf_lower_margin_z);
							}
							if (read_optional_length_from_group<double>(PWsettings,"tfsf_upper_margin_z",tfsf_upper_margin_z,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWMarginUpperZ = (int)round(tfsf_upper_margin_z/dx);
							}
							else
							{
								PWData.PWMarginUpperZ = (int)round(tfsf_upper_margin_z);
							}
//						}
//						/** TODO: Better organize the defaults in PWData **/
//						catch (AngoraSettingNotFoundException&)
//						{//default values
//						}

						double pw_origin_x,pw_origin_y,pw_origin_z;
//						try{
							if (read_optional_length_from_group<double>(PWsettings,"pw_origin_x",pw_origin_x,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWOriginX = pw_origin_x/dx + OriginX;  //these are grid cell indices, not coordinates
							}
							else
							{
								PWData.PWOriginX = pw_origin_x + OriginX;  //these are grid cell indices, not coordinates
							}
							if (read_optional_length_from_group<double>(PWsettings,"pw_origin_y",pw_origin_y,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWOriginY = pw_origin_y/dx + OriginY;  //these are grid cell indices, not coordinates
							}
							else
							{
								PWData.PWOriginY = pw_origin_y + OriginY;  //these are grid cell indices, not coordinates
							}
							if (read_optional_length_from_group<double>(PWsettings,"pw_origin_z",pw_origin_z,given_in_meters))
							if (given_in_meters)
							{
								PWData.PWOriginZ = pw_origin_z/dx + OriginZ;  //these are grid cell indices, not coordinates
							}
							else
							{
								PWData.PWOriginZ = pw_origin_z + OriginZ;  //these are grid cell indices, not coordinates
							}
//						}
//						/** TODO: Better organize the defaults in PWDataType **/
//						catch (AngoraSettingNotFoundException&)
//						{//default values
//						}

						read_optional_value_from_group<double>(PWsettings,"min_cells_per_lambda",PWData.L_req);
						read_optional_value_from_group<bool>(PWsettings,"display_warnings",PWData.DisplayWarnings);

						//read the extra amplitude factor
						read_optional_value_from_group<double>(PWsettings,"pw_extra_amplitude",PWData.E0);

						//read the index of the plane-wave waveform in the Waveforms object
						string waveform_tag;
						read_value_from_group<string>(PWsettings,"waveform_tag",waveform_tag);
						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
						//If not in check mode, add the plane wave to the TFSF object
						if (!check_mode)
						{
							//the pointer to the right waveform object
							PWData.waveform = waveformptr;

							//Add the plane wave to the TFSF object
							TFSF.AddPlaneWave(PWData);
						}
					}
				}
			}

			//Read focused Hermite-Gaussian laser-beam settings
			string flb_setting_path = "FocusedLaserBeams";
			if (TFSFsettings.exists(flb_setting_path))
			{
				Setting& FLBCollec = read_list_from_group(TFSFsettings,flb_setting_path);

				int num_of_flbs = FLBCollec.getLength();
				for (int flbindex=0; flbindex<num_of_flbs; flbindex++)
				{
					Setting& FLBsettings = FLBCollec[flbindex];	//go to the flbindex'th focused-beam setting
					//check group for invalid settings
					CheckAngoraGroupSetting(FLBsettings,validsettings);

					if (SettingEnabledForGrid(FLBsettings))		//apply only if enabled for this grid
					{
						FLBDataType FLBData;	//call the constructor for each new focused beam

						read_value_from_group<double>(FLBsettings,"theta",FLBData.THETA);
						FLBData.THETA *= M_PI/180;	//convert to radians

						read_value_from_group<double>(FLBsettings,"phi",FLBData.PHI);
						FLBData.PHI *= M_PI/180;	//convert to radians

						read_value_from_group<double>(FLBsettings,"psi",FLBData.PSI);
						FLBData.PSI *= M_PI/180;	//convert to radians

						if (!read_optional_value_from_group<double>(FLBsettings,"alpha",FLBData.alpha))
						{
							FLBData.alpha = 0;
						}
						FLBData.alpha *= M_PI/180;	//convert to radians

						read_value_from_group<int>(FLBsettings,"x_order",FLBData.x_order);
						read_value_from_group<int>(FLBsettings,"y_order",FLBData.y_order);

						read_value_from_group<double>(FLBsettings,"ap_half_angle",FLBData.ap_half_angle);
						FLBData.ap_half_angle *= M_PI/180;	//convert to radians

						if (!read_length_from_group<double>(FLBsettings,"back_focal_length",FLBData.f1))
						{//returned false: given in grid cells
							FLBData.f1 *= dx;
						}

						read_value_from_group<double>(FLBsettings,"filling_factor",FLBData.filling_factor);

						if (!read_optional_value_from_group<double>(FLBsettings,"object_space_refr_index",FLBData.n_obj))
						{
							FLBData.n_obj = 1.0;
						}

						double tfsf_back_margin_x,tfsf_front_margin_x,tfsf_left_margin_y,tfsf_right_margin_y,tfsf_lower_margin_z,tfsf_upper_margin_z;
						bool given_in_meters;
//						try{
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_back_margin_x",tfsf_back_margin_x,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginBackX = (int)round(tfsf_back_margin_x/dx);
							}
							else
							{
								FLBData.FLBMarginBackX = (int)round(tfsf_back_margin_x);
							}
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_front_margin_x",tfsf_front_margin_x,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginFrontX = (int)round(tfsf_front_margin_x/dx);
							}
							else
							{
								FLBData.FLBMarginFrontX = (int)round(tfsf_front_margin_x);
							}
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_left_margin_y",tfsf_left_margin_y,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginLeftY = (int)round(tfsf_left_margin_y/dx);
							}
							else
							{
								FLBData.FLBMarginLeftY = (int)round(tfsf_left_margin_y);
							}
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_right_margin_y",tfsf_right_margin_y,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginRightY = (int)round(tfsf_right_margin_y/dx);
							}
							else
							{
								FLBData.FLBMarginRightY = (int)round(tfsf_right_margin_y);
							}
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_lower_margin_z",tfsf_lower_margin_z,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginLowerZ = (int)round(tfsf_lower_margin_z/dx);
							}
							else
							{
								FLBData.FLBMarginLowerZ = (int)round(tfsf_lower_margin_z);
							}
							if (read_optional_length_from_group<double>(FLBsettings,"tfsf_upper_margin_z",tfsf_upper_margin_z,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBMarginUpperZ = (int)round(tfsf_upper_margin_z/dx);
							}
							else
							{
								FLBData.FLBMarginUpperZ = (int)round(tfsf_upper_margin_z);
							}
//						}
//						/** TODO: Better organize the defaults in FLBData **/
//						catch (AngoraSettingNotFoundException&)
//						{//default values
//						}

						double flb_origin_x,flb_origin_y,flb_origin_z;
//						try{
							if (read_optional_length_from_group<double>(FLBsettings,"flb_origin_x",flb_origin_x,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBOriginX = flb_origin_x/dx + OriginX;  //these are grid cell indices, not coordinates
							}
							else
							{
								FLBData.FLBOriginX = flb_origin_x + OriginX;  //these are grid cell indices, not coordinates
							}
							if (read_optional_length_from_group<double>(FLBsettings,"flb_origin_y",flb_origin_y,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBOriginY = flb_origin_y/dx + OriginY;  //these are grid cell indices, not coordinates
							}
							else
							{
								FLBData.FLBOriginY = flb_origin_y + OriginY;  //these are grid cell indices, not coordinates
							}
							if (read_optional_length_from_group<double>(FLBsettings,"flb_origin_z",flb_origin_z,given_in_meters))
							if (given_in_meters)
							{
								FLBData.FLBOriginZ = flb_origin_z/dx + OriginZ;  //these are grid cell indices, not coordinates
							}
							else
							{
								FLBData.FLBOriginZ = flb_origin_z + OriginZ;  //these are grid cell indices, not coordinates
							}
//						}
//						/** TODO: Better organize the defaults in FLBDataType **/
//						catch (AngoraSettingNotFoundException&)
//						{//default values
//						}

						read_optional_value_from_group<double>(FLBsettings,"min_cells_per_lambda",FLBData.L_req);
						read_optional_value_from_group<bool>(FLBsettings,"display_warnings",FLBData.DisplayWarnings);

						//read the extra amplitude factor
						read_optional_value_from_group<double>(FLBsettings,"flb_extra_amplitude",FLBData.E0);

						//read the index of the focused-beam waveform in the Waveforms object
						string waveform_tag;
						read_value_from_group<string>(FLBsettings,"waveform_tag",waveform_tag);
						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
						//If not in check mode, add the focused beam to the TF/SF object
						if (!check_mode)
						{
							//the pointer to the right waveform object
							FLBData.waveform = waveformptr;

							//Add a focused beam to the TF/SF object
							TFSF.AddFocusedLaserBeam(FLBData);
						}
					}
				}
			}

//			//Read focused-beam settings
//			string fp_setting_path = "FocusedBeams";
//			if (TFSFsettings.exists(fp_setting_path))
//			{
//				Setting& FBCollec = read_list_from_group(TFSFsettings,fp_setting_path);
//
//				int num_of_fps = FBCollec.getLength();
//				for (int fpindex=0; fpindex<num_of_fps; fpindex++)
//				{
//					Setting& FBsettings = FBCollec[fpindex];	//go to the fpindex'th focused-beam setting
//					//check group for invalid settings
//					CheckAngoraGroupSetting(FBsettings,validsettings);
//
//					if (SettingEnabledForGrid(FBsettings))		//apply only if enabled for this grid
//					{
//						FBDataType FBData;	//call the constructor for each new focused beam
//
//						read_value_from_group<double>(FBsettings,"theta",FBData.THETA);
//						FBData.THETA *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(FBsettings,"phi",FBData.PHI);
//						FBData.PHI *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(FBsettings,"psi",FBData.PSI);
//						FBData.PSI *= M_PI/180;	//convert to radians
//
//						if(!read_optional_value_from_group<string>(FBsettings,"angular_discretization",FBData.angular_discretization))
//						{
//							FBData.angular_discretization = "cartesian";
//						}
//						else if ((FBData.angular_discretization!="cartesian")&&(FBData.angular_discretization!="radial"))
//						{
//							throw AngoraInvalidSettingValueException(FBsettings["angular_discretization"],"should be either \"cartesian\" or \"radial\"");
//						}
//
//						read_value_from_group<double>(FBsettings,"ap_half_angle",FBData.ap_half_angle);
//						FBData.ap_half_angle *= M_PI/180;	//convert to radians
//
//						if (FBData.angular_discretization=="radial")
//						{
//							read_value_from_group<int>(FBsettings,"n_1",FBData.n_1);
//							read_value_from_group<int>(FBsettings,"n_2",FBData.n_2);
//						}
//						else
//						{
//							//then these are optional (left at their default sentinel values)
//							read_optional_value_from_group<int>(FBsettings,"n_1",FBData.n_1);
//							read_optional_value_from_group<int>(FBsettings,"n_2",FBData.n_2);
//						}
//
//						read_value_from_group<double>(FBsettings,"focal_length",FBData.f);
//
//						read_optional_value_from_group<int>(FBsettings,"tfsf_back_margin_x_in_cells",FBData.FBMarginBackX);
//						read_optional_value_from_group<int>(FBsettings,"tfsf_front_margin_x_in_cells",FBData.FBMarginFrontX);
//						read_optional_value_from_group<int>(FBsettings,"tfsf_left_margin_y_in_cells",FBData.FBMarginLeftY);
//						read_optional_value_from_group<int>(FBsettings,"tfsf_right_margin_y_in_cells",FBData.FBMarginRightY);
//						read_optional_value_from_group<int>(FBsettings,"tfsf_lower_margin_z_in_cells",FBData.FBMarginLowerZ);
//						read_optional_value_from_group<int>(FBsettings,"tfsf_upper_margin_z_in_cells",FBData.FBMarginUpperZ);
//
//						if (read_optional_value_from_group<double>(FBsettings,"fb_origin_x_in_cells",FBData.FBOriginX))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							FBData.FBOriginX += OriginX;
//						}
//
//						if (read_optional_value_from_group<double>(FBsettings,"fb_origin_y_in_cells",FBData.FBOriginY))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							FBData.FBOriginY += OriginY;
//						}
//
//						if (read_optional_value_from_group<double>(FBsettings,"fb_origin_z_in_cells",FBData.FBOriginZ))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							FBData.FBOriginZ += OriginZ;
//						}
//
//						read_optional_value_from_group<double>(FBsettings,"min_cells_per_lambda",FBData.L_req);
//						read_optional_value_from_group<bool>(FBsettings,"display_warnings",FBData.DisplayWarnings);
//
//						//read the extra amplitude factor
//						read_optional_value_from_group<double>(FBsettings,"fb_extra_amplitude",FBData.E0);
//
//						//read the index of the focused-beam waveform in the Waveforms object
//						string waveform_tag;
//						read_value_from_group<string>(FBsettings,"waveform_tag",waveform_tag);
//						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
//						//If not in check mode, add the focused beam to the TF/SF object
//						if (!check_mode)
//						{
//							//the pointer to the right waveform object
//							FBData.waveform = waveformptr;
//
//							//Add a focused beam to the TF/SF object
//							TFSF.AddFocusedBeam(FBData);
//						}
//					}
//				}
//			}
//
//			//Read Hermite-Gaussian-beam settings
//			string hgp_setting_path = "HermiteGaussianBeams";
//			if (TFSFsettings.exists(hgp_setting_path))
//			{
//				Setting& HGBCollec = read_list_from_group(TFSFsettings,hgp_setting_path);
//
//				int num_of_hgps = HGBCollec.getLength();
//				for (int hgpindex=0; hgpindex<num_of_hgps; hgpindex++)
//				{
//					Setting& HGBsettings = HGBCollec[hgpindex];	//go to the hgpindex'th Hermite-Gaussian-beam setting
//					//check group for invalid settings
//					CheckAngoraGroupSetting(HGBsettings,validsettings);
//
//					if (SettingEnabledForGrid(HGBsettings))		//apply only if enabled for this grid
//					{
//						HGBDataType HGBData;	//call the constructor for each new Hermite-Gaussian beam
//
//						read_value_from_group<double>(HGBsettings,"theta",HGBData.THETA);
//						HGBData.THETA *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(HGBsettings,"phi",HGBData.PHI);
//						HGBData.PHI *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(HGBsettings,"psi",HGBData.PSI);
//						HGBData.PSI *= M_PI/180;	//convert to radians
//
//						if (!read_optional_value_from_group<double>(HGBsettings,"alpha",HGBData.alpha))
//						{
//							HGBData.alpha = 0;
//						}
//						HGBData.alpha *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(HGBsettings,"beam_half_width",HGBData.beam_half_width);
//
//						read_value_from_group<int>(HGBsettings,"x_order",HGBData.x_order);
//						read_value_from_group<int>(HGBsettings,"y_order",HGBData.y_order);
//
//	//					read_value_from_group<double>(HGBsettings,"pol_angle",HGBData.polarization);
//	//					HGBData.polarization *= M_PI/180;	//convert to radians
//
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_back_margin_x_in_cells",HGBData.HGBMarginBackX);
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_front_margin_x_in_cells",HGBData.HGBMarginFrontX);
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_left_margin_y_in_cells",HGBData.HGBMarginLeftY);
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_right_margin_y_in_cells",HGBData.HGBMarginRightY);
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_lower_margin_z_in_cells",HGBData.HGBMarginLowerZ);
//						read_optional_value_from_group<int>(HGBsettings,"tfsf_upper_margin_z_in_cells",HGBData.HGBMarginUpperZ);
//
//						if (read_optional_value_from_group<double>(HGBsettings,"hgb_origin_x_in_cells",HGBData.HGBOriginX))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							HGBData.HGBOriginX += OriginX;
//						}
//
//						if (read_optional_value_from_group<double>(HGBsettings,"hgb_origin_y_in_cells",HGBData.HGBOriginY))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							HGBData.HGBOriginY += OriginY;
//						}
//
//						if (read_optional_value_from_group<double>(HGBsettings,"hgb_origin_z_in_cells",HGBData.HGBOriginZ))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							HGBData.HGBOriginZ += OriginZ;
//						}
//
//						read_optional_value_from_group<double>(HGBsettings,"min_cells_per_lambda",HGBData.L_req);
//						read_optional_value_from_group<bool>(HGBsettings,"display_warnings",HGBData.DisplayWarnings);
//
//						//read the extra amplitude factor
//						read_optional_value_from_group<double>(HGBsettings,"hgb_extra_amplitude",HGBData.E0);
//
//						//read the index of the Hermite-Gaussian-beam waveform in the Waveforms object
//						string waveform_tag;
//						read_value_from_group<string>(HGBsettings,"waveform_tag",waveform_tag);
//						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
//						//If not in check mode, add the Hermite-Gaussian beam to the TF/SF object
//						if (!check_mode)
//						{
//							//the pointer to the right waveform object
//							HGBData.waveform = waveformptr;
//
//							//Add a Hermite-Gaussian beam to the TF/SF object
//							TFSF.AddHermiteGaussianBeam(HGBData);
//						}
//					}
//				}
//			}
//
//			//Read Gaussian-Schell-model-beam settings
//			string gsmp_setting_path = "GaussianSchellModelBeams";
//			if (TFSFsettings.exists(gsmp_setting_path))
//			{
//				Setting& GSMBCollec = read_list_from_group(TFSFsettings,gsmp_setting_path);
//
//				int num_of_gsmps = GSMBCollec.getLength();
//				for (int gsmpindex=0; gsmpindex<num_of_gsmps; gsmpindex++)
//				{
//					Setting& GSMBsettings = GSMBCollec[gsmpindex];	//go to the gsmpindex'th Gaussian-Schell-model-beam setting
//					//check group for invalid settings
//					CheckAngoraGroupSetting(GSMBsettings,validsettings);
//
//					if (SettingEnabledForGrid(GSMBsettings))		//apply only if enabled for this grid
//					{
//						GSMBDataType GSMBData;	//call the constructor for each new Gaussian Schell-model beam
//
//						read_value_from_group<double>(GSMBsettings,"theta",GSMBData.THETA);
//						GSMBData.THETA *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(GSMBsettings,"phi",GSMBData.PHI);
//						GSMBData.PHI *= M_PI/180;	//convert to radians
//
//						read_value_from_group<double>(GSMBsettings,"psi",GSMBData.PSI);
//						GSMBData.PSI *= M_PI/180;	//convert to radians
//
//						read_value_from_group<int>(GSMBsettings,"n_mode_x",GSMBData.n_mode_x);
//						read_value_from_group<int>(GSMBsettings,"n_mode_y",GSMBData.n_mode_y);
//
//						read_value_from_group<double>(GSMBsettings,"beam_half_width",GSMBData.beam_half_width);
//
//						read_value_from_group<double>(GSMBsettings,"correlation_length",GSMBData.correlation_length);
//
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_back_margin_x_in_cells",GSMBData.GSMBMarginBackX);
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_front_margin_x_in_cells",GSMBData.GSMBMarginFrontX);
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_left_margin_y_in_cells",GSMBData.GSMBMarginLeftY);
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_right_margin_y_in_cells",GSMBData.GSMBMarginRightY);
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_lower_margin_z_in_cells",GSMBData.GSMBMarginLowerZ);
//						read_optional_value_from_group<int>(GSMBsettings,"tfsf_upper_margin_z_in_cells",GSMBData.GSMBMarginUpperZ);
//
//						if (read_optional_value_from_group<double>(GSMBsettings,"gsmb_origin_x_in_cells",GSMBData.GSMBOriginX))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							GSMBData.GSMBOriginX += OriginX;
//						}
//
//						if (read_optional_value_from_group<double>(GSMBsettings,"gsmb_origin_y_in_cells",GSMBData.GSMBOriginY))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							GSMBData.GSMBOriginY += OriginY;
//						}
//
//						if (read_optional_value_from_group<double>(GSMBsettings,"gsmb_origin_z_in_cells",GSMBData.GSMBOriginZ))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							GSMBData.GSMBOriginZ += OriginZ;
//						}
//
//						read_optional_value_from_group<double>(GSMBsettings,"min_cells_per_lambda",GSMBData.L_req);
//						read_optional_value_from_group<bool>(GSMBsettings,"display_warnings",GSMBData.DisplayWarnings);
//
//						//read the extra amplitude factor
//						read_optional_value_from_group<double>(GSMBsettings,"gsmb_extra_amplitude",GSMBData.E0);
//
//						//read the index of the Gaussian-Schell-model-beam waveform in the Waveforms object
//						string waveform_tag;
//						read_value_from_group<string>(GSMBsettings,"waveform_tag",waveform_tag);
//						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
//						//If not in check mode, add the Gaussian Schell-model beam to the TF/SF object
//						if (!check_mode)
//						{
//							//the pointer to the right waveform object
//							GSMBData.waveform = waveformptr;
//
//							//Add a Gaussian Schell-model beam to the TF/SF object
//							TFSF.AddGaussianSchellModelBeam(GSMBData);
//						}
//					}
//				}
//			}
//
//			//Read Kohler-beam settings
//			string ki_setting_path = "KohlerBeams";
//			if (TFSFsettings.exists(ki_setting_path))
//			{
//				Setting& KBCollec = read_list_from_group(TFSFsettings,ki_setting_path);
//
//				int num_of_kis = KBCollec.getLength();
//				for (int kiindex=0; kiindex<num_of_kis; kiindex++)
//				{
//					Setting& KBsettings = KBCollec[kiindex];	//go to the kiindex'th Kohler-beam setting
//					//check group for invalid settings
//					CheckAngoraGroupSetting(KBsettings,validsettings);
//
//					if (SettingEnabledForGrid(KBsettings))	//apply only if enabled for this grid
//					{
//						KBDataType KBData;	//create the Kohler-beam data structure
//
//	//					int random_seed;	//common random seed among nodes
//	//					read_optional_value_from_group<int>(KBsettings,"random_seed",KBData.random_seed);
//	//
//	//					read_value_from_group<double>(KBsettings,"ap_half_angle",KBData.ap_half_angle);
//	//					KBData.ap_half_angle *= M_PI/180;	//convert to radians
//	//
//	//					read_value_from_group<int>(KBsettings,"n_theta",KBData.n_theta);
//	//					read_value_from_group<int>(KBsettings,"n_phi",KBData.n_phi);
//	//
//	//					read_value_from_group<double>(KBsettings,"f",KBData.f);
//
//						if(!read_optional_value_from_group<string>(KBsettings,"simulation_type",KBData.simulation_type))
//						{
//							KBData.simulation_type = "deterministic";
//						}
//						else if ((KBData.simulation_type!="deterministic")&&(KBData.simulation_type!="stochastic"))
//						{
//							throw AngoraInvalidSettingValueException(KBsettings["simulation_type"],"should be either \"deterministic\" or \"stochastic\"");
//						}
//
//						if(!read_optional_value_from_group<string>(KBsettings,"angular_discretization",KBData.angular_discretization))
//						{
//							KBData.angular_discretization = "cartesian";
//						}
//						else if ((KBData.angular_discretization!="cartesian")&&(KBData.angular_discretization!="radial"))
//						{
//							throw AngoraInvalidSettingValueException(KBsettings["angular_discretization"],"should be either \"cartesian\" or \"radial\"");
//						}
//
//						read_value_from_group<double>(KBsettings,"ap_half_angle",KBData.ap_half_angle);
//						KBData.ap_half_angle *= M_PI/180;	//convert to radians
//
//						if (KBData.angular_discretization=="radial")
//						{
//							read_value_from_group<int>(KBsettings,"n_1",KBData.n_1);
//							read_value_from_group<int>(KBsettings,"n_2",KBData.n_2);
//						}
//						else
//						{
//							//then these are optional (left at their default sentinel values)
//							read_optional_value_from_group<int>(KBsettings,"n_1",KBData.n_1);
//							read_optional_value_from_group<int>(KBsettings,"n_2",KBData.n_2);
//						}
//
//	//					read_value_from_group<double>(KBsettings,"f",KBData.f);
//
//						read_value_from_group<double>(KBsettings,"pw_pol_angle",KBData.pw_pol);
//						KBData.pw_pol *= M_PI/180;	//convert to radians
//
//						read_optional_value_from_group<int>(KBsettings,"tfsf_back_margin_x_in_cells",KBData.KBMarginBackX);
//						read_optional_value_from_group<int>(KBsettings,"tfsf_front_margin_x_in_cells",KBData.KBMarginFrontX);
//						read_optional_value_from_group<int>(KBsettings,"tfsf_left_margin_y_in_cells",KBData.KBMarginLeftY);
//						read_optional_value_from_group<int>(KBsettings,"tfsf_right_margin_y_in_cells",KBData.KBMarginRightY);
//						read_optional_value_from_group<int>(KBsettings,"tfsf_lower_margin_z_in_cells",KBData.KBMarginLowerZ);
//						read_optional_value_from_group<int>(KBsettings,"tfsf_upper_margin_z_in_cells",KBData.KBMarginUpperZ);
//
//						if (read_optional_value_from_group<double>(KBsettings,"kb_origin_x_in_cells",KBData.KBOriginX))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							KBData.KBOriginX += OriginX;
//						}
//
//						if (read_optional_value_from_group<double>(KBsettings,"kb_origin_y_in_cells",KBData.KBOriginY))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							KBData.KBOriginY += OriginY;
//						}
//
//						if (read_optional_value_from_group<double>(KBsettings,"kb_origin_z_in_cells",KBData.KBOriginZ))
//						{
//							//if given in the config file, it is with respect to the grid origin
//							KBData.KBOriginZ += OriginZ;
//						}
//
//						read_optional_value_from_group<double>(KBsettings,"min_cells_per_lambda",KBData.L_req);
//						read_optional_value_from_group<bool>(KBsettings,"display_warnings",KBData.DisplayWarnings);
//
//						//read the extra amplitude factor
//						read_optional_value_from_group<double>(KBsettings,"kb_extra_amplitude",KBData.E0);
//
//						//read the index of the Kohler-beam waveform in the Waveforms object
//						string waveform_tag;
//						read_value_from_group<string>(KBsettings,"waveform_tag",waveform_tag);
//						const_Cwf_shared_ptr waveformptr = Waveforms[waveform_tag];
//
//						//If not in check mode, add the Kohler beam to the TF/SF object
//						if (!check_mode)
//						{
//							//the pointer to the right waveform object
//							KBData.waveform = waveformptr;
//
//							//Add a Kohler beam to the TF/SF object
//							TFSF.AddKohlerBeam(KBData);
//						}
//					}
//				}
//			}
		}
	}
}
