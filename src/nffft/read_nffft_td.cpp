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

//Includes the routine that reads the TIME_DOMAIN near-field-to-far-field transform (NFFFT) definitions

#include "headers.h"

#include "read_nffft_td.h"

//definition of Cnffft_td needed
#include "td/Cnffft_td.h"

//definition of TrDataType_td needed
#include "td/Ctr_td.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;
extern string OutputDir;
extern string TimeDomainNFFFTOutputDir;
extern const string default_TimeDomainNFFFTOutputDir;
extern const string default_TimeDomainNFFFT_filename;
extern const string default_td_nffft_fileextension;

extern int OriginX,OriginY,OriginZ;

extern double dx;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);
extern int create_path(const string& path);


void read_nffft_td(Cnffft_td &NFFFT_td, const Config& fdtdconfig, const Config& validsettings, const Cpointsources &PointSources)
{
	//Read NFFFT output directory name
	try{read_value_from_group<string>(fdtdconfig.getRoot(),"td_nffft_output_dir",TimeDomainNFFFTOutputDir);}
	catch (AngoraSettingNotFoundException&)
	{
		TimeDomainNFFFTOutputDir = default_TimeDomainNFFFTOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(TimeDomainNFFFTOutputDir);
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(TimeDomainNFFFTOutputDir))
	{
		TimeDomainNFFFTOutputDir = OutputDir + TimeDomainNFFFTOutputDir;	//TimeDomainNFFFTOutputDir is relative to the output directory
	}
	//create NFFFT output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(TimeDomainNFFFTOutputDir)<0)
		{
			/** throw exception **/
			if (rank==0) cout << "Could not create path " << TimeDomainNFFFTOutputDir << endl;
		}
	}

	//Read time-domain NFFFT settings
	string nffft_td_setting_path = "TimeDomainNFFFT";
	if (fdtdconfig.exists(nffft_td_setting_path))
	{
		Setting& TimeDomainNFFFTsettings = read_list_from_group(fdtdconfig.getRoot(),nffft_td_setting_path);

		//read transformer data
		int num_of_td_transformers = TimeDomainNFFFTsettings.getLength();
		for (int tdtrindex=0; tdtrindex<num_of_td_transformers; tdtrindex++)
		{
			Setting& TimeDomainTransformersettings = TimeDomainNFFFTsettings[tdtrindex];	//go to the tdtrindex'th time-domain transformer setting
			//check group for invalid settings
			CheckAngoraGroupSetting(TimeDomainTransformersettings,validsettings);

			if (SettingEnabledForGrid(TimeDomainTransformersettings))		//apply only if enabled for this grid
			{
				TrDataType_td TimeDomainTrData;

				read_value_from_group<double>(TimeDomainTransformersettings,"theta",TimeDomainTrData.THETA);
				TimeDomainTrData.THETA *= M_PI/180;	//convert to radians

				read_value_from_group<double>(TimeDomainTransformersettings,"phi",TimeDomainTrData.PHI);
				TimeDomainTrData.PHI *= M_PI/180;	//convert to radians

				read_optional_value_from_group<bool>(TimeDomainTransformersettings,"write_hertzian_dipole_far_field",TimeDomainTrData.write_hertzian_dipole_far_field);

				double nffft_back_margin_x,nffft_front_margin_x,nffft_left_margin_y,nffft_right_margin_y,nffft_lower_margin_z,nffft_upper_margin_z;
				bool given_in_meters;
//				try{
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_back_margin_x",nffft_back_margin_x,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginBackX = (int)round(nffft_back_margin_x/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginBackX = (int)round(nffft_back_margin_x);
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_front_margin_x",nffft_front_margin_x,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginFrontX = (int)round(nffft_front_margin_x/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginFrontX = (int)round(nffft_front_margin_x);
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_left_margin_y",nffft_left_margin_y,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginLeftY = (int)round(nffft_left_margin_y/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginLeftY = (int)round(nffft_left_margin_y);
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_right_margin_y",nffft_right_margin_y,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginRightY = (int)round(nffft_right_margin_y/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginRightY = (int)round(nffft_right_margin_y);
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_lower_margin_z",nffft_lower_margin_z,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginLowerZ = (int)round(nffft_lower_margin_z/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginLowerZ = (int)round(nffft_lower_margin_z);
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"nffft_upper_margin_z",nffft_upper_margin_z,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTMarginUpperZ = (int)round(nffft_upper_margin_z/dx);
					}
					else
					{
						TimeDomainTrData.NFFFTMarginUpperZ = (int)round(nffft_upper_margin_z);
					}
//				}
//				/** TODO: Better organize the defaults in TrDataType_td **/
//				catch (AngoraSettingNotFoundException&)
//				{//default values
//				}

				double far_field_origin_x,far_field_origin_y,far_field_origin_z;
//				try{
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"far_field_origin_x",far_field_origin_x,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTOriginX = far_field_origin_x/dx + OriginX;  //these are grid cell indices, not coordinates
					}
					else
					{
						TimeDomainTrData.NFFFTOriginX = far_field_origin_x + OriginX;  //these are grid cell indices, not coordinates
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"far_field_origin_y",far_field_origin_y,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTOriginY = far_field_origin_y/dx + OriginY;  //these are grid cell indices, not coordinates
					}
					else
					{
						TimeDomainTrData.NFFFTOriginY = far_field_origin_y + OriginY;  //these are grid cell indices, not coordinates
					}
					if (read_optional_length_from_group<double>(TimeDomainTransformersettings,"far_field_origin_z",far_field_origin_z,given_in_meters))
					if (given_in_meters)
					{
						TimeDomainTrData.NFFFTOriginZ = far_field_origin_z/dx + OriginZ;  //these are grid cell indices, not coordinates
					}
					else
					{
						TimeDomainTrData.NFFFTOriginZ = far_field_origin_z + OriginZ;  //these are grid cell indices, not coordinates
					}
//				}
//				/** TODO: Better organize the defaults in TrDataType_td **/
//				catch (AngoraSettingNotFoundException&)
//				{//default values
//				}

				//read far-field file name
				ostringstream farfieldfilenamestream;
				string FarFieldFilePath,FarFieldFileName;
				//read path
				try{read_value_from_group<string>(TimeDomainTransformersettings,"far_field_dir",FarFieldFilePath);}
				catch (AngoraSettingNotFoundException&)
				{//do nothing if does not exist, but type is checked in the try block
				}
				//add slash to path if necessary
				add_slash_to_path(FarFieldFilePath);
				//if path is not absolute, prepend the base path to get full path
				if (!is_absolute_path(FarFieldFilePath))
				{
					FarFieldFilePath = TimeDomainNFFFTOutputDir + FarFieldFilePath;	//FarFieldFilePath is relative to the recorder output directory
				}
				//create directory if it does not exist
				if (!check_mode)
				{
					if (create_path(FarFieldFilePath)<0)
					{
						/** throw exception **/
						if (rank==0) cout << "Could not create path " << FarFieldFilePath << endl;
					}
				}

				//read file name
				try{read_value_from_group<string>(TimeDomainTransformersettings,"far_field_file_name",FarFieldFileName);}
				catch (AngoraSettingNotFoundException&)
				{
					FarFieldFileName = default_TimeDomainNFFFT_filename;
				}

				//read file extension
				string FarFieldFileExtension;
				if (!read_optional_value_from_group<string>(TimeDomainTransformersettings,"far_field_file_extension",FarFieldFileExtension))
				{
					FarFieldFileExtension = default_td_nffft_fileextension;
				}

				//construct full filename
				string FarFieldFullFileName;
				farfieldfilenamestream << FarFieldFilePath << FarFieldFileName << "_" << GridIndex;
				bool append_group_index_to_file_name;
				if (!read_optional_value_from_group<bool>(TimeDomainTransformersettings,"append_group_index_to_file_name",append_group_index_to_file_name))
				{
					append_group_index_to_file_name = true;
				}
				if (append_group_index_to_file_name)
				{
					farfieldfilenamestream << "_" << tdtrindex;

				}
				if (FarFieldFileExtension!="")
				{
					farfieldfilenamestream << "." << FarFieldFileExtension;
				}

				//get string from string stream
				FarFieldFullFileName = farfieldfilenamestream.str();

				//If not in check mode, add the transformer to the Cnffft_td object
				if (!check_mode)
				{
					TimeDomainTrData.PointSourcesPtr = &PointSources;//pointer to PointSources object is used in theoretical FF calculations
					//finally, add the transformer to the Cnffft_td object
					NFFFT_td.AddTransformer(TimeDomainTrData,FarFieldFullFileName);
				}
			}
		}
	}
}
