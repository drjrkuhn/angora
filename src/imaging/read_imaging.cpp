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

//Includes the routine that reads the optical image definitions into the Cimg object OpticalImage

#include "headers.h"

#include "read_imaging.h"

//definition of Cimgs needed
#include "Cimgs.h"

//definition of ImgDataType needed
#include "Cimg.h"

//For file-directory  manipulations
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern string config_filename;
extern string OutputDir;
//extern string OpticalImagingOutputDir;
extern const string default_OpticalImagingOutputDir;
extern const string default_image_filename;
extern const string default_image_fileextension;

extern double dx;

extern int GridIndex;
extern int rank;

extern void add_slash_to_path(string& path);
extern bool is_absolute_path(string& path);

extern int create_path(const string& path);
extern void MPI_exit(const int& exitcode);


void read_imaging(Cimgs &OpticalImages, const Config& fdtdconfig, const Config& validsettings, const Ctfsf &TFSF)
{
	//Read imaging output directory name
	string OpticalImagingOutputDir;
	if (!fdtdconfig.lookupValue("imaging_output_dir",OpticalImagingOutputDir))
	{
		OpticalImagingOutputDir = default_OpticalImagingOutputDir;
	}
	//if path is not absolute, prepend the base path to get full path
	if (!is_absolute_path(OpticalImagingOutputDir))
	{
		OpticalImagingOutputDir = OutputDir + OpticalImagingOutputDir;
	}
	//add slash to path if necessary
	add_slash_to_path(OpticalImagingOutputDir);
	//create optical imaging output directory if it does not exist
	if (!check_mode)
	{
		if (create_path(OpticalImagingOutputDir)<0)
		{
			/** throw exception **/
		}
	}

	//Read optical imaging settings
	string optical_imaging_setting_path = "OpticalImages";
	if (fdtdconfig.exists(optical_imaging_setting_path))
	{
		Setting& OpticalImagingsettings = read_list_from_group(fdtdconfig.getRoot(),optical_imaging_setting_path);

		//read individual optical image data
		int num_of_optical_images = OpticalImagingsettings.getLength();
		for (int imgindex=0; imgindex<num_of_optical_images; imgindex++)
		{
			Setting& Imagesettings = OpticalImagingsettings[imgindex];	//go to the imgindex'th optical image setting group
			//check group for invalid settings
			CheckAngoraGroupSetting(Imagesettings,validsettings);

			if (SettingEnabledForGrid(Imagesettings))		//apply only if enabled for this grid
			{
				//Add image(s) to the OpticalImages object
				int num_of_lambdas;
				read_value_from_group<int>(Imagesettings,"num_of_lambdas",num_of_lambdas);
				if (num_of_lambdas<=0)
				{
					throw AngoraInvalidSettingValueException(Imagesettings["num_of_lambdas"], "should be positive");
				}
				Array<double,1> lambda(num_of_lambdas);
				double lambda_min,lambda_max;
				if (!read_length_from_group<double>(Imagesettings,"lambda_min",lambda_min))
				{//returned false: given in grid cells
					lambda_min *= dx;
				}
				string lambda_spacing_type;
				if (num_of_lambdas>1)
				{
					if (!read_length_from_group<double>(Imagesettings,"lambda_max",lambda_max))
					{//returned false: given in grid cells
						lambda_max *= dx;
					}
					if (lambda_min>=lambda_max)
					{
						throw AngoraInvalidSettingValueException(Imagesettings["lambda_max"],"should be greater than \"lambda_min\"");
					}
					read_value_from_group<string>(Imagesettings,"lambda_spacing_type",lambda_spacing_type);
					if ((lambda_spacing_type!="log")&&(lambda_spacing_type!="lambda-linear")&&(lambda_spacing_type!="k-linear"))
					{
						throw AngoraInvalidSettingValueException(Imagesettings["lambda_spacing_type"],"should be either \"log\", \"lambda-linear\" or \"k-linear\"");
					}
				}
				bool do_not_include_first_lambda;
				try{read_value_from_group<bool>(Imagesettings,"do_not_include_first_lambda",do_not_include_first_lambda);}
				catch (AngoraSettingNotFoundException&)
				{
					do_not_include_first_lambda = false;
				}
				bool do_not_include_last_lambda;
				try{read_value_from_group<bool>(Imagesettings,"do_not_include_last_lambda",do_not_include_last_lambda);}
				catch (AngoraSettingNotFoundException&)
				{
					do_not_include_last_lambda = false;
				}
				if (lambda.size()==1)
				{
					lambda = lambda_min;
				}
				else if (lambda_spacing_type=="lambda-linear")
				{//place lambda values linearly between lambda_min and lambda_max
					double lambda_start,lambda_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min + lambda_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min + lambda_step;
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							lambda_step = (lambda_max-lambda_min)/lambda.size();
							lambda_start = lambda_min;
						}
						else
						{//BOTH the first and the last lambda are included
							lambda_step = (lambda_max-lambda_min)/(lambda.size()-1);
							lambda_start = lambda_min;
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//linearly-spaced lambda values
						lambda(i) = lambda_start + i*lambda_step;
					}
				}
				else if (lambda_spacing_type=="k-linear")
				{//place k values linearly between 2pi/lambda_max and 2pi/lambda_min
					double k,k_start,k_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max + k_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max;
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/lambda.size();
							k_start = 2*M_PI/lambda_max + k_step;
						}
						else
						{//BOTH the first and the last lambda are included
							k_step = (2*M_PI/lambda_min-2*M_PI/lambda_max)/(lambda.size()-1);
							k_start = 2*M_PI/lambda_max;
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//linearly-spaced k values
						k = k_start + i*k_step;
						//convert to lambda
						lambda(i) = 2*M_PI/k;
					}
				}
				else if (lambda_spacing_type=="log")
				{//place k values logarithmically between 2pi/lambda_max and 2pi/lambda_min
					double klog,klog_start,klog_step;
					if (do_not_include_first_lambda)
					{
						if (do_not_include_last_lambda)
						{//NEITHER the first nor the last lambda is included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max) + klog_step/2;
						}
						else
						{//first lambda is NOT included, last lambda is included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max);
						}
					}
					else
					{
						if (do_not_include_last_lambda)
						{//first lambda is included, last lambda is NOT included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/lambda.size();
							klog_start = log10(2*M_PI/lambda_max) + klog_step;
						}
						else
						{//BOTH the first and the last lambda are included
							klog_step = (log10(2*M_PI/lambda_min)-log10(2*M_PI/lambda_max))/(lambda.size()-1);
							klog_start = log10(2*M_PI/lambda_max);
						}
					}
					for (int i=0;i<lambda.size();i++)
					{
						//logarithmically-spaced k values
						klog = klog_start + i*klog_step;
						//convert to lambda
						lambda(i) = 2*M_PI/pow(10.0,klog);
					}
				}

				//the lambda array is built; now define the optical image data object
				ImgDataType ImgData(lambda);

				//read the image arrays to be computed and written into file
				string output_data_optionname = "output_data";
				Setting& output_data_array = read_array_from_group(Imagesettings,output_data_optionname);
				if (output_data_array.getLength()==0)
				{
					throw AngoraSettingEmptyException(output_data_array);
				}
				if (output_data_array[0].getType()!=Setting::TypeString)
				{
					throw AngoraInvalidSettingTypeException(output_data_array,"should be a string array");
				}
				for (int listindex=0; listindex<output_data_array.getLength(); listindex++)
				{
					string current_list_member = (const char*)output_data_array[listindex];
					ImgData.add_output_data_item_to_list(current_list_member);
				}

				//read the aperture half angle
				read_value_from_group<double>(Imagesettings,"ap_half_angle",ImgData.ap_half_angle);
				ImgData.ap_half_angle *= M_PI/180;	//convert to radians

				try{read_value_from_group<double>(Imagesettings,"magnification",ImgData.magnification);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.magnification = 1;
				}

				try{read_value_from_group<string>(Imagesettings,"coll_half_space",ImgData.coll_half_space);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.coll_half_space = "upper";
				}

				if ((ImgData.coll_half_space!="upper")&&(ImgData.coll_half_space!="lower"))
				{
					throw AngoraInvalidSettingValueException(Imagesettings["coll_half_space"],"should be either \"upper\" or \"lower\"");
				}

				try{read_value_from_group<double>(Imagesettings,"image_space_refr_index",ImgData.img_space_refr_index);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.img_space_refr_index = 1;
				}
				try{read_value_from_group<double>(Imagesettings,"image_expansion_factor_x",ImgData.img_expansionfactor_x);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.img_expansionfactor_x = 1;
				}
				try{read_value_from_group<double>(Imagesettings,"image_expansion_factor_y",ImgData.img_expansionfactor_y);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.img_expansionfactor_y = 1;
				}
				try{read_value_from_group<double>(Imagesettings,"image_oversampling_rate_x",ImgData.img_oversamplingrate_x);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.img_oversamplingrate_x = 1;
				}
				try{read_value_from_group<double>(Imagesettings,"image_oversampling_rate_y",ImgData.img_oversamplingrate_y);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgData.img_oversamplingrate_y = 1;
				}

				bool choose_smallest_power_of_2;
				read_optional_value_from_group<bool>(Imagesettings,"choose_smallest_power_of_2",ImgData.choose_smallest_power_of_2);

				//read NFFFT box margins
				double nffft_back_margin_x,nffft_front_margin_x,nffft_left_margin_y,nffft_right_margin_y,nffft_lower_margin_z,nffft_upper_margin_z;
				bool given_in_meters;
//				try{
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_back_margin_x",nffft_back_margin_x,given_in_meters))
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginBackX = (int)round(nffft_back_margin_x/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginBackX = (int)round(nffft_back_margin_x);
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_front_margin_x",nffft_front_margin_x,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginFrontX = (int)round(nffft_front_margin_x/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginFrontX = (int)round(nffft_front_margin_x);
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_left_margin_y",nffft_left_margin_y,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginLeftY = (int)round(nffft_left_margin_y/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginLeftY = (int)round(nffft_left_margin_y);
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_right_margin_y",nffft_right_margin_y,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginRightY = (int)round(nffft_right_margin_y/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginRightY = (int)round(nffft_right_margin_y);
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_lower_margin_z",nffft_lower_margin_z,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginLowerZ = (int)round(nffft_lower_margin_z/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginLowerZ = (int)round(nffft_lower_margin_z);
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"nffft_upper_margin_z",nffft_upper_margin_z,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.NFFFTMarginUpperZ = (int)round(nffft_upper_margin_z/dx);
                                            }
                                            else
                                            {
                                                    ImgData.NFFFTMarginUpperZ = (int)round(nffft_upper_margin_z);
                                            }
                                        }
//				}
//				/** TODO: Better organize the defaults in ImgDataType **/
//				catch (AngoraSettingNotFoundException&)
//				{//default values
//				}

				//read image origin
				double image_origin_x,image_origin_y,image_origin_z;
//				try{
					if (read_optional_length_from_group<double>(Imagesettings,"image_origin_x",image_origin_x,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.ImgOriginX = image_origin_x/dx + OriginX;  //these are grid cell indices, not coordinates
                                            }
                                            else
                                            {
                                                    ImgData.ImgOriginX = image_origin_x + OriginX;  //these are grid cell indices, not coordinates
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"image_origin_y",image_origin_y,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.ImgOriginY = image_origin_y/dx + OriginY;  //these are grid cell indices, not coordinates
                                            }
                                            else
                                            {
                                                    ImgData.ImgOriginY = image_origin_y + OriginY;  //these are grid cell indices, not coordinates
                                            }
                                        }
					if (read_optional_length_from_group<double>(Imagesettings,"image_origin_z",image_origin_z,given_in_meters)) 
                                        {
                                            if (given_in_meters)
                                            {
                                                    ImgData.ImgOriginZ = image_origin_z/dx + OriginZ;  //these are grid cell indices, not coordinates
                                            }
                                            else
                                            {
                                                    ImgData.ImgOriginZ = image_origin_z + OriginZ;  //these are grid cell indices, not coordinates
                                            }
                                        }
//				}
//				/** TODO: Better organize the defaults in TrDataType_pd **/
//				catch (AngoraSettingNotFoundException&)
//				{//default values
//				}

				//read far-field file name
				ostringstream Imgfilenamestream;
				string ImgFilePath,ImgFileName;
				//read path
				try{read_value_from_group<string>(Imagesettings,"image_dir",ImgFilePath);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgFilePath = ".";
				}
				//add slash to path if necessary
				add_slash_to_path(ImgFilePath);
				//add path to string stream
				Imgfilenamestream << OpticalImagingOutputDir << ImgFilePath;
				//get string from string stream
				string ImgFullFilePath = Imgfilenamestream.str();
				//create transformer output directory if it does not exist
				if (!check_mode)
				{
					if (create_path(ImgFullFilePath)<0)
					{
						/** throw exception **/
					}
				}

				//read file name
				try{read_value_from_group<string>(Imagesettings,"image_file_name",ImgFileName);}
				catch (AngoraSettingNotFoundException&)
				{
					ImgFileName = default_image_filename;
				}

				//read file extension
				string ImgFileExtension;
				if (!read_optional_value_from_group<string>(Imagesettings,"image_file_extension",ImgFileExtension))
				{
					ImgFileExtension = default_image_fileextension;
				}

				//construct full filename
				string ImgFullFileName;
				Imgfilenamestream << ImgFileName << "_" << GridIndex;
				bool append_group_index_to_file_name;
				if (!read_optional_value_from_group<bool>(Imagesettings,"append_group_index_to_file_name",append_group_index_to_file_name))
				{
					append_group_index_to_file_name = true;
				}
				if (append_group_index_to_file_name)
				{
					Imgfilenamestream << "_" << imgindex;

				}
				if (ImgFileExtension!="")
				{
					Imgfilenamestream << "." << ImgFileExtension;
				}

				//get string from string stream
				ImgFullFileName = Imgfilenamestream.str();

				//If not in check mode, add the transformer to the Cnffft_pd object
				if (!check_mode)
				{
					ImgData.TFSFPtr = &TFSF;
					//add the optical image to the Cimgs object
					OpticalImages.AddOpticalImage(ImgData,ImgFullFileName);
				}
			}
		}
	}
}
