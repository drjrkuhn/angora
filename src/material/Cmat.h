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

//Declaration of the Cmat object that represents a material

#ifndef CMAT_H
#define CMAT_H

#include "Cmat_excp.h"

#include "Cmat_types.h"

//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

typedef boost::shared_ptr<eps_x_type> eps_x_type_ptr;
typedef boost::shared_ptr<eps_y_type> eps_y_type_ptr;
typedef boost::shared_ptr<eps_z_type> eps_z_type_ptr;
typedef boost::shared_ptr<mu_x_type> mu_x_type_ptr;
typedef boost::shared_ptr<mu_y_type> mu_y_type_ptr;
typedef boost::shared_ptr<mu_z_type> mu_z_type_ptr;
typedef boost::shared_ptr<cond_e_x_type> cond_e_x_type_ptr;
typedef boost::shared_ptr<cond_e_y_type> cond_e_y_type_ptr;
typedef boost::shared_ptr<cond_e_z_type> cond_e_z_type_ptr;
typedef boost::shared_ptr<cond_h_x_type> cond_h_x_type_ptr;
typedef boost::shared_ptr<cond_h_y_type> cond_h_y_type_ptr;
typedef boost::shared_ptr<cond_h_z_type> cond_h_z_type_ptr;
typedef boost::shared_ptr<omega_p_x_type> omega_p_x_type_ptr;
typedef boost::shared_ptr<omega_p_y_type> omega_p_y_type_ptr;
typedef boost::shared_ptr<omega_p_z_type> omega_p_z_type_ptr;
typedef boost::shared_ptr<tau_r_x_type> tau_r_x_type_ptr;
typedef boost::shared_ptr<tau_r_y_type> tau_r_y_type_ptr;
typedef boost::shared_ptr<tau_r_z_type> tau_r_z_type_ptr;


class Cmat
{
public:
	void set_to_vacuum();
	void set_eps(const float& epsilon_r);
	void set_eps_x(const float& epsilon_r_x);
	void set_eps_y(const float& epsilon_r_y);
	void set_eps_z(const float& epsilon_r_z);
	void set_mu(const float& mu_r);
	void set_mu_x(const float& mu_r_x);
	void set_mu_y(const float& mu_r_y);
	void set_mu_z(const float& mu_r_z);
	void set_cond_e(const float& cond_e);
	void set_cond_e_x(const float& cond_e_x);
	void set_cond_e_y(const float& cond_e_y);
	void set_cond_e_z(const float& cond_e_z);
	void set_cond_h(const float& cond_h);
	void set_cond_h_x(const float& cond_h_x);
	void set_cond_h_y(const float& cond_h_y);
	void set_cond_h_z(const float& cond_h_z);
	void set_omega_p(const float& omega);
	void set_omega_p_x(const float& omega);
	void set_omega_p_y(const float& omega);
	void set_omega_p_z(const float& omega);
	void set_tau_r(const float& tau);
	void set_tau_r_x(const float& tau);
	void set_tau_r_y(const float& tau);
	void set_tau_r_z(const float& tau);

	bool eps_x_exists() const;
	bool eps_y_exists() const;
	bool eps_z_exists() const;
	bool mu_x_exists() const;
	bool mu_y_exists() const;
	bool mu_z_exists() const;
	bool cond_e_x_exists() const;
	bool cond_e_y_exists() const;
	bool cond_e_z_exists() const;
	bool cond_h_x_exists() const;
	bool cond_h_y_exists() const;
	bool cond_h_z_exists() const;
	bool omega_p_x_exists() const;
	bool omega_p_y_exists() const;
	bool omega_p_z_exists() const;
	bool tau_r_x_exists() const;
	bool tau_r_y_exists() const;
	bool tau_r_z_exists() const;

	eps_x_type eps_x_index() const;
	eps_y_type eps_y_index() const;
	eps_z_type eps_z_index() const;
	mu_x_type mu_x_index() const;
	mu_y_type mu_y_index() const;
	mu_z_type mu_z_index() const;
	cond_e_x_type cond_e_x_index() const;
	cond_e_y_type cond_e_y_index() const;
	cond_e_z_type cond_e_z_index() const;
	cond_h_x_type cond_h_x_index() const;
	cond_h_y_type cond_h_y_index() const;
	cond_h_z_type cond_h_z_index() const;
	omega_p_x_type omega_p_x_index() const;
	omega_p_y_type omega_p_y_index() const;
	omega_p_z_type omega_p_z_index() const;
	tau_r_x_type tau_r_x_index() const;
	tau_r_y_type tau_r_y_index() const;
	tau_r_z_type tau_r_z_index() const;

	float eps_x_value() const;
	float eps_y_value() const;
	float eps_z_value() const;
	float mu_x_value() const;
	float mu_y_value() const;
	float mu_z_value() const;
	float cond_e_x_value() const;
	float cond_e_y_value() const;
	float cond_e_z_value() const;
	float cond_h_x_value() const;
	float cond_h_y_value() const;
	float cond_h_z_value() const;
	float omega_p_x_value() const;
	float omega_p_y_value() const;
	float omega_p_z_value() const;
	float tau_r_x_value() const;
	float tau_r_y_value() const;
	float tau_r_z_value() const;

private:
	eps_x_type_ptr _eps_x_ptr;
	eps_y_type_ptr _eps_y_ptr;
	eps_z_type_ptr _eps_z_ptr;
	mu_x_type_ptr _mu_x_ptr;
	mu_y_type_ptr _mu_y_ptr;
	mu_z_type_ptr _mu_z_ptr;
	cond_e_x_type_ptr _cond_e_x_ptr;
	cond_e_y_type_ptr _cond_e_y_ptr;
	cond_e_z_type_ptr _cond_e_z_ptr;
	cond_h_x_type_ptr _cond_h_x_ptr;
	cond_h_y_type_ptr _cond_h_y_ptr;
	cond_h_z_type_ptr _cond_h_z_ptr;
	omega_p_x_type_ptr _omega_p_x_ptr;
	omega_p_y_type_ptr _omega_p_y_ptr;
	omega_p_z_type_ptr _omega_p_z_ptr;
	tau_r_x_type_ptr _tau_r_x_ptr;
	tau_r_y_type_ptr _tau_r_y_ptr;
	tau_r_z_type_ptr _tau_r_z_ptr;

};


#endif // CMAT_H
