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

//Definition of the Cmat object that represents a material

#include "Cmat_excp.h"

#include "Cmat.h"

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,tau_r_x,tau_r_y,tau_r_z;


/** MOVE INTO SEPARATE HEADER **/
eps_x_type add_distinct_eps_x(const float& epsilon_r_x)
{
	eps_x_type num_of_distinct_eps_x = eps_x.size()+1;
	eps_x.resizeAndPreserve(num_of_distinct_eps_x);
	eps_x(num_of_distinct_eps_x-1) = epsilon_r_x;
	return num_of_distinct_eps_x-1;
}

eps_y_type add_distinct_eps_y(const float& epsilon_r_y)
{
	eps_y_type num_of_distinct_eps_y = eps_y.size()+1;
	eps_y.resizeAndPreserve(num_of_distinct_eps_y);
	eps_y(num_of_distinct_eps_y-1) = epsilon_r_y;
	return num_of_distinct_eps_y-1;
}

eps_z_type add_distinct_eps_z(const float& epsilon_r_z)
{
	eps_z_type num_of_distinct_eps_z = eps_z.size()+1;
	eps_z.resizeAndPreserve(num_of_distinct_eps_z);
	eps_z(num_of_distinct_eps_z-1) = epsilon_r_z;
	return num_of_distinct_eps_z-1;
}

mu_x_type add_distinct_mu_x(const float& mu_r_x)
{
	mu_x_type num_of_distinct_mu_x = mu_x.size()+1;
	mu_x.resizeAndPreserve(num_of_distinct_mu_x);
	mu_x(num_of_distinct_mu_x-1) = mu_r_x;
	return num_of_distinct_mu_x-1;
}

mu_y_type add_distinct_mu_y(const float& mu_r_y)
{
	mu_y_type num_of_distinct_mu_y = mu_y.size()+1;
	mu_y.resizeAndPreserve(num_of_distinct_mu_y);
	mu_y(num_of_distinct_mu_y-1) = mu_r_y;
	return num_of_distinct_mu_y-1;
}

mu_z_type add_distinct_mu_z(const float& mu_r_z)
{
	mu_z_type num_of_distinct_mu_z = mu_z.size()+1;
	mu_z.resizeAndPreserve(num_of_distinct_mu_z);
	mu_z(num_of_distinct_mu_z-1) = mu_r_z;
	return num_of_distinct_mu_z-1;
}

cond_e_x_type add_distinct_cond_e_x(const float& sigma_e_x)
{
	cond_e_x_type num_of_distinct_cond_e_x = cond_e_x.size()+1;
	cond_e_x.resizeAndPreserve(num_of_distinct_cond_e_x);
	cond_e_x(num_of_distinct_cond_e_x-1) = sigma_e_x;
	return num_of_distinct_cond_e_x-1;
}

cond_e_y_type add_distinct_cond_e_y(const float& sigma_e_y)
{
	cond_e_y_type num_of_distinct_cond_e_y = cond_e_y.size()+1;
	cond_e_y.resizeAndPreserve(num_of_distinct_cond_e_y);
	cond_e_y(num_of_distinct_cond_e_y-1) = sigma_e_y;
	return num_of_distinct_cond_e_y-1;
}

cond_e_z_type add_distinct_cond_e_z(const float& sigma_e_z)
{
	cond_e_z_type num_of_distinct_cond_e_z = cond_e_z.size()+1;
	cond_e_z.resizeAndPreserve(num_of_distinct_cond_e_z);
	cond_e_z(num_of_distinct_cond_e_z-1) = sigma_e_z;
	return num_of_distinct_cond_e_z-1;
}

cond_h_x_type add_distinct_cond_h_x(const float& sigma_h_x)
{
	cond_h_x_type num_of_distinct_cond_h_x = cond_h_x.size()+1;
	cond_h_x.resizeAndPreserve(num_of_distinct_cond_h_x);
	cond_h_x(num_of_distinct_cond_h_x-1) = sigma_h_x;
	return num_of_distinct_cond_h_x-1;
}

cond_h_y_type add_distinct_cond_h_y(const float& sigma_h_y)
{
	cond_h_y_type num_of_distinct_cond_h_y = cond_h_y.size()+1;
	cond_h_y.resizeAndPreserve(num_of_distinct_cond_h_y);
	cond_h_y(num_of_distinct_cond_h_y-1) = sigma_h_y;
	return num_of_distinct_cond_h_y-1;
}

cond_h_z_type add_distinct_cond_h_z(const float& sigma_h_z)
{
	cond_h_z_type num_of_distinct_cond_h_z = cond_h_z.size()+1;
	cond_h_z.resizeAndPreserve(num_of_distinct_cond_h_z);
	cond_h_z(num_of_distinct_cond_h_z-1) = sigma_h_z;
	return num_of_distinct_cond_h_z-1;
}

omega_p_x_type add_distinct_omega_p_x(const float& omega)
{
	omega_p_x_type num_of_distinct_omega_p_x = omega_p_x.size()+1;
	omega_p_x.resizeAndPreserve(num_of_distinct_omega_p_x);
	omega_p_x(num_of_distinct_omega_p_x-1) = omega;
	return num_of_distinct_omega_p_x-1;
}

omega_p_y_type add_distinct_omega_p_y(const float& omega)
{
	omega_p_y_type num_of_distinct_omega_p_y = omega_p_y.size()+1;
	omega_p_y.resizeAndPreserve(num_of_distinct_omega_p_y);
	omega_p_y(num_of_distinct_omega_p_y-1) = omega;
	return num_of_distinct_omega_p_y-1;
}

omega_p_z_type add_distinct_omega_p_z(const float& omega)
{
	omega_p_z_type num_of_distinct_omega_p_z = omega_p_z.size()+1;
	omega_p_z.resizeAndPreserve(num_of_distinct_omega_p_z);
	omega_p_z(num_of_distinct_omega_p_z-1) = omega;
	return num_of_distinct_omega_p_z-1;
}

tau_r_x_type add_distinct_tau_r_x(const float& tau)
{
	tau_r_x_type num_of_distinct_tau_r_x = tau_r_x.size()+1;
	tau_r_x.resizeAndPreserve(num_of_distinct_tau_r_x);
	tau_r_x(num_of_distinct_tau_r_x-1) = tau;
	return num_of_distinct_tau_r_x-1;
}

tau_r_y_type add_distinct_tau_r_y(const float& tau)
{
	tau_r_y_type num_of_distinct_tau_r_y = tau_r_y.size()+1;
	tau_r_y.resizeAndPreserve(num_of_distinct_tau_r_y);
	tau_r_y(num_of_distinct_tau_r_y-1) = tau;
	return num_of_distinct_tau_r_y-1;
}

tau_r_z_type add_distinct_tau_r_z(const float& tau)
{
	tau_r_z_type num_of_distinct_tau_r_z = tau_r_z.size()+1;
	tau_r_z.resizeAndPreserve(num_of_distinct_tau_r_z);
	tau_r_z(num_of_distinct_tau_r_z-1) = tau;
	return num_of_distinct_tau_r_z-1;
}
/** MOVE INTO SEPARATE HEADER **/


void Cmat::set_eps(const float& epsilon_r)
{
	set_eps_x(epsilon_r);
	set_eps_y(epsilon_r);
	set_eps_z(epsilon_r);
}
void Cmat::set_eps_x(const float& epsilon_r)
{
	_eps_x_ptr.reset(new eps_x_type(add_distinct_eps_x(epsilon_r)));
}
void Cmat::set_eps_y(const float& epsilon_r)
{
	_eps_y_ptr.reset(new eps_y_type(add_distinct_eps_y(epsilon_r)));
}
void Cmat::set_eps_z(const float& epsilon_r)
{
	_eps_z_ptr.reset(new eps_z_type(add_distinct_eps_z(epsilon_r)));
}

void Cmat::set_mu(const float& mu_r)
{
	set_mu_x(mu_r);
	set_mu_y(mu_r);
	set_mu_z(mu_r);
}
void Cmat::set_mu_x(const float& mu_r)
{
	_mu_x_ptr.reset(new mu_x_type(add_distinct_mu_x(mu_r)));
}
void Cmat::set_mu_y(const float& mu_r)
{
	_mu_y_ptr.reset(new mu_y_type(add_distinct_mu_y(mu_r)));
}
void Cmat::set_mu_z(const float& mu_r)
{
	_mu_z_ptr.reset(new mu_z_type(add_distinct_mu_z(mu_r)));
}

void Cmat::set_cond_e(const float& sigma_e)
{
	set_cond_e_x(sigma_e);
	set_cond_e_y(sigma_e);
	set_cond_e_z(sigma_e);
}
void Cmat::set_cond_e_x(const float& sigma_e)
{
	_cond_e_x_ptr.reset(new cond_e_x_type(add_distinct_cond_e_x(sigma_e)));
}
void Cmat::set_cond_e_y(const float& sigma_e)
{
	_cond_e_y_ptr.reset(new cond_e_y_type(add_distinct_cond_e_y(sigma_e)));
}
void Cmat::set_cond_e_z(const float& sigma_e)
{
	_cond_e_z_ptr.reset(new cond_e_z_type(add_distinct_cond_e_z(sigma_e)));
}

void Cmat::set_cond_h(const float& sigma_h)
{
	set_cond_h_x(sigma_h);
	set_cond_h_y(sigma_h);
	set_cond_h_z(sigma_h);
}
void Cmat::set_cond_h_x(const float& sigma_h)
{
	_cond_h_x_ptr.reset(new cond_h_x_type(add_distinct_cond_h_x(sigma_h)));
}
void Cmat::set_cond_h_y(const float& sigma_h)
{
	_cond_h_y_ptr.reset(new cond_h_y_type(add_distinct_cond_h_y(sigma_h)));
}
void Cmat::set_cond_h_z(const float& sigma_h)
{
	_cond_h_z_ptr.reset(new cond_h_z_type(add_distinct_cond_h_z(sigma_h)));
}

void Cmat::set_omega_p(const float& omega)
{
	set_omega_p_x(omega);
	set_omega_p_y(omega);
	set_omega_p_z(omega);
}

void Cmat::set_omega_p_x(const float& omega)
{
	_omega_p_x_ptr.reset(new omega_p_x_type(add_distinct_omega_p_x(omega)));
}

void Cmat::set_omega_p_y(const float& omega)
{
	_omega_p_y_ptr.reset(new omega_p_y_type(add_distinct_omega_p_y(omega)));
}

void Cmat::set_omega_p_z(const float& omega)
{
	_omega_p_z_ptr.reset(new omega_p_z_type(add_distinct_omega_p_z(omega)));
}

void Cmat::set_tau_r(const float& tau)
{
	set_tau_r_x(tau);
	set_tau_r_y(tau);
	set_tau_r_z(tau);
}

void Cmat::set_tau_r_x(const float& tau)
{
	_tau_r_x_ptr.reset(new tau_r_x_type(add_distinct_tau_r_x(tau)));
}

void Cmat::set_tau_r_y(const float& tau)
{
	_tau_r_y_ptr.reset(new tau_r_y_type(add_distinct_tau_r_y(tau)));
}

void Cmat::set_tau_r_z(const float& tau)
{
	_tau_r_z_ptr.reset(new tau_r_z_type(add_distinct_tau_r_z(tau)));
}

bool Cmat::eps_x_exists() const
{
	return _eps_x_ptr;
}

bool Cmat::eps_y_exists() const
{
	return _eps_y_ptr;
}

bool Cmat::eps_z_exists() const
{
	return _eps_z_ptr;
}

bool Cmat::mu_x_exists() const
{
	return _mu_x_ptr;
}

bool Cmat::mu_y_exists() const
{
	return _mu_y_ptr;
}

bool Cmat::mu_z_exists() const
{
	return _mu_z_ptr;
}

bool Cmat::cond_e_x_exists() const
{
	return _cond_e_x_ptr;
}

bool Cmat::cond_e_y_exists() const
{
	return _cond_e_y_ptr;
}

bool Cmat::cond_e_z_exists() const
{
	return _cond_e_z_ptr;
}

bool Cmat::cond_h_x_exists() const
{
	return _cond_h_x_ptr;
}

bool Cmat::cond_h_y_exists() const
{
	return _cond_h_y_ptr;
}

bool Cmat::cond_h_z_exists() const
{
	return _cond_h_z_ptr;
}

bool Cmat::omega_p_x_exists() const
{
	return _omega_p_x_ptr;
}

bool Cmat::omega_p_y_exists() const
{
	return _omega_p_y_ptr;
}

bool Cmat::omega_p_z_exists() const
{
	return _omega_p_z_ptr;
}

bool Cmat::tau_r_x_exists() const
{
	return _tau_r_x_ptr;
}

bool Cmat::tau_r_y_exists() const
{
	return _tau_r_y_ptr;
}

bool Cmat::tau_r_z_exists() const
{
	return _tau_r_z_ptr;
}

eps_x_type Cmat::eps_x_index() const
{
	if (eps_x_exists())
	{
		return *_eps_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_x");
	}
}

eps_y_type Cmat::eps_y_index() const
{
	if (eps_y_exists())
	{
		return *_eps_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_y");
	}
}

eps_z_type Cmat::eps_z_index() const
{
	if (eps_z_exists())
	{
		return *_eps_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_z");
	}
}

mu_x_type Cmat::mu_x_index() const
{
	if (mu_x_exists())
	{
		return *_mu_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_x");
	}
}

mu_y_type Cmat::mu_y_index() const
{
	if (mu_y_exists())
	{
		return *_mu_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_y");
	}
}

mu_z_type Cmat::mu_z_index() const
{
	if (mu_z_exists())
	{
		return *_mu_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_z");
	}
}

cond_e_x_type Cmat::cond_e_x_index() const
{
	if (cond_e_x_exists())
	{
		return *_cond_e_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_x");
	}
}

cond_e_y_type Cmat::cond_e_y_index() const
{
	if (cond_e_y_exists())
	{
		return *_cond_e_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_y");
	}
}

cond_e_z_type Cmat::cond_e_z_index() const
{
	if (cond_e_z_exists())
	{
		return *_cond_e_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_z");
	}
}

cond_h_x_type Cmat::cond_h_x_index() const
{
	if (cond_h_x_exists())
	{
		return *_cond_h_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_x");
	}
}

cond_h_y_type Cmat::cond_h_y_index() const
{
	if (cond_h_y_exists())
	{
		return *_cond_h_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_y");
	}
}

cond_h_z_type Cmat::cond_h_z_index() const
{
	if (cond_h_z_exists())
	{
		return *_cond_h_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_z");
	}
}

omega_p_x_type Cmat::omega_p_x_index() const
{
	if (omega_p_x_exists())
	{
		return *_omega_p_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_x");
	}
}

omega_p_y_type Cmat::omega_p_y_index() const
{
	if (omega_p_y_exists())
	{
		return *_omega_p_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_y");
	}
}

omega_p_z_type Cmat::omega_p_z_index() const
{
	if (omega_p_z_exists())
	{
		return *_omega_p_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_z");
	}
}

tau_r_x_type Cmat::tau_r_x_index() const
{
	if (tau_r_x_exists())
	{
		return *_tau_r_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_x");
	}
}

tau_r_y_type Cmat::tau_r_y_index() const
{
	if (tau_r_y_exists())
	{
		return *_tau_r_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_y");
	}
}

tau_r_z_type Cmat::tau_r_z_index() const
{
	if (tau_r_z_exists())
	{
		return *_tau_r_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_z");
	}
}

float Cmat::eps_x_value() const
{
	if (eps_x_exists())
	{
		return eps_x(*_eps_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_x");
	}
}

float Cmat::eps_y_value() const
{
	if (eps_y_exists())
	{
		return eps_y(*_eps_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_y");
	}
}

float Cmat::eps_z_value() const
{
	if (eps_z_exists())
	{
		return eps_z(*_eps_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_z");
	}
}

float Cmat::mu_x_value() const
{
	if (mu_x_exists())
	{
		return mu_x(*_mu_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_x");
	}
}

float Cmat::mu_y_value() const
{
	if (mu_y_exists())
	{
		return mu_y(*_mu_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_y");
	}
}

float Cmat::mu_z_value() const
{
	if (mu_z_exists())
	{
		return mu_z(*_mu_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_z");
	}
}

float Cmat::cond_e_x_value() const
{
	if (cond_e_x_exists())
	{
		return cond_e_x(*_cond_e_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_x");
	}
}

float Cmat::cond_e_y_value() const
{
	if (cond_e_y_exists())
	{
		return cond_e_y(*_cond_e_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_y");
	}
}

float Cmat::cond_e_z_value() const
{
	if (cond_e_z_exists())
	{
		return cond_e_z(*_cond_e_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_z");
	}
}

float Cmat::cond_h_x_value() const
{
	if (cond_h_x_exists())
	{
		return cond_h_x(*_cond_h_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_x");
	}
}

float Cmat::cond_h_y_value() const
{
	if (cond_h_y_exists())
	{
		return cond_h_y(*_cond_h_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_y");
	}
}

float Cmat::cond_h_z_value() const
{
	if (cond_h_z_exists())
	{
		return cond_h_z(*_cond_h_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_z");
	}
}

float Cmat::omega_p_x_value() const
{
	if (omega_p_x_exists())
	{
		return omega_p_x(*_omega_p_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_x");
	}
}

float Cmat::omega_p_y_value() const
{
	if (omega_p_y_exists())
	{
		return omega_p_y(*_omega_p_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_y");
	}
}

float Cmat::omega_p_z_value() const
{
	if (omega_p_z_exists())
	{
		return omega_p_z(*_omega_p_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("omega_p_z");
	}
}

float Cmat::tau_r_x_value() const
{
	if (tau_r_x_exists())
	{
		return tau_r_x(*_tau_r_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_x");
	}
}

float Cmat::tau_r_y_value() const
{
	if (tau_r_y_exists())
	{
		return tau_r_y(*_tau_r_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_y");
	}
}

float Cmat::tau_r_z_value() const
{
	if (tau_r_z_exists())
	{
		return tau_r_z(*_tau_r_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("tau_r_z");
	}
}
