/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_craig_visc_h
#define _aspect_material_model_craig_visc_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     *Craig fucking around with mantle viscosity mostly
     * */
    template <int dim>
    class MantleHack : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Evaluate material properties.
         */
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::vector<double> densities;

        std::vector<double> thermal_expansivities;
        double thermal_diffusivity;
        double gamma; // Coefficient of yield stress increase with depth
        double heat_capacity;
        double ref_strain_rate;
        double tau_0; // cohesive strength of rocks at the surface
        double reference_T;
        double min_strain_rate;
        double visc_min;
        double visc_max;
        double ref_visc;
        std::vector<double> UM_activation_energies;
        std::vector<double> LM_activation_energies;
        std::vector<double> UM_B;  // Preexponential constant in the viscous rheology law B
        std::vector<double> nvs; // Stress exponent, viscous rheology
        std::vector<double> nps;//Stress exponent, plastic rheology
        double UM_activation_volume;
        double LM_activation_volume;
        double Brig_sol_a;
        double Brig_sol_b;
        double Brig_sol_c;
        double Brig_sol_d;
        double Brig_sol_e;
        double Brig_sol_f;

        double Ferr_sol_a;
        double Ferr_sol_b;
        double Ferr_sol_c;
        double Ferr_sol_d;
        double Ferr_sol_e;
        double Ferr_sol_f;
        double Ferr_sol_g;

	double UM_sol_a;
        double UM_sol_b;
        double UM_sol_c;
        double UM_sol_d;

        double Brig_D0;
        double Brig_g0;
        double Brig_gr;
        double Brig_om;
        double Ferr_D0;
        double Ferr_g0;
        double Ferr_gr;
        double Ferr_om;
        double Ji;
	double Bdg_Dens_Multiplier;
	double exps_decrease;

	MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;
	//double UM_B; // Preexponential constant in the viscous rheology law B
        double LM_B; // Preexponential constant in the viscous rheology law B
        double yield_max;

    };
  }
}

#endif
