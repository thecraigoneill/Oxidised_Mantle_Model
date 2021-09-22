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

#include <string>
#include "mantle_bridgmanite.h"
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/visco_plastic.h>


namespace aspect
{
  template <int dim>
  namespace MaterialModel
  {
    void
    MantleHack<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      const double R = 8.3145; // J mol-1 K-1
      /*double vmax = 1e20;
      double vmin = 1e20;
      double tauy_max = 0.0;
      double tauy_min = 0.0;
      double edot_min = 0.0;
      double edot_max = 0.0; */
      const double g = 9.81; // m s-2
      //printf("###########################################################################################\n");
      int ii=0;
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition);
          //printf("T = %.3e , %.3e %.3e \n",temperature, in.temperature[i],composition);
          SymmetricTensor<2,dim> strain_rate;
          if (in.strain_rate.size()) {
            strain_rate = in.strain_rate[i];
            //printf("%d %d Strain rate size %.2e ",ii,i,strain_rate[0][0]);  // Does the 8 values correspond to an FEM element? DoFs? Quadrature points?
          }
          // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
          // The first time this function is called (first iteration of first time step)
          // a specified "reference" strain rate is used as the returned value would
          // otherwise be zero.
          const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                               ?
                               ref_strain_rate
                               :
                               std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate) );

          const double depth = this->get_geometry_model().depth(position); // units: m
          { //
            // Calculate mantle viscosity
            double veffv = 0.0;
            double veffp = 0.0;
            /*double UM_activation_energy = 0.0;
            double LM_activation_energy = 0.0;*/

            double nv = 0.0;
            double np = 0.0;
	    // Note: the following is only need for incompressible formulations. Otherwise, the adiabat is calculated self-consistently.
	    // Comment out the next line if using compressible.
            double Tadia = temperature+ 0.5 * depth/1e3; // Add the adiabatic gradient to the temperature for the viscosity calculation

	    //double UM_B1 = 0.0;             // Calc viscosities of each material component, then average (arithmetic)
            std::vector<double> composition_viscosities(volume_fractions.size());
	    //MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;
	    //viscosity_averaging = MaterialUtilities::arithmetic;
            
            // YIELDING - Byerlee's law approach
            // 
            /* double tauy = tau_0; // Small extra amount to prevent it being set to zero by input file - division later
            tauy = tau_0 + gamma*std::max(pressure,0.0);  
            tauy = std::min(yield_max,std::max(0.0,tauy));
            // Calculate the yielding viscosity
            veffp = tauy/(2.0*edot_ii);
            veffp = std::max(visc_min,std::min(visc_max,veffp)); */
     
            // Drucker Prager approach!
            const double angle_internal_friction = std::atan(gamma);
            const double sin_phi = std::sin(angle_internal_friction);
            const double cos_phi = std::cos(angle_internal_friction);
            const double stress_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

            const double cohesion = tau_0;
            double yield_stress = ( (dim==3)
                                ?
                                ( 6.0 * cohesion * cos_phi + 6.0 * std::max(pressure,0.0) * sin_phi) * stress_inv_part
                                :
                                cohesion * cos_phi + std::max(pressure,0.0) * sin_phi);
            
            const double strain_rate_effective_inv = 1./(2.*edot_ii);
            const double yield_stress1 = std::min(yield_stress, yield_max);
            veffp = yield_stress1 * strain_rate_effective_inv;
            //printf("\tCalculated yield\n");
            // End Drucker Prager approach    

	    for (unsigned int j=0; j < volume_fractions.size(); ++j) {
              /*UM_activation_energy += volume_fractions[j] * UM_activation_energies[j]; // NO Converted to J/mol to make units work
              LM_activation_energy += volume_fractions[j] * LM_activation_energies[j]; // NO Converted to J/mol to make units work
              nv += volume_fractions[j] * nvs[j];*/
              np += volume_fractions[j] * nps[j];
              veffv = visc_max;
	      double T2 = 0.0;
	      const double depth1 = depth; // Depths are in m in fit
	     //printf("I am here \n");
              if (Tadia != 0.0 && depth > 670e3) {
                // Homologous scaling approach
                //const double depth1 = depth/1e3;
                //const double depth1 = depth; // Depths are in m in fit
                // Fitted coeffs for Stixrude solidus: [ 1.77157968e+03  4.58499793e-04  3.27527226e-09 -2.81678963e-15  9.41736767e-22 -1.12730720e-28]
                try { 
		double Brig_sol = Brig_sol_a + Brig_sol_b*depth1 + Brig_sol_c*depth1*depth1 + Brig_sol_d*std::pow(depth1,3) + Brig_sol_e*std::pow(depth1,4) + Brig_sol_f*std::pow(depth1,5);
                // Fitted coeffs for MgO solidus:
                double Ferr_sol = Ferr_sol_a + Ferr_sol_b*depth1 + Ferr_sol_c*depth1*depth1 + Ferr_sol_d*std::pow(depth1,3)+ Ferr_sol_e*std::pow(depth1,4)+ Ferr_sol_f*std::pow(depth1,5)+ Ferr_sol_g*std::pow(depth1,6);
		// UM solidus coeffs (Stixrude): [ 1.49419506e+03  3.58738385e-03 -5.24082197e-09  2.91793473e-15]  
		//printf("\tCalculated solidus %.2f %.2e \n",depth1,Ferr_sol);
	  	if (Brig_sol < Brig_sol_a) {
			Brig_sol = Brig_sol_a;}
		if (Brig_sol > 5500) {
		  	Brig_sol = 5500;
		}
		if (Ferr_sol < Ferr_sol_a) {
                        Ferr_sol = Ferr_sol_a;}
                if (Ferr_sol > 4800) {
                        Ferr_sol = 4800;
                }	
                // For Brigmanite
                double Brig_Dff = 0.0;
                double Brig_veff = 1.e21;
                if (Tadia != 0.0) {
		 T2 = std::max(std::min(Tadia,Brig_sol),900.); //Limit the temperature extremes to avoid mishaps 
                 Brig_Dff = Brig_D0 * std::exp(-Brig_g0*Brig_sol/T2);
                  if ((Brig_gr != 0.0) && (T2 > 0)) {
                  Brig_veff = 1.0 / (13.3 * (Brig_om/(R * T2)) * (Brig_Dff/(Brig_gr*Brig_gr))); 
                  }
                }
               //printf("\tCalculated B %.2e \n",Brig_veff);
                // For Ferropericlase
                double Ferr_Dff = 0.0;
                double Ferr_veff = 1.e21;
               if (Tadia != 0.0) {
		  T2 = std::max(std::min(Ferr_sol,Tadia),900.);     
                  Ferr_Dff = Ferr_D0 * std::exp(-Ferr_g0*Ferr_sol/T2);
                 if ( (Ferr_gr != 0.0) && (T2>0)) {
                  Ferr_veff = 1.0 / (13.3 * (Ferr_om/(R * T2)) * (Ferr_Dff/(Ferr_gr*Ferr_gr))); 
                  }
                }
                //printf("\tCalculated FP %.2e\n",Ferr_veff);
                // Mixing law - Ji 2004 GMR, Takeda and Griera, 2006
                double brig_pc = composition[0];
                brig_pc = std::min(1.,std::max(brig_pc,0.));  // Keep in bounds of 0->1.
                try {
			veffv = std::pow(((1. - brig_pc)*std::pow(Ferr_veff,Ji)  +  brig_pc*std::pow(Brig_veff,Ji)
),1./Ji);
		}
		catch (...) {
			veffv = 1e23;
		}
                //printf("\tCalculated mix %.2e %.1f \n",veffv,brig_pc);
              }
		catch (...) {
		  veffv = 1e23;
		}
		}
               // Calculate T-dep viscosity for upper and lower mantle
              if (Tadia != 0.0 && depth <= 670e3) {
               	const double UM_Tsol = UM_sol_a + UM_sol_b*depth1 + UM_sol_c*depth1*depth1 + UM_sol_d*std::pow(depth1,3);
		T2 = std::min(Tadia,UM_Tsol);
		veffv =  (UM_B[j] * std::exp((UM_activation_energies[j] + UM_activation_volume*std::max(pressure,0.0))/(nvs[j]*R*T2)));
              }

              //printf("\t Calculated UM %.2e \n",veffv);
              //if (temperature != 0.0 && depth > 670e3) {
              //  veffv =  (LM_B * std::exp((LM_activation_energies[j] + LM_activation_volume*std::max(pressure,0.0))/(nvs[j]*R*temperature)));
              //}
              veffv = std::max(visc_min,std::min(visc_max,veffv));
              double viscous_stress = 2. * veffv * edot_ii;

              //printf("j %d \n",j);
                       
              // Effective viscosity = harmonic mean of diffusion and dislocation creep.
              // = (1/(v_eff^v) + 1/(v_eff^p))^-1

              double veff = veffv;
              /*if (veffv != 0) {
               veff = std::pow((1.0/veffv + 1.0/veffp), -1.0);
              } */
              if (viscous_stress >= yield_stress) {
                veff = veffp;
              }
              veff = std::max(visc_min,std::min(visc_max,veff));
              composition_viscosities[j] = veff;
            }
            //printf("Done visc \n");
            // 
            out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, composition_viscosities, viscosity_averaging);
          }
          
          for (unsigned int c=0; c<in.composition[i].size(); ++c) {
                  out.reaction_terms[i][c] = 0.000000;
              }

          // Calculate density

          {
            double density = 0.0;
            //printf("Start density");
	    for (unsigned int j=0; j < volume_fractions.size(); ++j)
              {
                // not strictly correct if thermal expansivities are different, since we are interpreting
                // these compositions as volume fractions, but the error introduced should not be too bad.
                double brig_pc = composition[0];
		double Bdg_DenFac= 1.0;
                if (depth > 670.e3) {
                	brig_pc = std::min(1.,std::max(brig_pc,0.)); 
			Bdg_DenFac =  1 + Bdg_Dens_Multiplier * brig_pc; // Bridgmanite density multiplier - should be
                // den = Den_FeO + (Den_Bdg - Den_FeO)* Bdg_PC
                // den/den_FeO = 1 + dRho/den_Feo * Bdg_pc
                // Brd_Dens_Multiplier = (Den_Bdg - Den_FeO)/den_FeO
                // density = den_FeO (default density)
                }
		const double temperature_factor = (1.0 - thermal_expansivities[j] * (temperature - reference_T));
                density += volume_fractions[j] * densities[j] * temperature_factor * Bdg_DenFac;
              }
            out.densities[i] = density;
            //printf("\t density done\n");
	   }

          // Thermal expansion coefficients at the given positions.
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];
	  // Depth dependent thermal expansion
	  thermal_expansivity *= (2890.e3 - exps_decrease*depth)/2890.e3;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = heat_capacity;
          // Thermal conductivity at the given positions.
          double therm_increase = 1.0;
	  if (depth > 670.e3) {
		therm_increase = 1.5 * (2890.e3 + depth - 670.e3)/(2890e3 - 670e3);
		}
	  out.thermal_conductivities[i] = therm_increase * thermal_diffusivity * heat_capacity * out.densities[i];
          
		// Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
        //printf("Max Veffp %.2e tauy %.2e edot %.2e || Min Veffp %.2e tauy %.2e edot %.2e \n",vmax, tauy_max,edot_max,vmin,tauy_min,edot_min);
          
          ii++;

      }  // Closes for loop
    //printf(" ##################################################### Finished evaluate #############################################\n"); 
    } // Closes evaluate

    template <int dim>
    double
    MantleHack<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    bool
    MantleHack<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    MantleHack<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("MantleHack");
        {
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J/kg/K$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $1 / K$");
          // Viscosity terms
          prm.declare_entry ("Upper mantle Activation energies", "500",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kJ / mol$");
          prm.declare_entry ("Lower mantle Activation energies", "500",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $kJ / mol$");

          prm.declare_entry ("Stress exponents for viscous rheology", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_v$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: None");
          prm.declare_entry ("Stress exponents for plastic rheology", "30",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_p$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: None"); 
          prm.declare_entry ("Upper mantle Activation volume", "6.4e-6", Patterns::Double(0), "($V_a$). Units: $m^3 / mol$");
          prm.declare_entry ("Lower mantle Activation volume", "6.4e-6", Patterns::Double(0), "($V_a$). Units: $m^3 / mol$");
          //prm.declare_entry ("UM Preexponential constant for viscous rheology law", "1.24e14", Patterns::Double(0), "($UM_B$). Units: None");
          prm.declare_entry ("UM Preexponential constant for viscous rheology law", "1.24e24",
                             Patterns::List(Patterns::Double(0)),
                             "List of Viscosity pre-exponents $UM_B$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: None"); 


          prm.declare_entry ("LM Preexponential constant for viscous rheology law", "1.24e14", Patterns::Double(0), "($LM_B$). Units: None");
        
          prm.declare_entry ("Brig_sol_a", "0.8e-6", Patterns::Double(0), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Brig_sol_b", "0.8e-6", Patterns::Double(0), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Brig_sol_c", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Brig_sol_d", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Brig_sol_e", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4 + eT^5 + fT^-6");
          prm.declare_entry ("Brig_sol_f", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
 
          prm.declare_entry ("Ferr_sol_a", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT+ eT^5 + fT^-64");
          prm.declare_entry ("Ferr_sol_b", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Ferr_sol_c", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Ferr_sol_d", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Ferr_sol_e", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Ferr_sol_f", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("Ferr_sol_g", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");


	  prm.declare_entry ("UM_sol_a", "0.8e-6", Patterns::Double(0), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("UM_sol_b", "0.8e-6", Patterns::Double(0), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("UM_sol_c", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");
          prm.declare_entry ("UM_sol_d", "0.8e-6", Patterns::Double(), "Solidus for homologous scaling. Ts = a + bT + cT^3 + dT4+ eT^5 + fT^-6");

          prm.declare_entry ("Brig_D0", "8.0e-7", Patterns::Double(0), "D0 for diffusion equation D = D0 * np.exp((-g0*Tsol)/T)");
          prm.declare_entry ("Brig_g0", "28.0", Patterns::Double(0), "g0 for diffusion equation D = D0 * np.exp((-g0*Tsol)/T)");
          prm.declare_entry ("Brig_om", "43.67e-6", Patterns::Double(0), "om for viscosity equation visc = 1/(13.3 * om/(k*T) * (D/(gr*gr)))");
          prm.declare_entry ("Brig_gr", "1.0e-3", Patterns::Double(0), "gr for viscosity equation visc = 1/(13.3 * om/(k*T) * (D/(gr*gr)))");
          
          prm.declare_entry ("Ferr_D0", "8.0e-7", Patterns::Double(0), "D0 for diffusion equation D = D0 * np.exp((-g0*Tsol)/T)");
          prm.declare_entry ("Ferr_g0", "28.0", Patterns::Double(0), "g0 for diffusion equation D = D0 * np.exp((-g0*Tsol)/T)");
          prm.declare_entry ("Ferr_om", "43.67e-6", Patterns::Double(0), "om for viscosity equation visc = 1/(13.3 * om/(k*T) * (D/(gr*gr)))");
          prm.declare_entry ("Ferr_gr", "1.0e-3", Patterns::Double(0), "gr for viscosity equation visc = 1/(13.3 * om/(k*T) * (D/(gr*gr)))");

          prm.declare_entry ("Ji","-0.25",Patterns::Double(), "Ji mixing rule exponent");
          prm.declare_entry ("Bdg_Dens_multiplier", "1.0", Patterns::Double(), "Bridgmanite density multiplier");
                   
          prm.declare_entry ("visc_min", "1.0e19", Patterns::Double(0), "Viscosity minimum cutoff");
          prm.declare_entry ("visc_max", "1.0e22", Patterns::Double(0), "Viscosity maximum cutoff");
          
	  prm.declare_entry ("exps_decrease","0.9",Patterns::Double(), "Expansivity decrease with mantle depth (2890km) -> 0.9 goes from alpha to 0.1 alpha");

	  prm.declare_entry ("Viscosity averaging scheme", "arithmetic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");

          prm.declare_entry ("Coefficient of yield stress increase with pressure", "0.25", Patterns::Double(0), "($\\gamma$). Units: None");
          prm.declare_entry ("Cohesive strength of rocks at the surface", "117", Patterns::Double(0), "($\\tau_0$). Units: $Pa$");
          prm.declare_entry ("Yield stress max", "10e9", Patterns::Double(0), "($\\yield_max$). Units: $Pa$");
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0), "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0), "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Reference strain rate", "6.4e-16", Patterns::Double(0), "($\\dot{\\epsilon}_\\text{ref}$). Units: $1 / s$");
          prm.declare_entry ("Reference viscosity", "1e21", Patterns::Double(0), "Reference viscosity for nondimensionalization.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MantleHack<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("MantleHack");
        {
          gamma = prm.get_double("Coefficient of yield stress increase with pressure");
          tau_0 = prm.get_double("Cohesive strength of rocks at the surface");
          yield_max = prm.get_double("Yield stress max");
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_strain_rate = prm.get_double("Minimum strain rate");
          reference_T = prm.get_double("Reference temperature");

          ref_visc = prm.get_double ("Reference viscosity");
          visc_min = prm.get_double("visc_min");
          visc_max = prm.get_double("visc_max");
	
          Brig_sol_a = prm.get_double("Brig_sol_a");
          Brig_sol_b = prm.get_double("Brig_sol_b");
          Brig_sol_c = prm.get_double("Brig_sol_c");
          Brig_sol_d = prm.get_double("Brig_sol_d");
          Brig_sol_e = prm.get_double("Brig_sol_e");
          Brig_sol_f = prm.get_double("Brig_sol_f");


          Ferr_sol_a = prm.get_double("Ferr_sol_a");
          Ferr_sol_b = prm.get_double("Ferr_sol_b");
          Ferr_sol_c = prm.get_double("Ferr_sol_c");
          Ferr_sol_d = prm.get_double("Ferr_sol_d");
          Ferr_sol_e = prm.get_double("Ferr_sol_e");
          Ferr_sol_f = prm.get_double("Ferr_sol_f");
          Ferr_sol_g = prm.get_double("Ferr_sol_g");
 
	  UM_sol_a = prm.get_double("Brig_sol_a");
          UM_sol_b = prm.get_double("Brig_sol_b");
          UM_sol_c = prm.get_double("Brig_sol_c");
          UM_sol_d = prm.get_double("Brig_sol_d");

          Brig_D0 = prm.get_double("Brig_D0");
          Brig_g0 = prm.get_double("Brig_g0");
          Brig_om = prm.get_double("Brig_om");
          Brig_gr = prm.get_double("Brig_gr");

          Ferr_D0 = prm.get_double("Ferr_D0");
          Ferr_g0 = prm.get_double("Ferr_g0");
          Ferr_om = prm.get_double("Ferr_om");
          Ferr_gr = prm.get_double("Ferr_gr");

          Ji = prm.get_double("Ji");
          Bdg_Dens_Multiplier = prm.get_double("Bdg_Dens_multiplier");
	  exps_decrease = prm.get_double("exps_decrease");

	  viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                prm);

          std::vector<double> x_values;

          // Parse densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of density list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            densities.assign( n_fields , x_values[0]);
          else
            densities = x_values;
        
          // Viscosity stuff
          UM_activation_volume = prm.get_double("Upper mantle Activation volume");
          LM_activation_volume = prm.get_double("Lower mantle Activation volume");
          //UM_B = prm.get_double("UM Preexponential constant for viscous rheology law");
          LM_B = prm.get_double("LM Preexponential constant for viscous rheology law");
          // Parse UM_B
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("UM Preexponential constant for viscous rheology law")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of UM Preexponential constant list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            UM_B.assign( n_fields , x_values[0] );
          else
            UM_B = x_values; 

          // Parse activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Upper mantle Activation energies")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            UM_activation_energies.assign( n_fields , x_values[0] );
          else
            UM_activation_energies = x_values; 

          // Parse activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Lower mantle Activation energies")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            LM_activation_energies.assign( n_fields , x_values[0] );
          else
            LM_activation_energies = x_values;

          // Parse thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of thermal expansivity list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            thermal_expansivities.assign( n_fields , x_values[0]);
          else
            thermal_expansivities = x_values;

          // Parse stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for viscous rheology")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of nv list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            nvs.assign( n_fields , x_values[0]);
          else
            nvs = x_values;
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for plastic rheology")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of np list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            nps.assign( n_fields , x_values[0]);
          else
            nps = x_values;


        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::compositional_fields;
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MantleHack,
                                   "MantleHack",
                                   "Mostly Craig just &^*&ing around")
  }
}
