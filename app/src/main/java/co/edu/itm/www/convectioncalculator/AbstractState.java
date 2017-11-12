/*******************************************************************************
 * The MIT License (MIT)
 *
 * Convection calculator app
 * Copyright (C) 2017 Luis F. Cardona S.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.8
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package co.edu.itm.www.convectioncalculator;

public class AbstractState {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected AbstractState(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(AbstractState obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        CoolPropJNI.delete_AbstractState(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public static AbstractState factory(String backend, String fluid_names) {
    long cPtr = CoolPropJNI.AbstractState_factory__SWIG_0(backend, fluid_names);
    return (cPtr == 0) ? null : new AbstractState(cPtr, false);
  }

  public static AbstractState factory(String backend, StringVector fluid_names) {
    long cPtr = CoolPropJNI.AbstractState_factory__SWIG_1(backend, StringVector.getCPtr(fluid_names), fluid_names);
    return (cPtr == 0) ? null : new AbstractState(cPtr, false);
  }

  public void set_T(double T) {
    CoolPropJNI.AbstractState_set_T(swigCPtr, this, T);
  }

  public String backend_name() {
    return CoolPropJNI.AbstractState_backend_name(swigCPtr, this);
  }

  public boolean using_mole_fractions() {
    return CoolPropJNI.AbstractState_using_mole_fractions(swigCPtr, this);
  }

  public boolean using_mass_fractions() {
    return CoolPropJNI.AbstractState_using_mass_fractions(swigCPtr, this);
  }

  public boolean using_volu_fractions() {
    return CoolPropJNI.AbstractState_using_volu_fractions(swigCPtr, this);
  }

  public void set_reference_stateS(String reference_state) {
    CoolPropJNI.AbstractState_set_reference_stateS(swigCPtr, this, reference_state);
  }

  public void set_reference_stateD(double T, double rhomolar, double hmolar0, double smolar0) {
    CoolPropJNI.AbstractState_set_reference_stateD(swigCPtr, this, T, rhomolar, hmolar0, smolar0);
  }

  public void set_mole_fractions(DoubleVector mole_fractions) {
    CoolPropJNI.AbstractState_set_mole_fractions(swigCPtr, this, DoubleVector.getCPtr(mole_fractions), mole_fractions);
  }

  public void set_mass_fractions(DoubleVector mass_fractions) {
    CoolPropJNI.AbstractState_set_mass_fractions(swigCPtr, this, DoubleVector.getCPtr(mass_fractions), mass_fractions);
  }

  public void set_volu_fractions(DoubleVector volu_fractions) {
    CoolPropJNI.AbstractState_set_volu_fractions(swigCPtr, this, DoubleVector.getCPtr(volu_fractions), volu_fractions);
  }

  public DoubleVector mole_fractions_liquid() {
    return new DoubleVector(CoolPropJNI.AbstractState_mole_fractions_liquid(swigCPtr, this), true);
  }

  public DoubleVector mole_fractions_vapor() {
    return new DoubleVector(CoolPropJNI.AbstractState_mole_fractions_vapor(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_CoolPropDbl_t get_mole_fractions() {
    return new SWIGTYPE_p_std__vectorT_CoolPropDbl_t(CoolPropJNI.AbstractState_get_mole_fractions(swigCPtr, this), false);
  }

  public SWIGTYPE_p_std__vectorT_CoolPropDbl_t get_mass_fractions() {
    return new SWIGTYPE_p_std__vectorT_CoolPropDbl_t(CoolPropJNI.AbstractState_get_mass_fractions(swigCPtr, this), true);
  }

  public void update(input_pairs input_pair, double Value1, double Value2) {
    CoolPropJNI.AbstractState_update(swigCPtr, this, input_pair.swigValue(), Value1, Value2);
  }

  public void update_with_guesses(input_pairs input_pair, double Value1, double Value2, GuessesStructure guesses) {
    CoolPropJNI.AbstractState_update_with_guesses(swigCPtr, this, input_pair.swigValue(), Value1, Value2, GuessesStructure.getCPtr(guesses), guesses);
  }

  public boolean available_in_high_level() {
    return CoolPropJNI.AbstractState_available_in_high_level(swigCPtr, this);
  }

  public String fluid_param_string(String arg0) {
    return CoolPropJNI.AbstractState_fluid_param_string(swigCPtr, this, arg0);
  }

  public StringVector fluid_names() {
    return new StringVector(CoolPropJNI.AbstractState_fluid_names(swigCPtr, this), true);
  }

  public double get_fluid_constant(long i, parameters param) {
    return CoolPropJNI.AbstractState_get_fluid_constant(swigCPtr, this, i, param.swigValue());
  }

  public void set_binary_interaction_double(String CAS1, String CAS2, String parameter, double value) {
    CoolPropJNI.AbstractState_set_binary_interaction_double__SWIG_0(swigCPtr, this, CAS1, CAS2, parameter, value);
  }

  public void set_binary_interaction_double(long i, long j, String parameter, double value) {
    CoolPropJNI.AbstractState_set_binary_interaction_double__SWIG_1(swigCPtr, this, i, j, parameter, value);
  }

  public void set_binary_interaction_string(String CAS1, String CAS2, String parameter, String value) {
    CoolPropJNI.AbstractState_set_binary_interaction_string__SWIG_0(swigCPtr, this, CAS1, CAS2, parameter, value);
  }

  public void set_binary_interaction_string(long i, long j, String parameter, String value) {
    CoolPropJNI.AbstractState_set_binary_interaction_string__SWIG_1(swigCPtr, this, i, j, parameter, value);
  }

  public double get_binary_interaction_double(String CAS1, String CAS2, String parameter) {
    return CoolPropJNI.AbstractState_get_binary_interaction_double__SWIG_0(swigCPtr, this, CAS1, CAS2, parameter);
  }

  public double get_binary_interaction_double(long i, long j, String parameter) {
    return CoolPropJNI.AbstractState_get_binary_interaction_double__SWIG_1(swigCPtr, this, i, j, parameter);
  }

  public String get_binary_interaction_string(String CAS1, String CAS2, String parameter) {
    return CoolPropJNI.AbstractState_get_binary_interaction_string(swigCPtr, this, CAS1, CAS2, parameter);
  }

  public void apply_simple_mixing_rule(long i, long j, String model) {
    CoolPropJNI.AbstractState_apply_simple_mixing_rule(swigCPtr, this, i, j, model);
  }

  public void set_cubic_alpha_C(long i, String parameter, double c1, double c2, double c3) {
    CoolPropJNI.AbstractState_set_cubic_alpha_C(swigCPtr, this, i, parameter, c1, c2, c3);
  }

  public void set_fluid_parameter_double(long i, String parameter, double value) {
    CoolPropJNI.AbstractState_set_fluid_parameter_double(swigCPtr, this, i, parameter, value);
  }

  public boolean clear() {
    return CoolPropJNI.AbstractState_clear(swigCPtr, this);
  }

  public SimpleState get_reducing_state() {
    return new SimpleState(CoolPropJNI.AbstractState_get_reducing_state(swigCPtr, this), false);
  }

  public SimpleState get_state(String state) {
    return new SimpleState(CoolPropJNI.AbstractState_get_state(swigCPtr, this, state), false);
  }

  public double Tmin() {
    return CoolPropJNI.AbstractState_Tmin(swigCPtr, this);
  }

  public double Tmax() {
    return CoolPropJNI.AbstractState_Tmax(swigCPtr, this);
  }

  public double pmax() {
    return CoolPropJNI.AbstractState_pmax(swigCPtr, this);
  }

  public double Ttriple() {
    return CoolPropJNI.AbstractState_Ttriple(swigCPtr, this);
  }

  public phases phase() {
    return phases.swigToEnum(CoolPropJNI.AbstractState_phase(swigCPtr, this));
  }

  public void specify_phase(phases phase) {
    CoolPropJNI.AbstractState_specify_phase(swigCPtr, this, phase.swigValue());
  }

  public void unspecify_phase() {
    CoolPropJNI.AbstractState_unspecify_phase(swigCPtr, this);
  }

  public double T_critical() {
    return CoolPropJNI.AbstractState_T_critical(swigCPtr, this);
  }

  public double p_critical() {
    return CoolPropJNI.AbstractState_p_critical(swigCPtr, this);
  }

  public double rhomolar_critical() {
    return CoolPropJNI.AbstractState_rhomolar_critical(swigCPtr, this);
  }

  public double rhomass_critical() {
    return CoolPropJNI.AbstractState_rhomass_critical(swigCPtr, this);
  }

  public SWIGTYPE_p_std__vectorT_CoolProp__CriticalState_t all_critical_points() {
    return new SWIGTYPE_p_std__vectorT_CoolProp__CriticalState_t(CoolPropJNI.AbstractState_all_critical_points(swigCPtr, this), true);
  }

  public void build_spinodal() {
    CoolPropJNI.AbstractState_build_spinodal(swigCPtr, this);
  }

  public SpinodalData get_spinodal_data() {
    return new SpinodalData(CoolPropJNI.AbstractState_get_spinodal_data(swigCPtr, this), true);
  }

  public void criticality_contour_values(SWIGTYPE_p_double L1star, SWIGTYPE_p_double M1star) {
    CoolPropJNI.AbstractState_criticality_contour_values(swigCPtr, this, SWIGTYPE_p_double.getCPtr(L1star), SWIGTYPE_p_double.getCPtr(M1star));
  }

  public double tangent_plane_distance(double T, double p, DoubleVector w, double rhomolar_guess) {
    return CoolPropJNI.AbstractState_tangent_plane_distance__SWIG_0(swigCPtr, this, T, p, DoubleVector.getCPtr(w), w, rhomolar_guess);
  }

  public double tangent_plane_distance(double T, double p, DoubleVector w) {
    return CoolPropJNI.AbstractState_tangent_plane_distance__SWIG_1(swigCPtr, this, T, p, DoubleVector.getCPtr(w), w);
  }

  public double T_reducing() {
    return CoolPropJNI.AbstractState_T_reducing(swigCPtr, this);
  }

  public double rhomolar_reducing() {
    return CoolPropJNI.AbstractState_rhomolar_reducing(swigCPtr, this);
  }

  public double rhomass_reducing() {
    return CoolPropJNI.AbstractState_rhomass_reducing(swigCPtr, this);
  }

  public double p_triple() {
    return CoolPropJNI.AbstractState_p_triple(swigCPtr, this);
  }

  public String name() {
    return CoolPropJNI.AbstractState_name(swigCPtr, this);
  }

  public double dipole_moment() {
    return CoolPropJNI.AbstractState_dipole_moment(swigCPtr, this);
  }

  public double keyed_output(parameters key) {
    return CoolPropJNI.AbstractState_keyed_output(swigCPtr, this, key.swigValue());
  }

  public double trivial_keyed_output(parameters key) {
    return CoolPropJNI.AbstractState_trivial_keyed_output(swigCPtr, this, key.swigValue());
  }

  public double saturated_liquid_keyed_output(parameters key) {
    return CoolPropJNI.AbstractState_saturated_liquid_keyed_output(swigCPtr, this, key.swigValue());
  }

  public double saturated_vapor_keyed_output(parameters key) {
    return CoolPropJNI.AbstractState_saturated_vapor_keyed_output(swigCPtr, this, key.swigValue());
  }

  public double T() {
    return CoolPropJNI.AbstractState_T(swigCPtr, this);
  }

  public double rhomolar() {
    return CoolPropJNI.AbstractState_rhomolar(swigCPtr, this);
  }

  public double rhomass() {
    return CoolPropJNI.AbstractState_rhomass(swigCPtr, this);
  }

  public double p() {
    return CoolPropJNI.AbstractState_p(swigCPtr, this);
  }

  public double Q() {
    return CoolPropJNI.AbstractState_Q(swigCPtr, this);
  }

  public double tau() {
    return CoolPropJNI.AbstractState_tau(swigCPtr, this);
  }

  public double delta() {
    return CoolPropJNI.AbstractState_delta(swigCPtr, this);
  }

  public double molar_mass() {
    return CoolPropJNI.AbstractState_molar_mass(swigCPtr, this);
  }

  public double acentric_factor() {
    return CoolPropJNI.AbstractState_acentric_factor(swigCPtr, this);
  }

  public double gas_constant() {
    return CoolPropJNI.AbstractState_gas_constant(swigCPtr, this);
  }

  public double Bvirial() {
    return CoolPropJNI.AbstractState_Bvirial(swigCPtr, this);
  }

  public double dBvirial_dT() {
    return CoolPropJNI.AbstractState_dBvirial_dT(swigCPtr, this);
  }

  public double Cvirial() {
    return CoolPropJNI.AbstractState_Cvirial(swigCPtr, this);
  }

  public double dCvirial_dT() {
    return CoolPropJNI.AbstractState_dCvirial_dT(swigCPtr, this);
  }

  public double compressibility_factor() {
    return CoolPropJNI.AbstractState_compressibility_factor(swigCPtr, this);
  }

  public double hmolar() {
    return CoolPropJNI.AbstractState_hmolar(swigCPtr, this);
  }

  public double hmass() {
    return CoolPropJNI.AbstractState_hmass(swigCPtr, this);
  }

  public double hmolar_excess() {
    return CoolPropJNI.AbstractState_hmolar_excess(swigCPtr, this);
  }

  public double hmass_excess() {
    return CoolPropJNI.AbstractState_hmass_excess(swigCPtr, this);
  }

  public double smolar() {
    return CoolPropJNI.AbstractState_smolar(swigCPtr, this);
  }

  public double smass() {
    return CoolPropJNI.AbstractState_smass(swigCPtr, this);
  }

  public double smolar_excess() {
    return CoolPropJNI.AbstractState_smolar_excess(swigCPtr, this);
  }

  public double smass_excess() {
    return CoolPropJNI.AbstractState_smass_excess(swigCPtr, this);
  }

  public double umolar() {
    return CoolPropJNI.AbstractState_umolar(swigCPtr, this);
  }

  public double umass() {
    return CoolPropJNI.AbstractState_umass(swigCPtr, this);
  }

  public double umolar_excess() {
    return CoolPropJNI.AbstractState_umolar_excess(swigCPtr, this);
  }

  public double umass_excess() {
    return CoolPropJNI.AbstractState_umass_excess(swigCPtr, this);
  }

  public double cpmolar() {
    return CoolPropJNI.AbstractState_cpmolar(swigCPtr, this);
  }

  public double cpmass() {
    return CoolPropJNI.AbstractState_cpmass(swigCPtr, this);
  }

  public double cp0molar() {
    return CoolPropJNI.AbstractState_cp0molar(swigCPtr, this);
  }

  public double cp0mass() {
    return CoolPropJNI.AbstractState_cp0mass(swigCPtr, this);
  }

  public double cvmolar() {
    return CoolPropJNI.AbstractState_cvmolar(swigCPtr, this);
  }

  public double cvmass() {
    return CoolPropJNI.AbstractState_cvmass(swigCPtr, this);
  }

  public double gibbsmolar() {
    return CoolPropJNI.AbstractState_gibbsmolar(swigCPtr, this);
  }

  public double gibbsmass() {
    return CoolPropJNI.AbstractState_gibbsmass(swigCPtr, this);
  }

  public double gibbsmolar_excess() {
    return CoolPropJNI.AbstractState_gibbsmolar_excess(swigCPtr, this);
  }

  public double gibbsmass_excess() {
    return CoolPropJNI.AbstractState_gibbsmass_excess(swigCPtr, this);
  }

  public double helmholtzmolar() {
    return CoolPropJNI.AbstractState_helmholtzmolar(swigCPtr, this);
  }

  public double helmholtzmass() {
    return CoolPropJNI.AbstractState_helmholtzmass(swigCPtr, this);
  }

  public double helmholtzmolar_excess() {
    return CoolPropJNI.AbstractState_helmholtzmolar_excess(swigCPtr, this);
  }

  public double helmholtzmass_excess() {
    return CoolPropJNI.AbstractState_helmholtzmass_excess(swigCPtr, this);
  }

  public double volumemolar_excess() {
    return CoolPropJNI.AbstractState_volumemolar_excess(swigCPtr, this);
  }

  public double volumemass_excess() {
    return CoolPropJNI.AbstractState_volumemass_excess(swigCPtr, this);
  }

  public double speed_sound() {
    return CoolPropJNI.AbstractState_speed_sound(swigCPtr, this);
  }

  public double isothermal_compressibility() {
    return CoolPropJNI.AbstractState_isothermal_compressibility(swigCPtr, this);
  }

  public double isobaric_expansion_coefficient() {
    return CoolPropJNI.AbstractState_isobaric_expansion_coefficient(swigCPtr, this);
  }

  public double fugacity_coefficient(long i) {
    return CoolPropJNI.AbstractState_fugacity_coefficient(swigCPtr, this, i);
  }

  public double fugacity(long i) {
    return CoolPropJNI.AbstractState_fugacity(swigCPtr, this, i);
  }

  public double chemical_potential(long i) {
    return CoolPropJNI.AbstractState_chemical_potential(swigCPtr, this, i);
  }

  public double fundamental_derivative_of_gas_dynamics() {
    return CoolPropJNI.AbstractState_fundamental_derivative_of_gas_dynamics(swigCPtr, this);
  }

  public double PIP() {
    return CoolPropJNI.AbstractState_PIP(swigCPtr, this);
  }

  public void true_critical_point(SWIGTYPE_p_double T, SWIGTYPE_p_double rho) {
    CoolPropJNI.AbstractState_true_critical_point(swigCPtr, this, SWIGTYPE_p_double.getCPtr(T), SWIGTYPE_p_double.getCPtr(rho));
  }

  public void ideal_curve(String type, DoubleVector T, DoubleVector p) {
    CoolPropJNI.AbstractState_ideal_curve(swigCPtr, this, type, DoubleVector.getCPtr(T), T, DoubleVector.getCPtr(p), p);
  }

  public double first_partial_deriv(parameters Of, parameters Wrt, parameters Constant) {
    return CoolPropJNI.AbstractState_first_partial_deriv(swigCPtr, this, Of.swigValue(), Wrt.swigValue(), Constant.swigValue());
  }

  public double second_partial_deriv(parameters Of1, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) {
    return CoolPropJNI.AbstractState_second_partial_deriv(swigCPtr, this, Of1.swigValue(), Wrt1.swigValue(), Constant1.swigValue(), Wrt2.swigValue(), Constant2.swigValue());
  }

  public double first_saturation_deriv(parameters Of1, parameters Wrt1) {
    return CoolPropJNI.AbstractState_first_saturation_deriv(swigCPtr, this, Of1.swigValue(), Wrt1.swigValue());
  }

  public double second_saturation_deriv(parameters Of1, parameters Wrt1, parameters Wrt2) {
    return CoolPropJNI.AbstractState_second_saturation_deriv(swigCPtr, this, Of1.swigValue(), Wrt1.swigValue(), Wrt2.swigValue());
  }

  public double first_two_phase_deriv(parameters Of, parameters Wrt, parameters Constant) {
    return CoolPropJNI.AbstractState_first_two_phase_deriv(swigCPtr, this, Of.swigValue(), Wrt.swigValue(), Constant.swigValue());
  }

  public double second_two_phase_deriv(parameters Of, parameters Wrt1, parameters Constant1, parameters Wrt2, parameters Constant2) {
    return CoolPropJNI.AbstractState_second_two_phase_deriv(swigCPtr, this, Of.swigValue(), Wrt1.swigValue(), Constant1.swigValue(), Wrt2.swigValue(), Constant2.swigValue());
  }

  public double first_two_phase_deriv_splined(parameters Of, parameters Wrt, parameters Constant, double x_end) {
    return CoolPropJNI.AbstractState_first_two_phase_deriv_splined(swigCPtr, this, Of.swigValue(), Wrt.swigValue(), Constant.swigValue(), x_end);
  }

  public void build_phase_envelope(String type) {
    CoolPropJNI.AbstractState_build_phase_envelope__SWIG_0(swigCPtr, this, type);
  }

  public void build_phase_envelope() {
    CoolPropJNI.AbstractState_build_phase_envelope__SWIG_1(swigCPtr, this);
  }

  public PhaseEnvelopeData get_phase_envelope_data() {
    return new PhaseEnvelopeData(CoolPropJNI.AbstractState_get_phase_envelope_data(swigCPtr, this), false);
  }

  public boolean has_melting_line() {
    return CoolPropJNI.AbstractState_has_melting_line(swigCPtr, this);
  }

  public double melting_line(int param, int given, double value) {
    return CoolPropJNI.AbstractState_melting_line(swigCPtr, this, param, given, value);
  }

  public double saturation_ancillary(parameters param, int Q, parameters given, double value) {
    return CoolPropJNI.AbstractState_saturation_ancillary(swigCPtr, this, param.swigValue(), Q, given.swigValue(), value);
  }

  public double viscosity() {
    return CoolPropJNI.AbstractState_viscosity(swigCPtr, this);
  }

  public void viscosity_contributions(SWIGTYPE_p_CoolPropDbl dilute, SWIGTYPE_p_CoolPropDbl initial_density, SWIGTYPE_p_CoolPropDbl residual, SWIGTYPE_p_CoolPropDbl critical) {
    CoolPropJNI.AbstractState_viscosity_contributions(swigCPtr, this, SWIGTYPE_p_CoolPropDbl.getCPtr(dilute), SWIGTYPE_p_CoolPropDbl.getCPtr(initial_density), SWIGTYPE_p_CoolPropDbl.getCPtr(residual), SWIGTYPE_p_CoolPropDbl.getCPtr(critical));
  }

  public double conductivity() {
    return CoolPropJNI.AbstractState_conductivity(swigCPtr, this);
  }

  public void conductivity_contributions(SWIGTYPE_p_CoolPropDbl dilute, SWIGTYPE_p_CoolPropDbl initial_density, SWIGTYPE_p_CoolPropDbl residual, SWIGTYPE_p_CoolPropDbl critical) {
    CoolPropJNI.AbstractState_conductivity_contributions(swigCPtr, this, SWIGTYPE_p_CoolPropDbl.getCPtr(dilute), SWIGTYPE_p_CoolPropDbl.getCPtr(initial_density), SWIGTYPE_p_CoolPropDbl.getCPtr(residual), SWIGTYPE_p_CoolPropDbl.getCPtr(critical));
  }

  public double surface_tension() {
    return CoolPropJNI.AbstractState_surface_tension(swigCPtr, this);
  }

  public double Prandtl() {
    return CoolPropJNI.AbstractState_Prandtl(swigCPtr, this);
  }

  public void conformal_state(String reference_fluid, SWIGTYPE_p_CoolPropDbl T, SWIGTYPE_p_CoolPropDbl rhomolar) {
    CoolPropJNI.AbstractState_conformal_state(swigCPtr, this, reference_fluid, SWIGTYPE_p_CoolPropDbl.getCPtr(T), SWIGTYPE_p_CoolPropDbl.getCPtr(rhomolar));
  }

  public void change_EOS(long i, String EOS_name) {
    CoolPropJNI.AbstractState_change_EOS(swigCPtr, this, i, EOS_name);
  }

  public double alpha0() {
    return CoolPropJNI.AbstractState_alpha0(swigCPtr, this);
  }

  public double dalpha0_dDelta() {
    return CoolPropJNI.AbstractState_dalpha0_dDelta(swigCPtr, this);
  }

  public double dalpha0_dTau() {
    return CoolPropJNI.AbstractState_dalpha0_dTau(swigCPtr, this);
  }

  public double d2alpha0_dDelta2() {
    return CoolPropJNI.AbstractState_d2alpha0_dDelta2(swigCPtr, this);
  }

  public double d2alpha0_dDelta_dTau() {
    return CoolPropJNI.AbstractState_d2alpha0_dDelta_dTau(swigCPtr, this);
  }

  public double d2alpha0_dTau2() {
    return CoolPropJNI.AbstractState_d2alpha0_dTau2(swigCPtr, this);
  }

  public double d3alpha0_dTau3() {
    return CoolPropJNI.AbstractState_d3alpha0_dTau3(swigCPtr, this);
  }

  public double d3alpha0_dDelta_dTau2() {
    return CoolPropJNI.AbstractState_d3alpha0_dDelta_dTau2(swigCPtr, this);
  }

  public double d3alpha0_dDelta2_dTau() {
    return CoolPropJNI.AbstractState_d3alpha0_dDelta2_dTau(swigCPtr, this);
  }

  public double d3alpha0_dDelta3() {
    return CoolPropJNI.AbstractState_d3alpha0_dDelta3(swigCPtr, this);
  }

  public double alphar() {
    return CoolPropJNI.AbstractState_alphar(swigCPtr, this);
  }

  public double dalphar_dDelta() {
    return CoolPropJNI.AbstractState_dalphar_dDelta(swigCPtr, this);
  }

  public double dalphar_dTau() {
    return CoolPropJNI.AbstractState_dalphar_dTau(swigCPtr, this);
  }

  public double d2alphar_dDelta2() {
    return CoolPropJNI.AbstractState_d2alphar_dDelta2(swigCPtr, this);
  }

  public double d2alphar_dDelta_dTau() {
    return CoolPropJNI.AbstractState_d2alphar_dDelta_dTau(swigCPtr, this);
  }

  public double d2alphar_dTau2() {
    return CoolPropJNI.AbstractState_d2alphar_dTau2(swigCPtr, this);
  }

  public double d3alphar_dDelta3() {
    return CoolPropJNI.AbstractState_d3alphar_dDelta3(swigCPtr, this);
  }

  public double d3alphar_dDelta2_dTau() {
    return CoolPropJNI.AbstractState_d3alphar_dDelta2_dTau(swigCPtr, this);
  }

  public double d3alphar_dDelta_dTau2() {
    return CoolPropJNI.AbstractState_d3alphar_dDelta_dTau2(swigCPtr, this);
  }

  public double d3alphar_dTau3() {
    return CoolPropJNI.AbstractState_d3alphar_dTau3(swigCPtr, this);
  }

  public double d4alphar_dDelta4() {
    return CoolPropJNI.AbstractState_d4alphar_dDelta4(swigCPtr, this);
  }

  public double d4alphar_dDelta3_dTau() {
    return CoolPropJNI.AbstractState_d4alphar_dDelta3_dTau(swigCPtr, this);
  }

  public double d4alphar_dDelta2_dTau2() {
    return CoolPropJNI.AbstractState_d4alphar_dDelta2_dTau2(swigCPtr, this);
  }

  public double d4alphar_dDelta_dTau3() {
    return CoolPropJNI.AbstractState_d4alphar_dDelta_dTau3(swigCPtr, this);
  }

  public double d4alphar_dTau4() {
    return CoolPropJNI.AbstractState_d4alphar_dTau4(swigCPtr, this);
  }

}