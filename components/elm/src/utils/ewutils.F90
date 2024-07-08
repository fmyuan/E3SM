module ewutils
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing enhanced weathering shared utilities
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use elm_varcon  , only: log_keq_hco3, log_keq_co3
  use elm_varpar  , only: ncations
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: mass_to_mol
  public :: mass_to_meq
  public :: mass_to_logmol
  public :: logmol_to_mass
  public :: mol_to_mass
  public :: meq_to_mass
  public :: ph_to_hco3
  public :: hco3_to_co3
  public :: objective_solveq
  public :: solve_eq

contains

  !-----------------------------------------------------------------------
  ! Unit conversion utilities
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function mass_to_mol(mass_conc, molar_mass, h2o) result(mol_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from ELM-standard
    ! units (g m-3 soil or g m-3 soil s-1) to typical chemical reaction
    ! database units (mol kg-1 water) 
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: mass_conc ! g m-3 soil
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: h2o ! m3 m-3, volumetric soil water content
    real(r8) :: mol_conc ! mol kg-1 water

    mol_conc = mass_conc / molar_mass / h2o / 1000._r8
  end function mass_to_mol

  !-----------------------------------------------------------------------
  function mass_to_meq(mass_conc, valence, molar_mass, bd) result(meq_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from ELM-standard
    ! units (g m-3 soil or g m-3 soil s-1) to typical CEC reaction
    ! database units (meq 100 g-1 dry soil)
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: mass_conc ! g m-3 soil
    real(r8), intent(in) :: valence ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: bd ! kg cm-3, soil bulk density
    real(r8) :: meq_conc ! mol kg-1 water

    meq_conc = mass_conc * 1000._r8 * valence / molar_mass / 10._r8 / bd
  end function mass_to_meq

  !-----------------------------------------------------------------------
  function mass_to_logmol(mass_conc, molar_mass, h2o) result(log_mol_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from ELM-standard
    ! units (g m-3 soil or g m-3 soil s-1) to typical chemical reaction
    ! database units (mol kg-1 water). Log transformation. 
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: mass_conc ! g m-3 soil
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: h2o ! m3 m-3, volumetric soil water content
    real(r8) :: log_mol_conc ! mol kg-1 water

    log_mol_conc = log10(mass_conc) - log10(molar_mass) - log10(h2o) - 3._r8
  end function mass_to_logmol

  !-----------------------------------------------------------------------
  function logmol_to_mass(log_mol_conc, molar_mass, h2o) result(mass_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from typical chemical reaction
    ! database units (mol kg-1 water) to ELM-standard units (g m-3 soil or 
    ! g m-3 soil s-1). Log transformation. 
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: log_mol_conc ! mol kg-1 water
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: h2o ! m3 m-3, volumetric soil water content
    real(r8) :: mass_conc ! g m-3 soil

    mass_conc = 10**(log_mol_conc + 3._r8) * molar_mass * h2o
  end function logmol_to_mass

  !-----------------------------------------------------------------------
  function mol_to_mass(mol_conc, molar_mass, h2o) result(mass_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from typical chemical reaction
    ! database units (mol kg-1 water) to ELM-standard units (g m-3 soil or 
    ! g m-3 soil s-1). No log transformation. 
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: mol_conc ! mol kg-1 water
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: h2o ! m3 m-3, volumetric soil water content
    real(r8) :: mass_conc ! g m-3 soil

    mass_conc = mol_conc * molar_mass* h2o * 1e3_r8
  end function mol_to_mass

  !-----------------------------------------------------------------------
  function meq_to_mass(meq_conc, valence, molar_mass, bd) result(mass_conc)
    !
    ! !DESCRIPTION:
    ! Convert the concentration of cations/solids from typical CEC reaction
    ! database units (meq 100 g-1 dry soil) to ELM-standard
    ! units (g m-3 soil or g m-3 soil s-1)
    !
    ! !ARGUMENTS: 
    real(r8), intent(in) :: meq_conc ! mol kg-1 water
    real(r8), intent(in) :: valence ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: bd ! kg cm-3, soil bulk density
    real(r8) :: mass_conc ! g m-3 soil

    mass_conc = meq_conc * molar_mass * 10._r8 * bd / 1000._r8 / valence
  end function meq_to_mass

  !-----------------------------------------------------------------------
  ! CO2 dynamics
  ! -----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function ph_to_hco3(soil_ph, co2_atm) result(hco3_conc)
    !
    ! !DESCRIPTION:
    ! Use the CO2(g) +1.0000 H2O = + 1.0000 H+ + 1.0000 HCO3-
    ! reaction to calculate the HCO3- concentration in mol/kg water under
    ! given soil pH and gaseous CO2 concentration
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: soil_ph
    real(r8), intent(in) :: co2_atm ! partial pressure of CO2 in atm 
    real(r8) :: hco3_conc ! mol/kg water

    hco3_conc = 10**(log_keq_hco3 + soil_ph) * co2_atm
  end function ph_to_hco3

  !-----------------------------------------------------------------------
  function hco3_to_co3(hco3_conc, soil_ph) result(co3_conc)
    !
    ! !DESCRIPTION:
    ! Use the 1.0000 HCO3- = CO3-- +1.0000 H+
    ! reaction to calculate the HCO3- concentration in mol/kg water under
    ! given soil pH and gaseous CO2 concentration
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hco3_conc ! mol/kg water
    real(r8), intent(in) :: soil_ph
    real(r8) :: co3_conc ! mol/kg water

    co3_conc = 10**(log_keq_co3 + soil_ph) * hco3_conc

  end function hco3_to_co3

  !-----------------------------------------------------------------------
  ! Functions related to the solution of dynamic pH
  !-----------------------------------------------------------------------

  function objective_solveq(soil_ph, b0, co2_atm, beta_list, kex_list, cation_valence) result (pcterr)
    !
    ! !DESCRIPTION:
    ! Calculate whether a given pH value satisfies the following set of equations
    ! 
    ! eq1 = sp.Eq(h * hco3 / co2_atm, 10**(-7.8136))
    ! eq2 = sp.Eq(h * co3 / hco3, 10**(-10.3288))
    ! eq3 = sp.Eq(h * oh, 1e-14)
    ! eq4 = sp.Eq(h / beta_h * (beta1 / ca)**(1/valence_Ca2), kex1) # 10**(3.4*(1-beta_h)) *  
    ! eq5 = sp.Eq(h / beta_h * (beta2 / mg)**(1/valence_Mg2), kex2) # 10**(3.4*(1-beta_h)) *  
    ! eq6 = sp.Eq(h / beta_h * (beta3 / na)**(1/valence_Na), kex3) # 10**(3.4*(1-beta_h)) * 
    ! eq7 = sp.Eq(h / beta_h * (beta4 / k)**(1/valence_K), kex4) # 10**(3.4*(1-beta_h)) * 
    ! eq8 = sp.Eq(h / beta_h * (beta5 / al)**(1/valence_Al3), kex5) # 10**(3.4*(1-beta_h)) * 
    ! eq9 = sp.Eq(h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al, b0)
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: soil_ph ! 
    real(r8), intent(in) :: b0 ! net charge balance (mol/kg)
    real(r8), intent(in) :: co2_atm ! atmospheric CO2 partial pressure (unit: atm)
    real(r8), intent(in) :: beta_list(1:ncations) ! fraction of cation exchange locations occupied by Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: kex_list(1:ncations)  ! exchange coefficient between H+ and Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: cation_valence(1:ncations)  ! valence of Ca2+, Mg2+, Na+, K+, Al3+
    real(r8) :: pcterr ! percentage error

    ! 
    ! !LOCAL VARIABLES:
    real(r8) :: h, beta_h, al_RHS, al_LHS

    !--------------------------------------------------------------

    h = 10**(-soil_ph)
    beta_h = 1.0_r8
    beta_h = beta_h - beta_list(1) - beta_list(2) - beta_list(3) - beta_list(4) - beta_list(5)

    al_RHS =  0.333333333333333*b0 &
            - 0.666666666666667*beta_list(1)/(beta_h*kex_list(1)/h)**cation_valence(1) &
            - 0.666666666666667*beta_list(2)/(beta_h*kex_list(2)/h)**cation_valence(2) &
            - 0.333333333333333*beta_list(3)/(beta_h*kex_list(3)/h)**cation_valence(3) &
            - 0.333333333333333*beta_list(4)/(beta_h*kex_list(4)/h)**cation_valence(4) &
            + (10**log_keq_hco3)/3._r8*co2_atm/h &
            + 0.666666666666667*(10**(log_keq_co3+log_keq_hco3))*co2_atm/h**2 &
            - 0.333333333333333*h &
            + 3.33333333333333e-15/h

    al_LHS = beta_list(5)/(beta_h*kex_list(5)/h)**(cation_valence(5))

    ! al_RHS = 0.333333333333333*b0 - 0.666666666666667*beta1/(0.000398107170553497*10.0**(3.4*beta_h)*beta_h*kex1/h)**valence['Ca2+'] - 0.666666666666667*beta2/(0.000398107170553497*10.0**(3.4*beta_h)*beta_h*kex2/h)**valence['Mg2+'] - 0.333333333333333*beta3/(0.000398107170553497*10.0**(3.4*beta_h)*beta_h*kex3/h)**valence['Na+'] - 0.333333333333333*beta4/(0.000398107170553497*10.0**(3.4*beta_h)*beta_h*kex4/h)**valence['K+'] + 5.12010356128343e-9*co2_atm/h + 4.80295746943524e-19*co2_atm/h**2 - 0.333333333333333*h + 3.33333333333333e-15/h
    ! al_LHS = beta2/(10.0**(3.4*beta_h - 3.4)*beta_h*kex2/h)**valence['Al3+']

    pcterr = abs(al_RHS - al_LHS) / (0.5*abs(al_RHS) + 0.5*abs(al_LHS))

  end function objective_solveq


  !-----------------------------------------------------------------------
  function solve_eq(b0, co2_atm, beta_list, kex_list, valence) result (best_ph)
    !
    ! !DESCRIPTION:
    ! Calculate whether a given pH value satisfies the following set of equations
    ! 
    ! eq1 = sp.Eq(h * hco3 / co2_atm, 10**(-7.8136))
    ! eq2 = sp.Eq(h * co3 / hco3, 10**(-10.3288))
    ! eq3 = sp.Eq(h * oh, 1e-14)
    ! eq4 = sp.Eq(h / beta_h * (beta1 / ca)**(1/valence_Ca2), kex1) # 10**(3.4*(1-beta_h)) *  
    ! eq5 = sp.Eq(h / beta_h * (beta2 / mg)**(1/valence_Mg2), kex2) # 10**(3.4*(1-beta_h)) *  
    ! eq6 = sp.Eq(h / beta_h * (beta3 / na)**(1/valence_Na), kex3) # 10**(3.4*(1-beta_h)) * 
    ! eq7 = sp.Eq(h / beta_h * (beta4 / k)**(1/valence_K), kex4) # 10**(3.4*(1-beta_h)) * 
    ! eq8 = sp.Eq(h / beta_h * (beta5 / al)**(1/valence_Al3), kex5) # 10**(3.4*(1-beta_h)) * 
    ! eq9 = sp.Eq(h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al, b0)
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: b0 ! net charge balance (mol/kg)
    real(r8), intent(in) :: co2_atm ! atmospheric CO2 partial pressure (unit: atm)
    real(r8), intent(in) :: beta_list(1:ncations) ! fraction of cation exchange locations occupied by Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: kex_list(1:ncations)  ! exchange coefficient between H+ and Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: valence(1:ncations)   ! valence of Ca2+, Mg2+, Na+, K+, Al3+
    real(r8) :: pcterr ! percentage error
    ! 
    ! !LOCAL VARIABLES:
    real(r8) :: best_ph, curr_ph, curr_err, min_err
    integer  :: i,j ! index
    integer  :: best_i
    integer  :: search_n
    real(r8) :: search_start, search_end, search_step

    ! Search the linear space to find where the pH minimizes error
    ! do four passes; fortran accuracy seems a little too low
    search_n = 501
    search_start = 0.5
    search_end = 13.5
    min_err = 999._r8
    j = 0
    do while ((j < 4) .and. (min_err > 0.001))
      search_step = (search_end - search_start) / (search_n - 1)
      do i = 1, search_n
        curr_ph = search_start + search_step * (i-1)
        curr_err = objective_solveq(curr_ph, b0, co2_atm, beta_list, kex_list, valence)
        if (curr_err < min_err) then
          best_i = i
          best_ph = curr_ph
          min_err = curr_err
        end if
      end do
      search_start = search_start + search_step * (max(best_i - 5, 1)-1)
      search_end = search_start + search_step * (min(best_i + 5, search_n)-1)
      j = j + 1
    end do

  end function solve_eq


end module ewutils
