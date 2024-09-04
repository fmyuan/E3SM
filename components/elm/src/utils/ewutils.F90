module ewutils
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing enhanced weathering shared utilities
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use elm_varcon  , only: log_keq_hco3, log_keq_co3
  use elm_varpar  , only: ncations
  use elm_varctl  , only: iulog
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
  public :: advection_diffusion

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


  subroutine advection_diffusion(conc_trcr,adv_flux,diffus,source,surf_bc,dtime,vwc,conc_change_rate)
    ! From B. Sulman; edited layer depth; soil bulk concentration can use g/m3
    ! 
    ! Advection and diffusion for a single tracer in one column given diffusion coefficient, flow, and source-sink terms
    ! Based on SoilLittVertTranspMod, which implements S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
    ! Not sure if this belongs here or somewhere else. Is it bad to do this in the EMI subroutine?

    use elm_varpar       , only : mixing_layer
    use elm_varcon       , only : zsoi, zisoi, dzsoi_decomp
    use abortutils       , only : endrun

    real(r8), intent(in) :: conc_trcr(1:mixing_layer)  ! Bulk concentration (e.g. mol/m3)
    real(r8), intent(in) :: adv_flux(1:mixing_layer+1) ! (m/s), vertical into layer (down is negative)
    real(r8), intent(in) :: diffus(1:mixing_layer)  ! diffusivity (m2/s)
    real(r8), intent(in) :: source(1:mixing_layer)  ! Source term (mol/m3/s)

    real(r8), intent(in) :: surf_bc                 ! Surface boundary layer concentration (for infiltration)
    real(r8), intent(in) :: dtime                   ! Time step (s)
    real(r8), intent(in) :: vwc(1:mixing_layer)     ! Volumetric soil moisture in layer (m3/m3)
    real(r8), intent(out):: conc_change_rate(1:mixing_layer) ! Bulk concentration (e.g. mol/m3/s)

    ! Local variables
    real(r8) :: aaa                          ! "A" function in Patankar
    real(r8) :: pe                           ! Pe for "A" function in Patankar
    real(r8) :: w_m1, w_p1                   ! Weights for calculating harmonic mean of diffusivity
    real(r8) :: d_m1, d_p1                   ! Harmonic mean of diffusivity
    real(r8) :: vwc_m1, vwc_p1                ! Harmonic mean of soil moisture
    real(r8) :: a_tri(0:mixing_layer+1)      ! "a" vector for tridiagonal matrix
    real(r8) :: b_tri(0:mixing_layer+1)      ! "b" vector for tridiagonal matrix
    real(r8) :: c_tri(0:mixing_layer+1)      ! "c" vector for tridiagonal matrix
    real(r8) :: r_tri(0:mixing_layer+1)      ! "r" vector for tridiagonal solution
    real(r8) :: d_p1_zp1(1:mixing_layer+1)   ! diffusivity/delta_z for next j  (set to zero for no diffusion)
    real(r8) :: d_m1_zm1(1:mixing_layer+1)   ! diffusivity/delta_z for previous j (set to zero for no diffusion)
    real(r8) :: f_p1(1:mixing_layer+1)       ! water flux for next j
    real(r8) :: f_m1(1:mixing_layer+1)       ! water flux for previous j
    real(r8) :: pe_p1(1:mixing_layer+1)      ! Peclet # for next j
    real(r8) :: pe_m1(1:mixing_layer+1)      ! Peclet # for previous j
    real(r8) :: dz_node(1:mixing_layer+1)    ! difference between nodes
    real(r8) :: a_p_0
    real(r8) :: conc_after(0:mixing_layer+1)

    integer :: j, info

    ! Statement function
    aaa (pe) = max (0._r8, (1._r8 - 0.1_r8 * abs(pe))**5)  ! "A" function from Patankar, Table 5.2, pg 95

    ! Set the distance between the node and the one ABOVE it   
    dz_node(1) = zsoi(1)
    do j = 2,mixing_layer+1
      dz_node(j)= zsoi(j) - zsoi(j-1)
    enddo

    write(iulog,*) 'adv_flux',adv_flux(1:mixing_layer+1)
    write(iulog,*) 'diffus',diffus(1:mixing_layer)
    write(iulog,*) 'source',source(1:mixing_layer)

    ! Calculate the D and F terms in the Patankar algorithm
    ! d: diffusivity
    ! f: flow
    ! m: layer above
    ! p: layer below
    ! pe: Peclet number (ratio of convection to diffusion)
    do j = 1,mixing_layer
      if (j == 1) then
        d_m1_zm1(j) = 0._r8
        w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
        if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
          d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
        else
          d_p1 = 0._r8
        endif
        d_p1_zp1(j) = d_p1 / dz_node(j+1)
        vwc_m1 = vwc(j)
        vwc_p1 = 1._r8 / ((1._r8 - w_p1) / vwc(j) + w_p1 / vwc(j+1))
        f_m1(j) = adv_flux(j) / vwc_m1  ! Include infiltration here
        f_p1(j) = adv_flux(j+1) / vwc_p1
        pe_m1(j) = 0._r8
        pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      elseif (j == mixing_layer) then
          ! At the bottom, assume no gradient in d_z (i.e., they're the same)
          w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
          if ( diffus(j) > 0._r8 .and. diffus(j-1) > 0._r8) then
            d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
          else
            d_m1 = 0._r8
          endif
          d_m1_zm1(j) = d_m1 / dz_node(j)
          d_p1_zp1(j) = d_m1_zm1(j) ! Set to be the same
          vwc_m1 = 1. / ((1. - w_m1) / vwc(j-1) + w_m1 / vwc(j))
          f_m1(j) = adv_flux(j) / vwc_m1
          !f_p1(j) = adv_flux(j+1)
          f_p1(j) = 0._r8
          pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
          pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      else
          ! Use distance from j-1 node to interface with j divided by distance between nodes
          w_m1 = (zisoi(j-1) - zsoi(j-1)) / dz_node(j)
          if ( diffus(j-1) > 0._r8 .and. diffus(j) > 0._r8) then
            d_m1 = 1._r8 / ((1._r8 - w_m1) / diffus(j) + w_m1 / diffus(j-1)) ! Harmonic mean of diffus
          else
            d_m1 = 0._r8
          endif
          w_p1 = (zsoi(j+1) - zisoi(j)) / dz_node(j+1)
          if ( diffus(j+1) > 0._r8 .and. diffus(j) > 0._r8) then
            d_p1 = 1._r8 / ((1._r8 - w_p1) / diffus(j) + w_p1 / diffus(j+1)) ! Harmonic mean of diffus
          else
            d_p1 = (1._r8 - w_p1) * diffus(j) + w_p1 * diffus(j+1) ! Arithmetic mean of diffus
          endif
          d_m1_zm1(j) = d_m1 / dz_node(j)
          d_p1_zp1(j) = d_p1 / dz_node(j+1)
          vwc_m1 = 1. / ((1. - w_m1) / vwc(j-1) + w_m1 / vwc(j))
          vwc_p1 = 1. / ((1. - w_p1) / vwc(j) + w_p1 / vwc(j+1))
          f_m1(j) = adv_flux(j) / vwc_m1
          f_p1(j) = adv_flux(j+1) / vwc_p1
          pe_m1(j) = f_m1(j) / d_m1_zm1(j) ! Peclet #
          pe_p1(j) = f_p1(j) / d_p1_zp1(j) ! Peclet #
      end if
    enddo ! j; mixing_layer


    ! Calculate the tridiagonal coefficients
    ! Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
    ! Here, this is equivalent to Patankar equation 5.56 and 5.57 (but in one dimension):
    ! a_P*phi_P = a_E*phi_E + a_W*phi_W + b [phi is concentration, = x in tridiagonal]. Converting East/West to above/below
    ! -> -a_E*phi_E + a_P*phi_P - a_W+phi_W = b
    ! -a_tri = a_above = D_above*A(Pe)+max(-F_above,0); D_above=diffus_above/dz
    ! b_tri = a_above+a_below+rho*dz/dt
    ! -c_tri = D_below*A(Pe)+max(F_below,0); D_below = diffus_below/dz
    ! r_tri = b = source_const*dz + conc*rho*dz/dt
    do j = 0,mixing_layer +1

      if (j > 0 .and. j < mixing_layer+1) then
          a_p_0 =  dzsoi_decomp(j) / dtime / vwc(j) ! Should this be multiplied by layer water content (for vwc)?
      endif

      if (j == 0) then ! top layer (atmosphere)
          a_tri(j) = 0._r8
          b_tri(j) = 1._r8
          c_tri(j) = -1._r8
          r_tri(j) = 0._r8
      elseif (j == 1) then
          a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
          c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
          b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
          ! r_tri includes infiltration assuming same concentration as top layer. May want to change to either provide upper boundary condition or include in source term
          ! r_tri(j) = source(j) * dzsoi_decomp(j) + (a_p_0 - adv_flux(j)) * conc_trcr(j)
          r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j)
          if(adv_flux(j)<0) then ! downward flow (infiltration)
            r_tri(j) = r_tri(j) - adv_flux(j)*surf_bc
            !  write (iulog,*) __LINE__,adv_flux(j),surf_bc,adv_flux(j)*surf_bc
          else ! upward flow to the surface
            r_tri(j) = r_tri(j) - adv_flux(j)*conc_trcr(j)
            ! write (iulog,*) __LINE__,adv_flux(j),conc_trcr(j),adv_flux(j)*conc_trcr(j)
          endif
          
      elseif (j < mixing_layer+1) then
          a_tri(j) = -(d_m1_zm1(j) * aaa(pe_m1(j)) + max( f_m1(j), 0._r8)) ! Eqn 5.47 Patankar
          c_tri(j) = -(d_p1_zp1(j) * aaa(pe_p1(j)) + max(-f_p1(j), 0._r8))
          b_tri(j) = -a_tri(j) - c_tri(j) + a_p_0
          r_tri(j) = source(j) * dzsoi_decomp(j) + a_p_0 * conc_trcr(j) ! Eq. 5.57
      else ! j==mixing_layer+1; 0 concentration gradient at bottom
          a_tri(j) = -1._r8
          b_tri(j) = 1._r8
          c_tri(j) = 0._r8 
          r_tri(j) = 0._r8
      endif
    enddo ! j; mixing_layer

    ! write(iulog,'(11a18)'),'a','b','c','r','ap0','pe_m','pe_p','f_m','f_p','d_m','d_p'
    ! j=0
    ! write(iulog,'(i3,4e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j)
    ! do j=1,mixing_layer
    !   write(iulog,'(i3,11e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j),dzsoi_decomp(j) / dtime * rho(j) ,pe_m1(j),pe_p1(j),f_m1(j),f_p1(j),d_m1_zm1(j)*dz_node(j),d_p1_zp1(j)*dz_node(j+1)
    ! enddo
    ! j=mixing_layer+1
    ! write(iulog,'(i3,4e18.9)'),j,a_tri(j),b_tri(j),c_tri(j),r_tri(j)

    ! Solve for the concentration profile for this time step
    ! call Tridiagonal(0, mixing_layer+1, 0, a_tri, b_tri, c_tri, r_tri, conc_after)
    ! This is the LAPACK tridiagonal solver which gave more accurate results in my testing
    call dgtsv( mixing_layer+2, 1, c_tri(0:mixing_layer), b_tri, a_tri(1:mixing_layer+1),  & 
                r_tri, mixing_layer+2, info )

    if(info < 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: illegal argument')
    if(info > 0) call endrun(msg='dgtsv error in adv_diff line __LINE__: singular matrix')
    conc_after = r_tri

    write (iulog,*) 'conc_before',conc_trcr
    write (iulog,*) 'conc_after',conc_after
    write (iulog,*) 'Diff=',sum((conc_after(1:mixing_layer)-conc_trcr)*dzsoi_decomp)
    write (iulog,*) 'Flow',adv_flux(1:mixing_layer+1)
    write (iulog,*) 'Diffus',diffus
    write (iulog,*) 'dz',dzsoi_decomp
    write (iulog,*) 'dznode',dz_node

    conc_change_rate = (conc_after(1:mixing_layer)-conc_trcr)/dtime

  end subroutine advection_diffusion

end module ewutils
