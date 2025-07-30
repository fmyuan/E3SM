module ewutils
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing enhanced weathering shared utilities
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use elm_varcon  , only: log_keq_hco3, log_keq_co3
  use elm_varpar  , only: ncations
  use shr_sys_mod , only: shr_sys_flush
  use spmdMod     , only: iam
  use elm_varctl  , only : iulog
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
  public :: u_pdf
  public :: get_ssa
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: analytical_c
  private :: analytical_c_int
  private :: analytical_dt

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
    real(r8), intent(in) :: valence ! eq mol-1 (note this need x1000 to meq), valance of the cation or solid
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: bd ! kg m-3, soil bulk density
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
    real(r8), intent(in) :: log_mol_conc ! mol kg-1 waterfad
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
    real(r8), intent(in) :: valence ! eq mol-1 (note this need x1000 to meq), valance of the cation or solid
    real(r8), intent(in) :: molar_mass ! g mol-1, molar mass of the cation or solid
    real(r8), intent(in) :: bd ! kg m-3, soil bulk density
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
  function find_net_charge(soil_ph, co2_atm, beta_list, kex_list, cation_valence) result (charge)
    !
    ! !DESCRIPTION:
    ! Calculate the net charge in the system using the following set of equations
    !
    ! eq1 = sp.Eq(h * hco3 / co2_atm, 10**(-7.8136))
    ! eq2 = sp.Eq(h * co3 / hco3, 10**(-10.3288))
    ! eq3 = sp.Eq(h * oh, 1e-14)
    !
    ! eq4 = sp.Eq(h / beta_h * (beta1 / ca)**(1/valence_Ca2), kex1)
    ! eq5 = sp.Eq(h / beta_h * (beta2 / mg)**(1/valence_Mg2), kex2)
    ! eq6 = sp.Eq(h / beta_h * (beta3 / na)**(1/valence_Na), kex3)
    ! eq7 = sp.Eq(h / beta_h * (beta4 / k)**(1/valence_K), kex4)
    ! eq8 = sp.Eq(h / beta_h * (beta5 / al)**(1/valence_Al3), kex5)
    !
    ! Aluminum hydrolysis will not be considered because we did not have those extra
    ! Al-species's cation exchange coefficients.
    !! eq9 = sp.Eq(aloh * h / al, 10**(-5))
    !! eq10 = sp.Eq(aloh2 * h * h / al, 10**(-10.1))
    !! eq11 = sp.Eq(aloh3 * h * h * h / al, 10**(-16.9))
    !! eq12 = sp.Eq(aloh4 * h * h * h * h / al, 10**(-22.7))
    !
    !! eq13 = sp.Eq(h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al + 2*aloh + aloh2 - aloh4, b0)
    ! eq13 = sp.Eq(h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al, b0)
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: soil_ph ! 
    real(r8), intent(in) :: co2_atm ! atmospheric CO2 partial pressure (unit: atm)
    real(r8), intent(in) :: beta_list(1:ncations) ! fraction of cation exchange locations occupied by Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: kex_list(1:ncations)  ! exchange coefficient between H+ and Ca2+, Mg2+, Na+, K+, Al3+
    real(r8), intent(in) :: cation_valence(1:ncations)  ! valence of Ca2+, Mg2+, Na+, K+, Al3+
    real(r8) :: pcterr ! percentage error

    ! 
    ! !LOCAL VARIABLES:
    real(r8) :: h, beta_h, charge, oh, hco3, co3, ca, mg, na, k, al !!, aloh, aloh2, aloh3, aloh4

    !--------------------------------------------------------------

    h = 10**(-soil_ph)
    beta_h = 1.0_r8 - beta_list(1) - beta_list(2) - beta_list(3) - beta_list(4) - beta_list(5)

    oh = 10**(soil_ph-14_r8)
    hco3 = 1.53603106838503_r8*10**(-8_r8 + soil_ph)*co2_atm
    co3 = 7.20443620415286_r8*10**(-19_r8 + 2_r8*soil_ph)*co2_atm
    ca = beta_list(1)/(beta_h*kex_list(1)/h)**cation_valence(1)
    mg = beta_list(2)/(beta_h*kex_list(2)/h)**cation_valence(2)
    na = beta_list(3)/(beta_h*kex_list(3)/h)**cation_valence(3)
    k = beta_list(4)/(beta_h*kex_list(4)/h)**cation_valence(4)
    al = beta_list(5)/(beta_h*kex_list(5)/h)**cation_valence(5)
    !aloh = 10**(soil_ph-5)*al
    !aloh2 = 10**(2_r8*soil_ph-10.1_r8)*al
    !aloh3 = 10**(3_r8*soil_ph-16.9_r8)*al
    !aloh4 = 10**(4_r8*soil_ph-22.7_r8)*al

    charge = h - oh - hco3 - 2*co3 + 2*ca + 2*mg + na + k + 3*al ! + 2*aloh + aloh2 - aloh4

  end function find_net_charge


  function objective_solveq(soil_ph, b0, co2_atm, beta_list, kex_list, cation_valence) result (pcterr)
    !
    ! !DESCRIPTION:
    ! Calculate whether a given pH value satisfies the listed set of equations in find_net_charge
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
    real(r8) :: charge

    !--------------------------------------------------------------

    charge = find_net_charge(soil_ph, co2_atm, beta_list, kex_list, cation_valence)
    pcterr = abs(charge - b0) / abs(b0)

  end function objective_solveq


  !-----------------------------------------------------------------------
  function solve_eq(b0, co2_atm, beta_list, kex_list, valence) result (best_ph)
    !
    ! !DESCRIPTION:
    ! Calculate whether a given pH value satisfies the set of equations provided in
    ! objective_solveq
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
    search_n = 161
    search_start = 2 ! only acid mines reach this level
    search_end = 10 ! car wash level
    min_err = 1e10
    j = 0
    do while ((j < 8) .and. (min_err > 0.01))
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
      ! (make sure the upper bound is not beyond previous round)
      search_end = min(search_end, search_start + search_step * (min(best_i + 5, search_n)-1))
      j = j + 1
    end do

  end function solve_eq


  pure function analytical_c(C0, r, k, t1, t2) result(c)
    !------------------------------------------------------------------
    ! Concentration at t2, starting from C0 at t1,  dC/dt = k C + r
    !------------------------------------------------------------------
    real(r8), intent(in) :: C0, r, k, t1, t2
    real(r8)             :: c
    c = (C0 + r/k) * exp(k*(t2 - t1)) - r/k
  end function analytical_c


  pure function analytical_c_int(C0, r, k, t1, t2) result(c_int)
    !------------------------------------------------------------------
    ! Integral of C over [t1, t2] for the same ODE
    !------------------------------------------------------------------
    real(r8), intent(in) :: C0, r, k, t1, t2
    real(r8)             :: c_int
    c_int = (C0 + r/k)/k * (exp(k*(t2 - t1)) - 1.0_r8) - (r/k) * (t2 - t1)
  end function analytical_c_int


  pure function analytical_dt(c1, c2, r, k) result(dt)
    !------------------------------------------------------------------
    ! Time Δt needed for C to evolve from c1 to c2
    !------------------------------------------------------------------
    real(r8), intent(in) :: c1, c2, r, k
    real(r8)             :: dt
    dt = log((c2 + r/k) / (c1 + r/k)) / k
  end function analytical_dt


  subroutine advection_diffusion(conc_in, rain_conc, q_int, theta, watsat, srcsk, &
                                 d0, dt, dz, nlevbed, dcdt)
    !------------------------------------------------------------------
    ! Use the explicit solution at each time step
    ! 
    ! If we ignore spatial dependency, the equation becomes very simple
    ! 
    ! Consider a generic cell C_i with adjacent cells
    !     C_{i-1} <--- Δ x_i ---> C_i, D_{eff,i} <--- Δ x_{i+1} ---> C_{i+1}
    ! 
    ! Total mass change within C_i due to outflow
    !     Δ x  * θ * \frac{ dC_i }{ dt }
    ! 
    ! Diffusion to the upper cell, if C_i > C_{i-1}
    !     dup = D_{eff,i} * I[C_i > C_{i-1}] * \frac{C_i - C_{i-1}}{Δ x_i}
    ! 
    ! Diffusion to the lower cell, if C_i > C_{i+1}
    !     dlow = D_{eff,i} * I[C_i > C_{i+1}] * \frac{C_i - C_{i+1}}{Δ x_{i+1}}
    ! 
    ! Advection to upper cell, if q_{in,i} < 0
    !     aup = I[q_{in,i} < 0] * abs( q_{in,i} ) * C_i
    ! 
    ! Advection to lower cell, if q_{out,i} > 0
    !     alow = I[q_{out,i} > 0] * q_{out,i} * C_i
    ! 
    ! Rain inflow
    !     rain = I[q_{in,i} > 0] * q_{in,0} * C_{rain}
    ! 
    ! The relationship between self-change and the four outflow and internal source/sink, 
    !     after simplification, is
    ! 
    !     Δ x * θ * \frac{ dC_i }{ dt } = 
    !       - ( \frac{ D_{eff,i} * I[C_i > C_{i-1}] }{Δ x_i} + 
    !           \frac{ D_{eff,i} * I[C_i > C_{i+1}] }{Δ x_{i+1}} + 
    !           I[q_{in,i} < 0] * abs( q_{in,i} ) + I[q_{out,i} > 0] * q_{out,i} ) * C_i
    !       + ( D_{eff,i} * I[C_i > C_{i-1}] \frac{C_{i-1}}{Δ x_i} )
    !       + ( D_{eff,i} * I[C_i > C_{i+1}] \frac{C_{i+1}}{Δ x_{i+1}} )
    !       + I[q_{in,i} > 0] * q_{in,0} * C_{rain}
    !       + R
    ! 
    !     , which is of the form dC/dt = kC + r, and the analytical solution is just
    !          C = (C0 + r/k)e^{kt} - r/k
    !          on the interval [t, t + Δ t]
    ! 
    ! The equations means when t → ∞, C → -r/k
    ! 
    ! It is guaranteed that k < 0. Therefore, 
    ! 
    ! if r > 0, the concentration will not be negative for large t (we can be sure )
    ! 
    ! If r < 0, the concentration may be negative because the equilibrium value is
    !     negative. This can only happen when r < 0. In this case, the analytical 
    !     form may need to be assessed piecewise to insure non-negativity:
    ! 
    !     Step 1. We evaluate the analytical form at (t + Δ t),
    !             if the resulting C_t' > C_{i-1} and C_t' > C_{i+1}, then C_t = C_t', 
    !             otherwise, we need to go to step 2
    !     Step 2. We assess piecewise: 
    !             first, find the dt' such that C(t + dt') = max(C_{i-1}, C_{i+1})
    !             then, on the interval [t+dt', t + Δ t], set the appropriate I[⋅]
    !                    term to 0, re-evaluate to t + Δ t
    !             third, if the result gives C(t + Δ t) > min(C_{i-1}, C_{i}), stop
    !                    otherwise, need a third piecewise, go to step 3
    !     Step 3. find the dt" such that C(t + dt") = min(C_{i-1}, C_{i+1}), on the
    !             interval [t+dt", t + Δ t], all the diffusion terms need to be 0.
    !             Re-integrate the pure-advection term to t + Δ t
    !------------------------------------------------------------------
    use elm_varcon       , only : zsoi, zisoi
    use elm_varpar       , only : nlevsoi
    use abortutils       , only : endrun

    ! Input-output variables
    real(r8), intent(in) :: conc_in(1:nlevsoi) ! start concentration (mol/m3-soil)
    real(r8), intent(in) :: rain_conc ! upper boundary condition (rain chemistry) (mol/m3-water)
    real(r8), intent(in) :: q_int(1:nlevsoi+1) ! water flux at grid boundaries, positive downwards (m/s)
    real(r8), intent(in) :: theta(1:nlevsoi) ! soil moisture (m3/m3)
    real(r8), intent(in) :: watsat(1:nlevsoi) ! soil porosity (m3/m3)
    real(r8), intent(in) :: srcsk(1:nlevsoi) ! source/sink strength (mol/m3-soil/s)
    real(r8), intent(in) :: d0 ! diffusion coefficient in water (m2/s)
    real(r8), intent(in) :: dt ! time step size (s)
    real(r8), intent(in) :: dz(1:nlevsoi) ! soil layer thickness (m)
    integer , intent(in) :: nlevbed ! number of hydrologically active layers
    real(r8), intent(out) :: dcdt(1:nlevsoi) ! rate of change of concentration (mol/m3-soil/s)

    ! Local variables
    real(r8) :: c_prev(1:nlevsoi) ! start concentration, converted to (mol/m3-water)
    real(r8) :: c_next(1:nlevsoi) ! end concentration (mol/m3-water)
    real(r8) :: Deff(1:nlevsoi) ! effective diffusion coefficient in porous media (m2/s)
    real(r8) :: dx(1:nlevsoi+2) ! padded dz; the first & last elements don't matter so long as >0
    real(r8) :: scaler(1:nlevsoi) ! \Delta x * 0 on the LHS
    real(r8) :: sourcesink(1:nlevsoi) ! source sink, converted from mol/m3-soil/s to mol/m2/s; with upper boundary
    integer  :: i, i1, i2, i3, i4 ! indices and boolean indicators
    integer  :: i1_keep, i2_keep  ! side of diffusion
    real(r8) :: k, r ! temporary variables to save the ODE parameters
    real(r8) :: c, c_p, cmax, cmin, c_int, c_int_p, c_int_pp, dt_p, dt_pp ! temporary variables for ODE integration

    ! Local variables, can be converted to diagnostic outputs if needed
    real(r8) :: niter(1:nlevsoi) ! keep track how many piecewise integrations are done
    real(r8) :: dc_up(1:nlevsoi) ! flow upward due to diffusion and advection (mol/m2-ground/s)
    real(r8) :: dc_down(1:nlevsoi) ! flow downward due to diffusion and advection (mol/m2-ground/s)

    ! the algorithm generates negative values when sourcesink gets too negative
    ! because my sourcesink is guaranteeded not to directly cause negative
    ! concentrations, the best approach is to update the cation concentrations
    ! first, and then apply vertial movement

    ! convert concentration from per m3-soil to per m3-water
    c_prev = (conc_in + srcsk*dt) / theta

    ! convert diffusion coefficient from in soil to in water
    Deff = d0 * theta**(10.0_r8/3.0_r8) / watsat**2

    ! pad the dz's
    dx(2:nlevsoi+1) = dz(1:nlevsoi)
    dx(1) = 100._r8
    dx(nlevsoi+2) = 100._r8

    ! LHS scaler
    scaler = dz * theta

    sourcesink(1:nlevsoi) = 0._r8
    if (q_int(1) > 0._r8) then
      ! convert rainfall input from mol/m3-water/s to mol/m2/s
      sourcesink(1) = sourcesink(1) + q_int(1) * rain_conc
    end if

    do i = 1,nlevbed
      if (i == 1) then
        ! top layer: no diffusion to above
        i1 = 0
        i2 = merge(1, 0, c_prev(i) > c_prev(i+1)) ! boolean to integer

      elseif (i < nlevbed) then
        ! middle layers
        i1 = merge(1, 0, c_prev(i) > c_prev(i-1))
        i2 = merge(1, 0, c_prev(i) > c_prev(i+1))

      else
        ! bottom layer: no diffusion to below
        i1 = merge(1, 0, c_prev(i) > c_prev(i-1))
        i2 = 0
      end if

      i3 = merge(1, 0, q_int(i) < 0)
      i4 = merge(1, 0, q_int(i+1) > 0)

      k = - (i1*Deff(i)/dx(i) + i2*Deff(i)/dx(i+1) + i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)

      if (i == 1) then
        r = (Deff(i)*i2*c_prev(i+1)/dx(i+1) + sourcesink(i)) / scaler(i)
      elseif (i < nlevbed) then
        r = (Deff(i)*i1*c_prev(i-1)/dx(i) + Deff(i)*i2*c_prev(i+1)/dx(i+1) \
                + sourcesink(i)) / scaler(i)
      else
        r = (Deff(i)*i1*c_prev(i-1)/dx(i) + sourcesink(i)) / scaler(i)
      end if

      !!write (iulog, *) 'advection_diffusion', i, i1, i2, i3, i4, k, r

      ! if there is no flow out of this cell, degrades to linear source
      if (k == 0._r8) then
        c_next(i) = c_prev(i) + r*dt
        dc_up(i) = 0
        dc_down(i) = 0
        niter(i) = 0

      else
        ! otherwise, do actual calculations

        ! step 1: find end-of-step solution
        c = analytical_c(c_prev(i), r, k, 0._r8, dt)

        ! top layer
        if (i == 1) then
          cmax = c_prev(i+1)

          ! end-of-step solution works, or there is no diffusion to begin with
          if ((c > cmax) .or. (i2 == 0)) then
              c_next(i) = c
              niter(i) = 1

              ! integration of c over [0, dt]
              c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

              ! use c_int to calculate the fluxes
              dc_up(i) = i3*abs(q_int(i))*c_int
              dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int - c_prev(i+1)*dt) + i4*q_int(i+1)*c_int

          else
              ! step 2: piecewise integrate to cmax, then update k, r and continue
              dt_p = analytical_dt(c_prev(i), cmax, r, k)

              ! integration of c over [0, dt']
              c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

              ! use c_int to calculate the fluxes during [0, dt']
              dc_up(i) = i3*abs(q_int(i))*c_int_p
              dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int_p - c_prev(i+1)*dt_p) &
                            + i4*q_int(i+1)*c_int_p

              ! update k & r (no more diffusion)
              k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
              r = sourcesink(i) / scaler(i)

              ! itegration of c over [dt', dt]
              c_int_pp = analytical_c_int(cmax, r, k, dt_p, dt)

              ! add to the fluxes
              dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
              dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

              ! find final value
              c_next(i) = analytical_c(cmax, r, k, dt_p, dt)
              niter(i) = 2

          end if

        ! intermediate layer
        else if (i < nlevbed) then

          if ((i1 == 0) .and. (i2 == 0)) then
            ! if there is no diffusion to begin with, end of solution works

            ! final value found
            c_next(i) = c
            niter(i) = 1

            ! integration of c over (0, dt)
            c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

            ! calculate the fluxes during (0, dt) (all diffusion)
            dc_up(i) = i3*abs(q_int(i))*c_int
            dc_down(i) = i4*q_int(i+1)*c_int

          else

            ! if there is only one side diffusion
            if (i1 == 0) then
                cmax = c_prev(i+1) ! diffuse down
                cmin = 0._r8
                i1_keep = 0 ! once cmax is hit, no more diffusion
                i2_keep = 0

            else if (i2 == 0) then
                cmax = c_prev(i-1) ! diffuse up
                cmin = 0._r8
                i1_keep = 0 ! once cmax is hit, no more diffusion
                i2_keep = 0

            else
                ! both sides have diffusion, decide which side to begin with

                if (c_prev(i-1) > c_prev(i+1)) then
                    cmax = c_prev(i-1)
                    cmin = c_prev(i+1)
                    i1_keep = 0
                    i2_keep = 1 ! once cmax is hit, cmin still diffuse
                else
                    cmax = c_prev(i+1)
                    cmin = c_prev(i-1)
                    i1_keep = 1 ! once cmax is hit, cmin still diffuse
                    i2_keep = 0
                end if

            end if

            ! if end of solution works
            if (c > cmax) then

                ! final value found
                c_next(i) = c
                niter(i) = 1

                ! integration of c over [0, dt]
                c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

                ! calculate the fluxes during [0, dt] (all diffusion)
                dc_up(i) = i1*Deff(i)/dx(i) * (c_int - c_prev(i-1)*dt) + i3*abs(q_int(i))*c_int
                dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int - c_prev(i+1)*dt) + i4*q_int(i+1)*c_int

            else
              ! step 2: piecewise integrate to cmax, then update k, r and continue
              dt_p = analytical_dt(c_prev(i), cmax, r, k)

              ! integration of c over (0, dt')
              c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

              ! calculate the fluxes during (0, dt') (all diffusion)
              dc_up(i) = i1*Deff(i)/dx(i) * (c_int_p - c_prev(i-1)*dt_p) &
                          + i3*abs(q_int(i))*c_int_p
              dc_down(i) = i2*Deff(i)/dx(i+1) * (c_int_p - c_prev(i+1)*dt_p) &
                            + i4*q_int(i+1)*c_int_p

              !DEBUG
              !write (iulog, *) '1 dc_up', dc_up(i), i, i1, Deff(i), dx(i), c_int_p, c_prev(i-1), dt_p, i3, abs(q_int(i))
              !write (iulog, *) '1 dc_down', dc_down(i), i, i2, Deff(i), dx(i+1), c_int_p, c_prev(i+1), dt_p, i4, q_int(i+1)

              ! update k & r (drop one side of diffusion)
              k = - (i1_keep*i1*Deff(i)/dx(i) + i2_keep*i2*Deff(i)/dx(i+1) + &
                     i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
              r = (Deff(i)*i1_keep*i1*c_prev(i-1)/dx(i) + &
                    Deff(i)*i2_keep*i2*c_prev(i+1)/dx(i+1) + sourcesink(i)) / scaler(i)

              ! continue to end, but we need another check
              c_p = analytical_c(cmax, r, k, dt_p, dt)

              if (c_p > cmin) then

                ! itegration of c over [dt', dt]
                c_int_p = analytical_c_int(cmax, r, k, dt_p, dt)

                ! add the fluxes during (dt', dt) (drop one side of diffusion)
                dc_up(i) = dc_up(i) &
                    + i1_keep*i1*Deff(i)/dx(i) * (c_int_p - c_prev(i-1)*(dt-dt_p)) &
                    + i3*abs(q_int(i))*c_int_p
                dc_down(i) = dc_down(i) &
                    + i2_keep*i2*Deff(i)/dx(i+1) * (c_int_p - c_prev(i+1)*(dt-dt_p)) &
                    + i4*q_int(i+1)*c_int_p

                !DEBUG
                !write (iulog, *) '2 dc_up', dc_up(i), i, i1_keep, i1, Deff(i), dx(i), c_int_p, cmax, dt_p, dt, i3, abs(q_int(i))
                !write (iulog, *) '2 dc_down', dc_down(i), i, i2_keep, i2, Deff(i), dx(i+1), c_int_p, cmax, dt_p, dt, i4, q_int(i+1)

                ! final value at dt
                c_next(i) = c_p
                niter(i) = 2

              else
                ! step 3: do another piecewise integration

                ! second stop point
                dt_pp = analytical_dt(cmax, cmin, r, k)

                ! itegration of c over [dt', dt"]
                c_int_pp = analytical_c_int(c_p, r, k, dt_p, dt_pp)

                ! add the fluxes during [dt', dt"] (drop one side of diffusion)
                dc_up(i) = dc_up(i) &
                    + i1_keep*i1*Deff(i)/dx(i) * (c_int_pp - c_prev(i-1)*(dt_pp-dt_p)) &
                    + i3*abs(q_int(i))*c_int_pp
                dc_down(i) = dc_down(i) &
                    + i2_keep*i2*Deff(i)/dx(i+1) * (c_int_pp - c_prev(i+1)*(dt_pp-dt_p)) &
                    + i4*q_int(i+1)*c_int_pp

                ! update the k & r (drop all diffusion)
                k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
                r = sourcesink(i) / scaler(i)

                ! itegration of c over [dt", dt]
                c_int_pp = analytical_c_int(cmin, r, k, dt_pp, dt)

                ! add the fluxes during [dt", dt]
                dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
                dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

                ! final value
                c_next(i) = analytical_c(cmin, r, k, dt_pp, dt)
                niter(i) = 3
              end if
            end if
          end if

        ! last soil layer
        else
          cmax = c_prev(i-1)

          ! end-of-step solution works, or there is no diffusion to begin with 
          if ((c > cmax) .or. (i1 == 0)) then
            c_next(i) = c
            niter(i) = 1

            ! integration of c over (0, dt)
            c_int = analytical_c_int(c_prev(i), r, k, 0._r8, dt)

            ! use c_int to calculate the fluxes
            dc_up(i) = i1*Deff(i)/dx(i) * (c_int - c_prev(i-1)*dt) + i3*abs(q_int(i))*c_int
            dc_down(i) = i4*q_int(i+1)*c_int

          else
            ! step 2: piecewise integrate to cmax, then update k, r and continue
            dt_p = analytical_dt(c_prev(i-1), cmax, r, k)

            ! integration of c over [0, dt']
            c_int_p = analytical_c_int(c_prev(i), r, k, 0._r8, dt_p)

            ! calculate the fluxes during [0, dt'] (all diffusion)
            dc_up(i) = i1*Deff(i)/dx(i) * (c_int_p - c_prev(i-1)*dt_p) + i3*abs(q_int(i))*c_int_p
            dc_down(i) = i4*q_int(i+1)*c_int_p

            ! update k & r (no more diffusion)
            k = - (i3*abs(q_int(i)) + i4*q_int(i+1)) / scaler(i)
            r = sourcesink(i) / scaler(i)

            ! if non-negative, proceed to end of time step

            ! itegration of c over [dt', dt]
            c_int_pp = analytical_c_int(cmax, r, k, dt_p, dt)

            ! add to the fluxes
            dc_up(i) = dc_up(i) + i3*abs(q_int(i))*c_int_pp
            dc_down(i) = dc_down(i) + i4*q_int(i+1)*c_int_pp

            ! find final value
            c_next(i) = analytical_c(cmax, r, k, dt_p, dt)
            niter(i) = 2

          end if

        end if ! end distinction of first, last, and intermediate layers

      end if ! end distinction between having diffusion-advection or not

    end do ! end iteration thru soil layers

    ! we still need to catch the case when r is simply too negative
    ! in that case, we really need to reduce r (secondary mineral
    ! precipitation, etc.)
    ! (TBD)

    !write (iam+100, *) 'c_prev', c_prev(1:nlevbed)
    !write (iam+100, *) 'q_int', q_int(1:nlevbed+1)
    !write (iam+100, *) 'dc_up', dc_up(1:nlevbed)
    !write (iam+100, *) 'dc_down', dc_down(1:nlevbed)
    !write (iam+100, *) 'c_next', c_next(1:nlevbed)
    !write (iam+100, *) 'niter', niter(1:nlevbed)

    ! calculate the net between self-outflow and inflow fluxes
    ! note the inflow fluxes need to be scalerd by soil moisture
    ! to get the correct concentration implications
    !c_next(1:nlevbed-1) = c_next(1:nlevbed-1) + dc_up(2:nlevbed)/theta(1:nlevbed-1)/dz(1:nlevbed-1)
    !c_next(2:nlevbed) = c_next(2:nlevbed) + dc_down(:nlevbed-1)/theta(2:nlevbed)/dz(2:nlevbed)

    !write (iam+100, *) 'c_next2', c_next(1:nlevbed)

    ! convert end concentration from per m3-water to per m3-soil, and perform delta
    dcdt = (c_next - c_prev) * theta / dt

  end subroutine advection_diffusion


  function u_pdf(u, theta) result(f)
    ! Gamma PDF with k = 0.5 transformed by x = u**2
    !   pdf(x)dx = pdf(u^2) * 2u du
    !    note gamma(0.5) = np.sqrt(pi)
    real(r8), intent(in) :: u, theta
    real(r8) :: f
    f = 2 / sqrt(theta * 3.1415926535_r8) * exp(-u**2 / theta)
  end function u_pdf


  function get_ssa(gs) result(a)
    ! calculate the specific surface area given grain size in um
    real(r8), intent(in) :: gs
    real(r8) :: a
    a = 69.18_r8 * (gs ** (-1.24_r8)) ! unit: m^2 g-1
  end function get_ssa

  function get_r(r0, m0, m) result(r)
    ! calculate the change in grain size after dissolution
    ! per time step
    !
    ! assume general density of rock rho = 2.9 g/cm3
    ! note: 1 meter = 1e6 um
    ! 
    ! V - volume [m3], m - mass [g], r - radius [m]
    !
    !! dV = SSA * m * dr
    !! => dm = SSA * rho * m * dr
    !! => dm = 69.18 * (2 * r * 1e6) ^{-1.24} * rho * m * dr
    !! => ln(m/m0) = (69.18 * 2 * 1e6^{-1.24} * rho / (-0.24)) * (r^{-1.24} - r0^{-1.24})
    !! => ln(m/m0) = -60.70 (r^{-1.24} - r0^{-1.24})
    !! => r = [ - 0.0164744 * ln(m/m0) + r0^{-1.24} ]^{-1/1.24}
    !
    ! dV = 4 pi r^2 * dr
    ! => dr = (36*pi)^{-1/3} V^{-2/3} dV
    ! => r - r0 = (36*pi)^{-1/3} * 3 * (V^{1/3} - V0^{1/3})

    real(r8), intent(in) :: r0 ! start grain size (um)
    real(r8), intent(in) :: m0 ! start mineral mass (g)
    real(r8), intent(in) :: m  ! end mineral mass (g)
    real(r8) :: r ! end mineral grain size
    real(r8) :: rho = 2.9e6

    ! r = (-0.03295_r8 * log(m/m0) + r0**(-1.24_r8))**(-1/1.24_r8)
    r = r0 + (36*3.1415926)**(-1/3) * 3 * ((m/rho)**(1/3_r8) - (m0/rho)**(1/3_r8))
  end function get_r

end module ewutils