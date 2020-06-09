! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
!
! ***********************************************************************

      module atm_def
      
      use const_def, only: dp
      
      implicit none    

      integer, parameter :: atm_simple_photosphere = 1
      integer, parameter :: atm_Eddington_grey = 2 ! Eddington T-tau integration
      integer, parameter :: atm_Krishna_Swamy = 3 ! Krishna Swamy T-tau integration
      integer, parameter :: atm_solar_Hopf_grey = 4
         ! solar calibrated Hopf-function T-tau integration
         ! T^4 = 3/4 Teff^4 (tau + q(tau))
         ! q(tau) = q1 + q2 exp(-q3 tau) + q4 exp(-q5 tau)
         ! solar calibrated q's (from Jorgen Christensen-Dalsgaard) are
         !     q1 = 1.0361
         !     q2 = -0.3134 
         !     q3 = 2.44799995
         !     q4 = -0.29589999
         !     q5 = 30.0
         ! tau_photoshere is tau s.t. tau + q(tau) = 4/3 => tau_photosphere = 0.4116433502
         
      integer, parameter :: atm_tau_100_tables = 5
         ! use model atmosphere tables for Pgas and T at tau=100; solar Z only.
      integer, parameter :: atm_tau_10_tables = 6
         ! use model atmosphere tables for Pgas and T at tau=10; solar Z only.
      integer, parameter :: atm_tau_1_tables = 7 
         ! use model atmosphere tables for Pgas and T at tau=1; solar Z only.
      integer, parameter :: atm_tau_1m1_tables = 8
         ! use model atmosphere tables for Pgas and T at tau=1e-1; solar Z only.
                  
      integer, parameter :: atm_photosphere_tables = 9 
         ! use model atmosphere tables for photosphere Pgas; [Z/Z_SOLAR] from -4.0 to +0.5
      integer, parameter :: atm_grey_and_kap = 10 ! find consistent P, T, and kap at surface
      integer, parameter :: atm_grey_irradiated = 11  
         ! based on Guillot, T, and Havel, M., A&A 527, A20 (2011). see equation 6.
      integer, parameter :: atm_Paczynski_grey = 12
         ! integrate an atmosphere for given base conditions.
         ! inspired by B. Paczynski, 1969, Acta Astr., vol. 19, 1.
         ! takes into account dilution when tau < 2/3,
         ! and calls mlt to get gradT allowing for convection.
         ! note: only available from mesa/star since requires star lib information
      integer, parameter :: atm_WD_tau_25_tables = 13
         ! hydrogen atmosphere tables for cool white dwarfs
         ! giving Pgas and T at log10(tau) = 1.4 (tau = 25.11886)
         ! Teff goes from 40,000 K down to 2,000K with step of 100 K
         ! Log10(g) goes from 9.5 down to 5.5 with step of 0.1 
         ! reference
            ! R.D. Rohrmann, L.G. Althaus, and S.O. Kepler,
            ! Lyman α wing absorption in cool white dwarf stars,
            ! Mon. Not. R. Astron. Soc. 411, 781–791 (2011)
      integer, parameter :: atm_fixed_Teff = 14
         ! set Tsurf from Eddington T-tau relation for given Teff and tau
         ! set Psurf = Radiation_Pressure(Tsurf)
      integer, parameter :: atm_fixed_Tsurf = 15
         ! set Teff from Eddington T-tau relation for given Tsurf and tau=2/3
         ! set Psurf = Radiation_Pressure(Tsurf)
      integer, parameter :: atm_fixed_Psurf = 16
         ! set Tsurf from L and R using L = 4*pi*R^2*boltz_sigma*T^4.
         ! set Teff using Eddington T-tau relation for tau=2/3 and T=Tsurf.
      integer, parameter :: atm_fixed_Psurf_and_Tsurf = 17
         ! set Teff using Eddington T-tau relation for tau=2/3 and T=Tsurf.

      integer, parameter :: min_atm_option = 1 
      integer, parameter :: max_atm_option = 17
      
      
      ! info about structure of atmosphere
      integer, parameter :: atm_xm = 1 ! mass of atm exterior to this point (g)
      integer, parameter :: atm_delta_r = atm_xm+1 ! radial distance above base of envelope (cm)
      integer, parameter :: atm_lnP = atm_delta_r+1
      integer, parameter :: atm_lnd = atm_lnP+1
      integer, parameter :: atm_lnT = atm_lnd+1
      integer, parameter :: atm_gradT = atm_lnT+1
      integer, parameter :: atm_kap = atm_gradT+1
      integer, parameter :: atm_gamma1 = atm_kap+1
      integer, parameter :: atm_grada = atm_gamma1+1
      integer, parameter :: atm_chiT = atm_grada+1
      integer, parameter :: atm_chiRho = atm_chiT+1
      integer, parameter :: atm_cv = atm_chiRho+1
      integer, parameter :: atm_cp = atm_cv+1
      integer, parameter :: atm_lnfree_e = atm_cp+1
      integer, parameter :: atm_dlnkap_dlnT = atm_lnfree_e+1
      integer, parameter :: atm_dlnkap_dlnd = atm_dlnkap_dlnT+1
      integer, parameter :: atm_lnPgas = atm_dlnkap_dlnd+1
      integer, parameter :: atm_tau = atm_lnPgas+1
      integer, parameter :: atm_gradr = atm_tau+1

      integer, parameter :: num_results_for_create_atm = atm_gradr 
      
      
      ! tables
      integer, parameter :: table_atm_version = 6
      
      type Atm_Info
         integer :: which_atm_option, nZ, ng, nT, ilinT, iling
         real(dp), pointer :: Teff_array(:), logg_array(:), Teff_bound(:)
         real(dp), pointer :: logZ(:), alphaFe(:)
         real(dp), pointer :: Pgas_interp1(:), T_interp1(:)
         real(dp), pointer :: Pgas_interp(:,:,:,:), T_interp(:,:,:,:)
         character(len=8), pointer :: atm_mix(:)
         character(len=40), pointer :: table_atm_files(:)
         logical, pointer :: have_atm_table(:)
      end type Atm_Info
      
      type (Atm_Info), target :: &
         ai_two_thirds_info, ai_100_info, ai_10_info, ai_1_info, &
         ai_1m1_info, ai_wd_25_info
      type (Atm_Info), pointer :: &
         ai_two_thirds, ai_100, ai_10, ai_1, ai_1m1, ai_wd_25


      ! pairs of x and log10[ExpIntegralE[2,x]] from Mathematica
      integer, parameter :: npairs = 571
      real(dp), target :: E2_x(npairs)
      real(dp) :: E2_pairs(2*npairs)
      real(dp), target :: E2_f_ary(4*npairs)
      real(dp), pointer :: E2_f1(:), E2_f(:,:)
      logical :: have_E2_interpolant = .false.

      logical :: table_atm_is_initialized = .false.
      

      type Int_Atm_Info
         logical :: save_atm_structure_info
         integer :: atm_structure_num_pts
         real(dp), pointer :: atm_structure(:,:) ! will be allocated if necessary
            ! (num_results_for_create_atm, num_atm_points)
         ! bookkeeping
         integer :: handle
         logical :: in_use
      end type Int_Atm_Info
      
      logical :: int_atm_is_initialized = .false.
      integer, parameter :: max_atm_handles = 10
      type (Int_Atm_Info), target :: atm_handles(max_atm_handles)


      contains
      
      subroutine set_E2_pairs 
         include 'e2_pairs.dek'
      end subroutine set_E2_pairs

      end module atm_def

