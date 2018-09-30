! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras
      
      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      use rates_def
      use net_def
      
      implicit none
      
      real(dp) :: original_diffusion_dt_limit
      real(dp) :: burn_check = 0.0
      logical :: debug_use_other_torque = .false.
      logical :: debug_get_cz_info = .false.

      type conv_zone_info
       real(dp) :: &
         top_radius, &
         bot_radius, &
         d_radius, &
         top_mass, &
         bot_mass, &
         d_mass
       integer :: &
         top_zone, &
         bot_zone, &
         d_zone
      end type conv_zone_info

      type (conv_zone_info), dimension(:), allocatable :: conv_zone_list

      
!     these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++ROQUE++++++++++++++++++++++++++++++++++++++++++++++++++++'
         
         original_diffusion_dt_limit = s% diffusion_dt_limit
         !s% other_wind => Reimers_then_VW
         s% other_wind => Reimers_then_Blocker

         s% other_torque => tfm_other_torque

         !debug flags
         debug_use_other_torque = s% x_logical_ctrl(1)
         debug_get_cz_info = s% x_logical_ctrl(2)

      
      ! Once you have set the function pointers you want,
      ! then uncomment this (or set it in your star_job inlist)
      ! to disable the printed warning message,
         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  


      s% job% warn_run_star_extras =.false.             
      end subroutine extras_controls

      subroutine tfm_other_torque(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: r_st, m_st, i_st
         real(dp) :: j_dot, omega_surf, m_dot, eta_surf, v_inf, v_esc, B
         type (conv_zone_info), target :: cz_info
         type (conv_zone_info), pointer :: cz_info_ptr
         real(dp), dimension(:), pointer :: mag_brk_jdot;
         
         cz_info_ptr => cz_info

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         allocate(mag_brk_jdot(s% nz))
         s% extra_jdot(:) = 0
         !s% extra_omegadot(:) = 0

         !Magnetic braking according to MESA school 2012 assignment by Cantiello
         !j_dot = 2/3*m_dot*omega*alfven_r*alfven_r
         ! j_dot = angular momentum lost
         ! omega_surf = surface angular velocity (rad/s)
         ! m_dot = mass lost rate (Msun/year)
         ! alfven_r = Alfven radius

         !Wind-confinmenet parameter eta_surf
         !eta_surf = ((r_st*r_st)/(B*B)) / (m_dot*v_inf)
         ! r_st = stellar surface radius (cm)
         ! B = magnetic field torque (G)
         ! v_inf = terminal velocity of the stellar wind (cm/s)

         !v_inf = 1.92*v_esc
         ! v_esc = photospheric escape velocity (cm/s)

         !v_esc = 618*((r_sol/r_st)*(m_st/m_sol))^(1/2)

         !j_dot = 2/3*m_dot*omega_surf*r_st*r_st*eta_surf

         if ((s% use_other_torque) .and. (s% star_mdot < 0.0)) then

            call get_convective_info(s, cz_info_ptr)

            !Star data
            r_st = s% r(1)
            m_st = s% m(1)
            !r_st = Rsun
            !m_st = Msun
            i_st = s% i_rot(1)
            !omega_surf = s% omega_avg_surf
            omega_surf = s% omega(1)

            v_esc = (618 * ((Rsun/r_st)*(m_st/Msun))**0.5) * 100000 !100000 transform from km/s to cm/s

            v_inf = 1.92 * v_esc

            B = s% x_ctrl(6)

            !m_dot = s% star_mdot
            m_dot = s% mstar_dot


            !eta = ( s% photosphere_r * rsun * bfield )**2.0 / (abs( s% mstar_dot ) * vinf) 
            !with vinf in cm/s (the wind terminal velocity) rsun in cm, bfield in Gauss 
            !rest are MESA vars accessible through the star structure in run_star_extras.f
            !eta_surf = abs(((r_st/Rsun)**2/B**2)/(m_dot * v_inf)) !100000 pass from km/s to cm/s
            eta_surf = ((r_st * B)**2)/(abs(m_dot) * v_inf)
                        

            !j_dot = two_thirds * m_dot * omega_surf * (r_st/Rsun)**2 * eta_surf
            j_dot = two_thirds * m_dot * omega_surf * (r_st**2) * eta_surf

            !call distribute_j_dot(s, j_dot, cz_info_ptr, mag_brk_jdot)
            call distribute_simple_j_dot(s, j_dot, cz_info_ptr, mag_brk_jdot)
            !It happens that s% extra_jdot is longer than s% nz but mag_brk_jdot is just defined
            !for s% nz elements
            s% extra_jdot(1:s% nz) = mag_brk_jdot

            if (debug_use_other_torque) then
               write(*,*) "Rsun=", Rsun, "Msun=", Msun, "r_st=", r_st, "m_st=", m_st, &
                  "v_esc=", v_esc, "v_inf", v_inf, "B", B, "m_dot", m_dot, "eta_surf", eta_surf, &
                  "omega_surf", omega_surf, "j_dot", j_dot, "i_st=", i_st
                  do k=1, size(mag_brk_jdot)
                  !      k=30 !sacamos solo el valor de la zona 30 para acortar la salida
                        write(*,*) "jdot(k)=", s% extra_jdot(k), "omega(k)*i_rot(k)=", s% omega(k) * s% i_rot(k), &
                        "omega(k)=", s% omega(k), "i_rot(k)=", s% i_rot(k)
                  end do
                  !s% extra_jdot(:) = 0.0
                  !write(*,*) "Size s% extra_jdot = ", size(s% extra_jdot), "Size s% nz=", s% nz, &
                  !"Size mag_brk_jdot=", size(mag_brk_jdot)
            end if



!            s% x_ctrl(7) = v_esc
!            s% x_ctrl(8) = v_inf
!            s% x_ctrl(9) = eta_surf
!            s% x_ctrl(10) = j_dot
!            s% x_ctrl(11) = m_dot

         end if

         deallocate(mag_brk_jdot)

      end subroutine tfm_other_torque

      ! Collect information about the outermost convection zone
      subroutine get_convective_info(s, cz_info)
         type (star_info), pointer, intent(in) :: s
         integer :: i, k, nz, n_conv_bdy
         type (conv_zone_info), pointer, intent(out) :: cz_info

         nz = s% nz
         ! boundaries of regions with mixing_type = convective_mixing
         n_conv_bdy = s% num_conv_boundaries

         ! convection regions
         i = s% n_conv_regions


         ! Reset bottom and top zone, mass and radius values
         cz_info% top_zone = 0
         cz_info% bot_zone = 0

         cz_info% top_mass = 0.0
         cz_info% bot_mass = 0.0

         cz_info% top_radius = 0.0
         cz_info% bot_radius = 0.0
         
         ! mixing regions (from surface inward)
         ! check the outermost convection zone
         ! if dM_convenv/M < 1d-8, there's no conv env.
         if (s% n_conv_regions > 0) then
             if ((s% cz_top_mass(i) / s% mstar > 0.99d0) .and. &
             ((s% cz_top_mass(i) - s% cz_bot_mass(i)) / s% mstar > 1d-11)) then 
             
                 cz_info% bot_mass = s% cz_bot_mass(i)
                 cz_info% top_mass = s% cz_top_mass(i)                 

                 !get top radius information
                 !start from k=2 (second most outer zone) in order to access k-1
                 !iterate till the mass of zone k is smaller than the mass of cz limit
                 do k=2,nz
                     if (s% m(k) < cz_info% top_mass) then 
                         cz_info% top_radius = s% r(k-1)
                         cz_info% top_zone = k-1
                         exit
                     end if
                 end do

                 !get bottom radius information
                 do k=2,nz 
                     if (s% m(k) < cz_info% bot_mass) then 
                         cz_info% bot_radius = s% r(k-1)
                         cz_info% bot_zone = k-1
                         exit
                     end if
                 end do

                 
                 !if the star is fully convective, then the bottom boundary is the center
                 if ((cz_info% bot_zone == 0) .and. (cz_info% top_zone > 0)) then
                     cz_info% bot_zone = nz
                 end if

                 !calculate deltas
                 cz_info% d_mass = cz_info% top_mass - cz_info% bot_mass
                 cz_info% d_radius = cz_info% top_radius - cz_info% bot_radius
                 cz_info% d_zone = cz_info% bot_zone - cz_info% top_zone
             end if
         endif

         if (debug_get_cz_info) then
            write(*,*) "nz", nz, "num_cz=", s% n_conv_regions, &
               "bot_zone=", cz_info% bot_zone, "top_zone=", cz_info% top_zone, "d_zone=", cz_info% d_zone, &
               "bot_mass=", cz_info% bot_mass/msun, "top_mass=", cz_info% top_mass/msun, "d_mass=", cz_info% d_mass/msun, &
               "bot_radius=", cz_info% bot_radius/rsun, "top_radius=", cz_info% top_radius/rsun, "d_radius=", cz_info% d_radius/rsun
         end if

         
         !names(1) = 'conv_env_top_mass'
         !vals(1) = ocz_top_mass/msun
         !names(2) = 'conv_env_bot_mass'
         !vals(2) = ocz_bot_mass/msun
         !names(3) = 'conv_env_top_radius'
         !vals(3) = ocz_top_radius/rsun
         !names(4) = 'conv_env_bot_radius'
         !vals(4) = ocz_bot_radius/rsun
         !names(5) = 'conv_env_turnover_time_g'
         !vals(5) = ocz_turnover_time_g
         !names(6) = 'conv_env_turnover_time_l_b'
         !vals(6) = ocz_turnover_time_l_b
         !names(7) = 'conv_env_turnover_time_l_t'
         !vals(7) = ocz_turnover_time_l_t

      end subroutine get_convective_info

      subroutine distribute_j_dot(s, total_j_dot,cz_info, mb_jdot_list)
         type (star_info), pointer, intent(in) :: s
         real(dp), intent(in) :: total_j_dot
         type (conv_zone_info), pointer, intent(in) :: cz_info
         real(dp), dimension(:), pointer, intent(out) :: mb_jdot_list
         integer :: i
         real(dp) :: sum_jdot

         !By default, no lost of angular moment
         mb_jdot_list(:) = 0.0
         sum_jdot = 0.0
      
         !write(*,*) "bot_zone=", cz_info% bot_zone, "top_zone=", cz_info% top_zone
         !write(*,*) "sum_jdot=", sum_jdot
         do i = cz_info% top_zone, cz_info% bot_zone, 1
            !Here the jdot distribution strategy is defined
            !TODO externalize to a method, this will isolate the strategy implementation
            !Simple rule of three distribution based on the mass of the zone vs total cz mass
            mb_jdot_list(i) = (s% dm(i) * total_j_dot) / cz_info% d_mass
            sum_jdot = sum_jdot + mb_jdot_list(i)
            !write(*,*) "mb_jdot_list=", mb_jdot_list(i)
         end do
         !write(*,*) "sum_jdot=", sum_jdot
      end subroutine distribute_j_dot


      subroutine distribute_simple_j_dot(s, total_j_dot,cz_info, mb_jdot_list)
         type (star_info), pointer, intent(in) :: s
         real(dp), intent(in) :: total_j_dot
         type (conv_zone_info), pointer, intent(in) :: cz_info
         real(dp), dimension(:), pointer, intent(out) :: mb_jdot_list
         integer :: i

         !By default, no lost of angular moment
         mb_jdot_list(:) = 0.0
      
         !Just reduce a 1% for current angular moment
         do i = cz_info% top_zone, cz_info% bot_zone, 1
            mb_jdot_list(i) = -1.0 *s% omega(i) * s% i_rot(i) * 0.000000000000001
         end do

      end subroutine distribute_simple_j_dot

      
      integer function extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: core_ov_full_on, core_ov_full_off, frac, rot_full_off, rot_full_on, frac2, vct30, vct100
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_startup = 0
      if (.not. restart) then
         call alloc_extra_info(s)
      else                      ! it is a restart
         call unpack_extra_info(s)
      end if
      
!     set OPACITIES: Zbase for Type 2 Opacities automatically to the Z for the star
      s% Zbase = 1.0 - (s% job% initial_h1 + s% job% initial_h2 + &
      s% job% initial_he3 + s% job% initial_he4)
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(*,*) 'Zbase for Type 2 Opacities: ', s% Zbase
      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      
!     set ROTATION: extra param are set in inlist: star_job
      !rot_full_off = s% job% extras_rpar(1) !1.2
      !rot_full_on = s% job% extras_rpar(2) !1.8
      
      !if (s% job% extras_rpar(3) > 0.0) then
      !   if (s% star_mass < rot_full_off) then
      !      frac2 = 0
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !      write(*,*) 'no rotation'
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !   else if (s% star_mass >= rot_full_off .and. s% star_mass <= rot_full_on) then
      !      frac2 = (s% star_mass - rot_full_off) / &
      !      (rot_full_on - rot_full_off)
      !      frac2 = 0.5d0*(1 - cos(pi*frac2))
      !      s% job% set_near_zams_omega_div_omega_crit_steps = 10
      !      s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !      write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !   else
      !      frac2 = 1.0
      !      s% job% set_near_zams_omega_div_omega_crit_steps = 10
      !      s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2 !nominally 0.4
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !      write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
      !      write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !   end if
      !else
      !   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !   write(*,*) 'no rotation'
      !   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      !end if
      
      
      end function extras_startup
      
!     returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
      integer, intent(in) :: id, id_extra
      integer :: ierr, r, burn_category
      real(dp) :: envelope_mass_fraction, L_He, L_tot, orig_eta, target_eta, min_center_h1_for_diff, critmass, feh
      real(dp) :: category_factors(num_categories)
      real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
      real(dp), parameter :: new_varcontrol_target = 1d-3
      real(dp), parameter :: Zsol = 0.0142
      type (star_info), pointer :: s
      type (Net_General_Info), pointer :: g
      character (len=strlen) :: photoname
      
      ierr = 0	 
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going
      
      ierr = 0
      call get_net_ptr(s% net_handle, g, ierr)
      if (ierr /= 0) stop 'bad handle'	  
      
      
!     check DIFFUSION: to determine whether or not diffusion should happen
!     no diffusion for fully convective, post-MS, and mega-old models 
	  s% diffusion_dt_limit = 3.15d7
      if(abs(s% mass_conv_core - s% star_mass) < 1d-2) then ! => fully convective
         s% diffusion_dt_limit = huge_dt_limit
      end if
      if (s% star_age > 5d10) then !50 Gyr is really old
         s% diffusion_dt_limit = huge_dt_limit
      end if
      min_center_h1_for_diff = 1d-10
      if (s% center_h1 < min_center_h1_for_diff) then
         s% diffusion_dt_limit = huge_dt_limit
      end if
      
	  end function extras_check_model
      
      
      integer function how_many_extra_history_columns(id, id_extra)
      integer, intent(in) :: id, id_extra
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
      integer, intent(in) :: id, id_extra, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      end subroutine data_for_extra_history_columns
      
      
      integer function how_many_extra_profile_columns(id, id_extra)
      use star_def, only: star_info
      integer, intent(in) :: id, id_extra
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
      use star_def, only: star_info, maxlen_profile_column_name
      use const_def, only: dp
      integer, intent(in) :: id, id_extra, n, nz
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz,n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns
      
      
!     returns either keep_going or terminate.
!     note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
      integer, intent(in) :: id, id_extra
      integer :: ierr, k
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going
      call store_extra_info(s)
      
!     set BC: change to tables after running on simple photosphere
      if (s% model_number == 100) then
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'switching from simple photosphere to ', s% job% extras_cpar(1)
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         s% which_atm_option = s% job% extras_cpar(1)
      endif

!         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!         write(*,*) 'check extra_jdot '
!         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
!            do k=1, s% nz
!                  write(*,*) s% extra_jdot(k)
!            end do

      end function extras_finish_step
      
	  
	  subroutine Reimers_then_Blocker(id, Lsurf, Msurf, Rsurf, Tsurf, w, ierr)
      use star_def
      use chem_def, only: ih1, ihe4
      integer, intent(in) :: id
      real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf ! surface values (cgs)
!     NOTE: surface is outermost cell. not necessarily at photosphere.
!     NOTE: don't assume that vars are set at this point.
!     so if you want values other than those given as args,
!     you should use values from s% xh(:,:) and s% xa(:,:) only.
!     rather than things like s% Teff or s% lnT(:) which have not been set yet.
      real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
      integer, intent(out) :: ierr
      integer :: h1, he4
      real(dp) :: plain_reimers, reimers_w, blocker_w, center_h1, center_he4
	  type (star_info), pointer :: s
	  ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
	  
	  plain_reimers = 4d-13*(Lsurf*Rsurf/Msurf)/(Lsun*Rsun/Msun)
	  
	  reimers_w = plain_reimers * s% Reimers_scaling_factor
	  blocker_w = plain_reimers * s% Blocker_scaling_factor * &
               4.83d-9 * pow_cr(Msurf/Msun,-2.1d0) * pow_cr(Lsurf/Lsun,2.7d0)

          h1 = s% net_iso(ih1)
          he4 = s% net_iso(ihe4)
          center_h1 = s% xa(h1,s% nz)
          center_he4 = s% xa(he4,s% nz)

          !prevent the low mass RGBs from using Blocker
          if (center_h1 < 0.01d0 .and. center_he4 > 0.1d0) then
             w = reimers_w
          else 
             w = max(reimers_w, blocker_w)
          end if
	  
	  end subroutine Reimers_then_Blocker
      
      subroutine extras_after_evolve(id, id_extra, ierr)
      integer, intent(in) :: id, id_extra
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      subroutine alloc_extra_info(s)
      integer, parameter :: extra_info_alloc = 1
      type (star_info), pointer :: s
      call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info
      
      
      subroutine unpack_extra_info(s)
      integer, parameter :: extra_info_get = 2
      type (star_info), pointer :: s
      call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info
      
      
      subroutine store_extra_info(s)
      integer, parameter :: extra_info_put = 3
      type (star_info), pointer :: s
      call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info
      
      
      subroutine move_extra_info(s,op)
      integer, parameter :: extra_info_alloc = 1
      integer, parameter :: extra_info_get = 2
      integer, parameter :: extra_info_put = 3
      type (star_info), pointer :: s
      integer, intent(in) :: op
      
      integer :: i, j, num_ints, num_dbls, ierr
      
      i = 0
!     call move_int or move_flg    
      num_ints = i
      
      i = 0
!     call move_dbl       
      
      num_dbls = i
      
      if (op /= extra_info_alloc) return
      if (num_ints == 0 .and. num_dbls == 0) return
      
      ierr = 0
      call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
      if (ierr /= 0) then
         write(*,*) 'failed in star_alloc_extras'
         write(*,*) 'alloc_extras num_ints', num_ints
         write(*,*) 'alloc_extras num_dbls', num_dbls
         stop 1
      end if
      
      contains
      
      subroutine move_dbl(dbl)
      real(dp) :: dbl
      i = i+1
      select case (op)
      case (extra_info_get)
         dbl = s% extra_work(i)
      case (extra_info_put)
         s% extra_work(i) = dbl
      end select
      end subroutine move_dbl
      
      subroutine move_int(int)
      integer :: int
      i = i+1
      select case (op)
      case (extra_info_get)
         int = s% extra_iwork(i)
      case (extra_info_put)
         s% extra_iwork(i) = int
      end select
      end subroutine move_int
      
      subroutine move_flg(flg)
      logical :: flg
      i = i+1
      select case (op)
      case (extra_info_get)
         flg = (s% extra_iwork(i) /= 0)
      case (extra_info_put)
         if (flg) then
            s% extra_iwork(i) = 1
         else
            s% extra_iwork(i) = 0
         end if
      end select
      end subroutine move_flg
      
      end subroutine move_extra_info
      
      
      end module run_star_extras
