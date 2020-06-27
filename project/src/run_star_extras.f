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
      use chem_def, only: ipp, icno
      
      implicit none
      
      real(dp) :: original_diffusion_dt_limit
      real(dp) :: burn_check = 0.0
      real(dp) :: eps_threshold = 0.1d-100
      real(dp) :: disk_lt
      real(dp) :: disk_omega
      logical :: debug_use_other_torque = .false.
      logical :: debug_reset_other_torque = .false.
      logical :: debug_get_cz_info = .false.
      logical :: debug_get_core_info = .false.
      logical :: keep_on_rad_core = .false.
      logical :: rad_core_developed = .false.
      integer :: idx_low_x_ctrl = 20
      integer :: idx_high_x_ctrl = 46

      type star_zone_info
       real(dp) :: &
         top_radius, &
         bot_radius, &
         d_radius, &
         top_mass, &
         bot_mass, &
         d_mass, &
         top_vrot, &
         bot_vrot, &
         d_vrot, &
         top_omega, &
         bot_omega, &
         d_omega
       integer :: &
         top_zone, &
         bot_zone, &
         d_zone
      end type star_zone_info

      type (star_zone_info), dimension(:), allocatable :: conv_zone_list

      
!     these routines are called by the standard run_star check_model
      contains
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         write(*,*) '*******  Version: 2.2.0'
         write(*,*) '*******  09.06.2020: Disk locking & Magnetic braking routines'
         write(*,*) '*******  Roque Caballero'
         
         original_diffusion_dt_limit = s% diffusion_dt_limit
         !s% other_wind => Reimers_then_VW
         s% other_wind => Reimers_then_Blocker
         !s% other_torque => other_torque_mag_brk
         s% other_torque => other_torque_hook

         !debug flags
         debug_use_other_torque = s% x_logical_ctrl(1)
         debug_reset_other_torque = s% x_logical_ctrl(2)
         debug_get_cz_info = s% x_logical_ctrl(3)
         debug_get_core_info = s% x_logical_ctrl(4)

         !If true, once the radiative core is developed, report always true
         !in is_radiative_core function
         keep_on_rad_core = s% x_logical_ctrl(6)

         !eps thershold
         eps_threshold = s% x_ctrl(2)

         !disk locking
         disk_lt = s% x_ctrl(3)
         disk_omega = s% x_ctrl(4)


      
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

      subroutine other_torque_hook(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         !real(dp) :: disk_lt

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !disk_lt = s% x_ctrl(3)

         if ((disk_lt > 0) .and. (s% star_age < disk_lt) .and. (s% omega(1) > disk_omega)) then
            call other_torque_disk_lock(id, ierr)
         else
            call other_torque_mag_brk(id, ierr)
         endif
         
      end subroutine other_torque_hook

      ! This routine implements a disk locking effect based on 
      ! Marc Pinsonneault Evolution of low mass stars MESA 2019 
      ! Summer school lecture.
      subroutine other_torque_disk_lock(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: j

         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         do j=1,s% nz
            s% omega(j) = disk_omega
         end do

         s% extra_omegadot(:) = 0
         s% extra_jdot(:) = 0

         if (debug_use_other_torque) then
               write(*,*) "other_torque_disk_lock invoked:"
               do j=1, s% nz
                     write(*,*) "omega(", j, ")=", s% omega(j), "i_rot(", j, ")=", s% i_rot(j)
               end do
         end if
      end subroutine other_torque_disk_lock

      ! This routine implements a magnetic braking effect based on
      ! Matteos MESA star summer school. Additionally, this implementation
      ! distribute among the zones which conforms the convective shell the
      ! loss of angular momentum
      subroutine other_torque_mag_brk(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: r_st, m_st, i_st, alfven_r
         real(dp) :: j_dot, omega_surf, m_dot, eta_surf, v_inf, v_esc, B, jdot_routine
         type (star_zone_info), target :: sz_info, core_info
         type (star_zone_info), pointer :: sz_info_ptr, core_info_ptr
         real(dp), dimension(:), pointer :: mag_brk_jdot
         logical :: only_cz !controls if jdot distribution must only affect the convective zone
         logical :: wait_rad_core !controls if jdot distribution must wait till a radiative core is develop
         integer :: activated !signals when the jdot routine is activated
         
         !Pointer to structure which conveys information about the convectice and core zones
         sz_info_ptr => sz_info
         core_info_ptr => core_info

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         allocate(mag_brk_jdot(s% nz))
         s% extra_jdot(:) = 0
         s% extra_omegadot(:) = 0
         activated = 0
         call reset_x_ctrl(s, idx_low_x_ctrl, idx_high_x_ctrl)
         call reset_core_info(core_info_ptr)
         call reset_convective_info(sz_info_ptr)

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
         !v_esc = photospheric escape velocity (cm/s)
         !v_esc = 618*((r_sol/r_st)*(m_st/m_sol))^(1/2)
         !loss of angular momentum
         !j_dot = 2/3*m_dot*omega_surf*r_st*r_st*eta_surf

         !Wait till radiative core is develop?
         wait_rad_core = s% x_logical_ctrl(5)

         !Get information about the convective zone and core

         !Loss of angular momentum distribution method
         jdot_routine = s% x_integer_ctrl(1)
         only_cz = .true.
         if (jdot_routine /= 0) then
               only_cz = .false.
         end if

         !Get information about the convective zone
         if (only_cz) then
            !Get information about the outter convective zone
            call get_convective_info(s, sz_info_ptr)
         else
            !Get information about the outter convective zone till star surface
            call get_convective_to_surf_info(s, sz_info_ptr)
         end if

         !Get information about the core
         call get_core_info(s, core_info_ptr)


         !Calculate amount of loss of angular moment

         !TODO move out the reading of inlist parameters
         !Magentic field intensity
         B = s% x_ctrl(1)
         !B_new = B * (Rsun/r_st)**3 !q=3 para un dipolo simple

         !Star data
         r_st = s% r(1)
         m_st = s% m(1)
         omega_surf = s% omega(1)
         !omega_surf = s% omega_avg_surf

         ! escape and infinite velocities
         v_esc = (618 * ((Rsun/r_st)*(m_st/Msun))**0.5) * 100000 !100000 transform from km/s to cm/s
         v_inf = 1.92 * v_esc


         !m_dot = s% star_mdot !This gives the mass loss rate in Mstar/year
         m_dot = s% mstar_dot !This in g/s

         !eta = ( s% photosphere_r * rsun * bfield )**2.0 / (abs( s% mstar_dot ) * vinf)
         !with vinf in cm/s (the wind terminal velocity) rsun in cm, bfield in Gauss
         !rest are MESA vars accessible through the star structure in run_star_extras.f
         eta_surf = ((r_st * B)**2)/(abs(m_dot) * v_inf)
         !eta_surf = ((s% photosphere_r * rsun * B)**2)/(abs(m_dot) * v_inf)

         j_dot = two_thirds * m_dot * omega_surf * (r_st**2) * eta_surf !Formula 2.3 Cantiello's MESA assigment
         !alfven_r = 50.0 * Rsun
         !j_dot = two_thirds * m_dot * omega_surf * (alfven_r**2) * eta_surf !Formula 2.1 Cantiello's MESA assigment

         ! The magnetic braking routine is activated under the following conditions:
         ! - use_other_torque flag is activated in inlist
         ! AND
         ! - the star is losing mass
         ! AND
         ! - the magnetic field intensitive is bigger than 0.0 (allow to execute the routine if use_other_torque=.true.)
         ! AND
         ! (
         !   - a radiative core developed isn't required (configured in inlist)
         !   OR
         !   - a radiative core is required AND this was developed
         !)
         if ((s% use_other_torque) .and. (s% mstar_dot < 0.0) .and. (B > 0.0) .and. &
            (.not. wait_rad_core .or. (wait_rad_core .and. is_core_rad(s)))) then
            activated = 1

            !Distribute the loss of angular momentum
            call distribute_j_dot(s, j_dot, sz_info_ptr, mag_brk_jdot)

            !It happens that s% extra_jdot is longer than s% nz but mag_brk_jdot is just defined
            !for s% nz elements
            s% extra_jdot(1:s% nz) = mag_brk_jdot

            if (debug_use_other_torque) then
                  !write(*,*) "photosphere_r", s% photosphere_r*rsun, "r(1)", s% r(1), "r(nz)", s% r(s% nz), &
                  !"j_dot", j_dot, "j_dot(nz)", mag_brk_jdot(s% nz), "j_dot(1)", mag_brk_jdot(1)

                  do k=1, size(mag_brk_jdot)
                        if (s% extra_jdot(k) < 0.0) then
                              write(*,*) "jdot(k)=", s% extra_jdot(k), "omega(k)*i_rot(k)=", s% omega(k) * s% i_rot(k), &
                        "omega(k)=", s% omega(k), "i_rot(k)=", s% i_rot(k)
                        end if
                  end do

                  if (debug_reset_other_torque) then
                        !Just only print out the values, don't pass them back to MESA
                        s% extra_jdot(:) = 0.0
                  end if
            end if
         end if

         !We abuse the x_ctrl array for storing temporary the values to be reported in history file
         s% x_ctrl(20) = v_esc
         s% x_ctrl(21) = v_inf
         s% x_ctrl(22) = eta_surf
         s% x_ctrl(23) = j_dot
         s% x_ctrl(24) = m_dot

         !Routine activated
         s% x_ctrl(25) = activated

         !Convective zone info
         s% x_ctrl(26) = sz_info_ptr% top_radius / r_st
         s% x_ctrl(27) = sz_info_ptr% bot_radius / r_st
         s% x_ctrl(28) = sz_info_ptr% top_mass / m_st
         s% x_ctrl(29) = sz_info_ptr% bot_mass/ m_st
         s% x_ctrl(30) = sz_info_ptr% top_zone
         s% x_ctrl(31) = sz_info_ptr% bot_zone
         s% x_ctrl(32) = sz_info_ptr% top_vrot
         s% x_ctrl(33) = sz_info_ptr% bot_vrot
         s% x_ctrl(34) = sz_info_ptr% top_omega
         s% x_ctrl(35) = sz_info_ptr% bot_omega

         !Core info
         s% x_ctrl(36) = eps_threshold
         s% x_ctrl(37) = core_info_ptr% top_radius / r_st
         s% x_ctrl(38) = core_info_ptr% bot_radius / r_st
         s% x_ctrl(39) = core_info_ptr% top_mass / m_st
         s% x_ctrl(40) = core_info_ptr% bot_mass/ m_st
         s% x_ctrl(41) = core_info_ptr% top_zone
         s% x_ctrl(42) = core_info_ptr% bot_zone
         s% x_ctrl(43) = core_info_ptr% top_vrot
         s% x_ctrl(44) = core_info_ptr% bot_vrot
         s% x_ctrl(45) = core_info_ptr% top_omega
         s% x_ctrl(46) = core_info_ptr% bot_omega


         deallocate(mag_brk_jdot)
      end subroutine other_torque_mag_brk

      subroutine reset_x_ctrl(s, k_ini, k_end)
         type (star_info), pointer, intent(in) :: s
         integer, intent(in) :: k_ini, k_end
         integer :: k

         do k=k_ini, k_end
            s% x_ctrl(k) = 0
         end do
      end subroutine reset_x_ctrl

      subroutine reset_core_info(core_info)
         type (star_zone_info), pointer, intent(out) :: core_info

         call reset_star_zone_info(core_info)
      end subroutine reset_core_info

      subroutine reset_convective_info(conv_info)
         type (star_zone_info), pointer, intent(out) :: conv_info

         call reset_star_zone_info(conv_info)
      end subroutine reset_convective_info

      subroutine reset_star_zone_info(sz_info)
         type (star_zone_info), pointer, intent(out) :: sz_info

         ! Reset bottom and top zone, mass and radius values
         sz_info% top_zone = 0
         sz_info% bot_zone = 0
         sz_info% d_zone = 0

         sz_info% top_mass = 0.0
         sz_info% bot_mass = 0.0
         sz_info% d_mass = 0.0

         sz_info% top_radius = 0.0
         sz_info% bot_radius = 0.0
         sz_info% d_radius = 0.0

         sz_info% top_vrot = 0.0
         sz_info% bot_vrot = 0.0

         sz_info% top_omega = 0.0
         sz_info% bot_omega = 0.0
         sz_info% d_omega = 0.0

      end subroutine reset_star_zone_info



      ! We define that a radiative core is developed if a convective
      ! one isn't present.
      function is_core_rad(s) result(flag)
         type (star_info), pointer, intent(in) :: s
         logical :: flag

         ! Once developed, return always true
         if (rad_core_developed) then
            flag = .true.
         ! Don't check if the star is too young
         else
            if (s% star_age > 1.0e5) then
                if (s% mass_conv_core > 0.0) then
                    flag = .false.
                else
                    if (keep_on_rad_core) then
                        rad_core_developed = .true.
                    end if
                    flag = .true.
                end if
            else
                flag = .false.
            end if
         end if
      end function

      ! Retrive information about the outermost convective region of the star
      subroutine get_outermost_conv_info(s, sz_info, only_cz)
         type (star_info), pointer, intent(in) :: s
         logical, intent(in) :: only_cz
         integer :: i, j, k, nz, n_conv_bdy
         type (star_zone_info), pointer, intent(out) :: sz_info
         real(dp) zone_top_mass

         nz = s% nz
         ! boundaries of regions with mixing_type = convective_mixing
         n_conv_bdy = s% num_conv_boundaries

         ! convection regions
         i = s% n_conv_regions


         ! mixing regions (from surface inward)
         ! check the outermost convection zone
         ! if dM_convenv/M < 1d-8, there's no conv env.
         if (s% n_conv_regions > 0) then

            ! Check if only calculate for the convective zone or get till the surface
            if (only_cz) then
                  zone_top_mass = s% cz_top_mass(i)
            else 
                  zone_top_mass = s% m(1)
            end if
         

             ! with 'i', we peek the outermost convection zone and check if most of the
             ! star's mass is below it. Basically it checks that the convection zone
             ! is close to the surface
             if ((zone_top_mass / s% mstar > 0.99d0) .and. &
             ! additionally, the convection zone must have a minimum amount of mass
             ((zone_top_mass - s% cz_bot_mass(i)) / s% mstar > 1d-11)) then 
             
                 sz_info% bot_mass = s% cz_bot_mass(i)
                 sz_info% top_mass = zone_top_mass

                 !get top radius information
                 !start from k=2 (second most outer zone) in order to access k-1
                 !iterate till the mass of zone k is smaller than the mass of cz limit
                 if (only_cz) then
                        do k=2,nz
                              if (s% m(k) < sz_info% top_mass) then 
                                    sz_info% top_radius = s% r(k-1)
                                    sz_info% top_zone = k-1
                                    exit
                              end if
                        end do
                  else 
                        sz_info% top_radius = s% r(1)
                        sz_info% top_zone = 1
                  end if


                 !get bottom radius information
                 do k=2,nz 
                     if (s% m(k) < sz_info% bot_mass) then 
                         sz_info% bot_radius = s% r(k-1)
                         sz_info% bot_zone = k-1
                         exit
                     end if
                 end do

                 
                 !if the star is fully convective, then the bottom boundary is the center
                 !if sz_info% bot_zone remains equal to 0 means that the bottom limit were not yet found
                 if ((sz_info% bot_zone == 0) .and. (sz_info% top_zone > 0)) then
                     sz_info% bot_zone = nz
                 end if

                 !retrieve angular velocity: omega
                 sz_info% top_omega = s% omega(sz_info% top_zone) ! rad/sec
                 sz_info% bot_omega = s% omega(sz_info% bot_zone) ! rad/sec

                 !retrieve rotational velocities: omega * radius
                 sz_info% top_vrot = s% omega(sz_info% top_zone)* s% r(sz_info% top_zone)*1d-5 ! km/sec
                 sz_info% bot_vrot = s% omega(sz_info% bot_zone)* s% r(sz_info% bot_zone)*1d-5 ! km/sec

                 !calculate deltas
                 sz_info% d_mass = sz_info% top_mass - sz_info% bot_mass
                 sz_info% d_radius = sz_info% top_radius - sz_info% bot_radius
                 sz_info% d_zone = sz_info% bot_zone - sz_info% top_zone
                 sz_info% d_omega = sz_info% top_omega - sz_info% bot_omega
                 sz_info% d_vrot = sz_info% top_vrot - sz_info% bot_vrot
             end if
         endif

         if (debug_get_cz_info) then
            write(*,*) "nz", nz, "num_cz=", s% n_conv_regions, &
               "bot_zone=", sz_info% bot_zone, "top_zone=", sz_info% top_zone, "d_zone=", sz_info% d_zone, &
               "bot_mass=", sz_info% bot_mass/msun, "top_mass=", sz_info% top_mass/msun, "d_mass=", sz_info% d_mass/msun, &
               "bot_radius=", sz_info% bot_radius/rsun, "top_radius=", sz_info% top_radius/rsun, &
               "d_radius=", sz_info% d_radius/rsun, &
               "bot_omega=", sz_info% bot_omega, "top_omega=", sz_info% top_omega, "d_omega=", sz_info% d_omega, &
               "bot_vrot=", sz_info% bot_vrot, "top_vrot=", sz_info% top_vrot, "d_vrot=", sz_info% d_vrot

         end if
      end subroutine get_outermost_conv_info
       

      ! Collect information about the outermost convection zone
      subroutine get_convective_info(s, sz_info)
         type (star_info), pointer, intent(in) :: s
         type (star_zone_info), pointer, intent(out) :: sz_info

            call get_outermost_conv_info(s, sz_info, .true.)
      end subroutine get_convective_info

      ! Collect information from the botton of the outermost convection zone till surface
      subroutine get_convective_to_surf_info(s, sz_info)
         type (star_info), pointer, intent(in) :: s
         type (star_zone_info), pointer, intent(out) :: sz_info

            call get_outermost_conv_info(s, sz_info, .false.)
      end subroutine get_convective_to_surf_info


      subroutine distribute_j_dot(s, total_j_dot, sz_info, mb_jdot_list)
         type (star_info), pointer, intent(in) :: s
         real(dp), intent(in) :: total_j_dot
         type (star_zone_info), pointer, intent(in) :: sz_info
         real(dp), dimension(:), pointer, intent(out) :: mb_jdot_list
         integer :: k
         real(dp) :: sum_jdot, dm_jdot, dm_bar_jdot, factor

         !By default, no lost of angular moment
         mb_jdot_list(:) = 0.0
      
         do k = sz_info% top_zone, sz_info% bot_zone, 1
            !Here the jdot distribution strategy is defined
            !TODO externalize to a method, this will isolate the strategy implementation
            !Simple rule of three distribution loss of angular momentum based on the 
            !angular momentum of the zone vs total angular momentum of the convective zone

            !IMPORTANT: Don't forget to divide by dm(k) in order to get an "specific" jdot
            mb_jdot_list(k) = ((s% dm(k) * s% r(k)**2 * total_j_dot) / &
                  (sz_info% d_mass * sz_info% d_radius**2)) / s% dm(k)

         end do
      end subroutine distribute_j_dot

      subroutine distribute_all_zones_j_dot(s, total_j_dot, mb_jdot_list)
         type (star_info), pointer, intent(in) :: s
         real(dp), intent(in) :: total_j_dot
         real(dp), dimension(:), pointer, intent(out) :: mb_jdot_list
         integer :: k

         !By default, no lost of angular moment
         mb_jdot_list(:) = 0.0
      
         do k = 1, s% nz, 1
            mb_jdot_list(k) = ((s% dm(k) * s% r(k)**2 * total_j_dot) / (s% m(1) * s% r(1)**2)) / s% dm(k)
            !Using the photsosphere radius
            !mb_jdot_list(k) = (s% dm(k) * s% r(k)**2 * total_j_dot) / (s% m(1) * (s% photosphere_r * rsun)**2)
         end do
      end subroutine distribute_all_zones_j_dot


      ! all the zones in my star in which H fusion occurs.
      subroutine get_core_info(s, core_info)
         type (star_info), pointer, intent(in) :: s
               integer :: i, j, k, nz
         type (star_zone_info), pointer, intent(out) :: core_info
         real(dp) zone_top_mass
         logical fusion_h

         nz = s% nz

         fusion_h = .false.
         do k=nz, 1, -1
            ! When the relese nuclear energy is below than the value in the if comparison, it means that there is not H fusion
            ! The top zone of the core was found.
            ! NB: nz is the inner most zone and 1 the outter most one.
            if (s% eps_nuc_categories(ipp, k) + s% eps_nuc_categories(icno, k) < eps_threshold) then
               if (k /= nz) then
                  core_info% top_radius = s% r(k+1)
                  core_info% bot_radius = s% r(nz)
                  core_info% top_zone = k+1
                  core_info% bot_zone = nz
                  core_info% top_mass = s% m(k+1)
                  core_info% bot_mass = s% m(nz)

                  !retrieve angular velocity: omega
                  core_info% top_omega = s% omega(core_info% top_zone) ! rad/sec
                  core_info% bot_omega = s% omega(core_info% bot_zone) ! rad/sec

                  !retrieve rotational velocities
                  core_info% top_vrot = s% omega(core_info% top_zone)*s% r(core_info% top_zone)*1d-5 ! km/sec
                  core_info% bot_vrot = s% omega(core_info% bot_zone)*s% r(core_info% bot_zone)*1d-5 ! km/sec

                  !calculate deltas
                  core_info% d_mass = core_info% top_mass - core_info% bot_mass
                  core_info% d_radius = core_info% top_radius - core_info% bot_radius
                  core_info% d_zone = core_info% bot_zone - core_info% top_zone
                  core_info% d_omega = core_info% top_omega - core_info% bot_omega
                  core_info% d_vrot = core_info% top_vrot - core_info% bot_vrot

                  fusion_h = .true.
               end if

               exit
            end if
         end do

         if (debug_get_core_info) then
            write(*,*) "nz", nz, "fusion_h=", fusion_h, &
               "bot_zone=", core_info% bot_zone, "top_zone=", core_info% top_zone, "d_zone=", core_info% d_zone, &
               "bot_mass=", core_info% bot_mass/Msun, "top_mass=", core_info% top_mass/Msun, "d_mass=", core_info% d_mass/Msun, &
               "bot_radius=", core_info% bot_radius/Rsun, "top_radius=", core_info% top_radius/Rsun, & 
               "d_radius=", core_info% d_radius/Rsun, &
               "bot_omega=", core_info% bot_omega, "top_omega=", core_info% top_omega, "d_omega=", core_info% d_omega, &
               "bot_vrot=", core_info% bot_vrot, "top_vrot=", core_info% top_vrot, "d_vrot=", core_info% d_vrot
               
         end if
      end subroutine get_core_info
      
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

            !if (s% use_other_torque .eqv. .true.) then
                  how_many_extra_history_columns = 28
            !end if
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
      integer, intent(in) :: id, id_extra, n
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      integer :: i
      type (star_info), pointer :: s

      ierr = 0
      i = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      !if (s% use_other_torque .eqv. .true.) then
        !routine activated
        names(1) = 'activated'
        vals(1) = s% x_ctrl(25) 

        !magnetic field torque
        names(2) = 'B'
        vals(2) = s% x_ctrl(1) 

        !photospheric escape velocity
        names(3) = 'v_esc'
        vals(3) = s% x_ctrl(20)
                  
        !terminal velocity
        names(4) = 'v_inf'
        vals(4) = s% x_ctrl(21)
                  
        !Wind-confinmenet parameter
        names(5) = 'eta_surf'
        vals(5) = s% x_ctrl(22)
                  
        !lost fo angular momentum
        names(6) = 'j_dot'
        vals(6) = s% x_ctrl(23)

        !mass lost
        names(7) = 'm_dot'                  
        vals(7) = s% x_ctrl(24)

        !sz_top_radius
        names(8) = 'sz_top_radius'
        vals(8) = s% x_ctrl(26)

        !sz_bottom_radius
        names(9) = 'sz_bot_radius'
        vals(9) = s% x_ctrl(27)

        !sz_top_mass
        names(10) = 'sz_top_mass'
        vals(10) = s% x_ctrl(28)

        !sz_bottom_mass
        names(11) = 'sz_bot_mass'                  
        vals(11) = s% x_ctrl(29)

        !sz_top_zone
        names(12) = 'sz_top_zone'                  
        vals(12) = s% x_ctrl(30)

        !sz_bottom_zone
        names(13) = 'sz_bot_zone'                  
        vals(13) = s% x_ctrl(31)

        !sz_top_vrot
        names(14) = 'sz_top_vrot'                  
        vals(14) = s% x_ctrl(32)

        !sz_bot_vrot
        names(15) = 'sz_bot_vrot'                  
        vals(15) = s% x_ctrl(33)

        !sz_top_omega
        names(16) = 'sz_top_omega'
        vals(16) = s% x_ctrl(34)

        !sz_bot_omega
        names(17) = 'sz_bot_omega'
        vals(17) = s% x_ctrl(35)

        !core_top_radius
        names(18) = 'eps_threshold'
        vals(18) = s% x_ctrl(36)

        !core_top_radius
        names(19) = 'core_top_radius'
        vals(19) = s% x_ctrl(37)

        !core_bottom_radius
        names(20) = 'core_bot_radius'
        vals(20) = s% x_ctrl(38)

        !core_top_mass
        names(21) = 'core_top_mass'
        vals(21) = s% x_ctrl(39)

        !core_bottom_mass
        names(22) = 'core_bot_mass'
        vals(22) = s% x_ctrl(40)

        !core_top_zone
        names(23) = 'core_top_zone'
        vals(23) = s% x_ctrl(41)

        !core_bottom_zone
        names(24) = 'core_bot_zone'
        vals(24) = s% x_ctrl(42)

        !core_top_vrot
        names(25) = 'core_top_vrot'
        vals(25) = s% x_ctrl(43)

        !core_bot_vrot
        names(26) = 'core_bot_vrot'
        vals(26) = s% x_ctrl(44)

        !core_top_omega
        names(27) = 'core_top_omega'
        vals(27) = s% x_ctrl(45)

        !core_bot_omega
        names(28) = 'core_bot_omega'
        vals(28) = s% x_ctrl(46)

      !end if

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

      if (s% use_other_torque .eqv. .true.) then 
            how_many_extra_profile_columns = 1
      end if

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

         names(1) = 'extra_jdot'
         do k = 1, nz
               vals(k,1) = s% extra_jdot(k)
         end do
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
