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
      
      implicit none


      
      ! these routines are called by the standard run_star check_model
      contains

      subroutine tfm_other_torque(id, ierr)
         use const_def
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: r_sol, m_sol, r_st, m_st
         real(dp) :: j_dot, omega_surf, m_dot, eta_surf, v_inf, v_esc, B
         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% extra_jdot(:) = 0
         s% extra_omegadot(:) = 0

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

         if (s% use_other_torque) then
            !Constants (cgs)
            r_sol = Rsun
            m_sol = Msun


            !Star data
            r_st = s% r(1)
            m_st = s% m(1)
            omega_surf = s% omega_avg_surf

            v_esc = 618 * (((r_sol/r_st)*(m_st/m_sol)))**0.5
            !v_esc = vesc / 100000 ! convert from cm/s to km/s.

            v_inf = 1.92 * v_esc

            B = s% x_ctrl(6)
            m_dot = s% star_mdot
            eta_surf = (r_st**2/B**2)/(m_dot * v_inf)

            j_dot = (2/3) * m_dot * omega_surf * r_st**2 * eta_surf

            write(*,*) "r_sol=", r_sol, "m_sol=", m_sol, "r_st=", r_st, "m_st=", m_st, &
               "v_esc=", v_esc, "v_inf", v_inf, "B", B, "m_dot", m_dot, "eta_surf", eta_surf, &
               "omega_surf", omega_surf, "j_dot", j_dot

            !s% x_ctrl(6) = B
            s% x_ctrl(7) = v_esc
            s% x_ctrl(8) = v_inf
            s% x_ctrl(9) = eta_surf
            s% x_ctrl(10) = j_dot
               


         end if

      end subroutine tfm_other_torque

      subroutine tfm_other_cgrav(id, ierr)
         use const_def, only: standard_cgrav
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, size_cgrav, size_nz
         real(dp) :: tfm_mp, tfm_dp, tfm_s_age, tfm_perc_dp, tfm_min_dp
         real(dp) :: tfm_lz_d, tfm_new_cgrav

         ierr = 0

         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ! star age, time in years
         tfm_s_age = s% x_ctrl(1)

         ! min distance, relative to tfm_dp, at which the planet must
         ! situated before activating its influence
         tfm_perc_dp = s% x_ctrl(2)  
         tfm_dp = s% x_ctrl(3)*100000 ! convert from km to cm
         tfm_mp = m_jupiter * s% x_ctrl(4) ! convert from Jupiter mass to grams
         tfm_min_dp = tfm_perc_dp * tfm_dp



         ! Get the size of the outtest star zone.
         tfm_lz_d = s% r(1)
         !write(*,*) "dist min=", tfm_min_dp, "dist lz", tfm_lz_d, "current dist=", tfm_dp - tfm_lz_d
         ! Activate planet influence only when star age reaches tfm_s_age
         ! and the planet is at least at tfm_min_dp distance
         if (s% use_other_cgrav .and. &
            (s% star_age  > tfm_s_age) .and. &
            ((tfm_dp - tfm_lz_d) >= tfm_min_dp)) then
            ! geff = g - G * Mp / (D-r)²
            ! geff = G' * m / r²
            ! G' * m / r² = g - G * Mp / (D-r)²
            ! g = G * m / r²
            !
            ! G' =  G (1 - ((Mp/m)*r²/(D-r)²))
            !PDTE de comprobar si utilzamos m(k) o dm(k)
            !masa contenida desde el centro de la estrella hasta el borde
            !de celda k o sólo la masa contenida en celda k.
            !Utilizamos la primera opción de las comentadas
            !size_cgrav = size(s% cgrav)
            
            !write(*,*) "num cells=", s% nz, "cgrav size", size_cgrav
            !write(*,*) "modified gravity constant"
            ! Get the number of star zones in the current model
            size_nz = s% nz
            do k = 1, size_nz 
               tfm_new_cgrav = standard_cgrav * (1 - ( (tfm_mp / s% m(k)) * s% r(k)**2 / (tfm_dp - s% r(k))**2 ) )
               s% cgrav(k) = tfm_new_cgrav
               !write(*,*) k,"m=", s% m(k), "r=",s% r(k),"G=", s% cgrav(k),  &
               !   "Mp/m=",(tfm_mp / s% m(k)),"r²/(D-r)²=",s% r(k)**2 / (tfm_dp - s% r(k))**2 
            end do
         else
            s% cgrav(:) = standard_cgrav
            !write(*,*) "standard gravity constant"
         end if

      end subroutine tfm_other_cgrav

      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0

          ! before we can use controls associated with the star we need to get access
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         s% other_cgrav => tfm_other_cgrav
         s% other_torque => tfm_other_torque
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         
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

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.       
            
      end subroutine extras_controls
      
      ! None of the following functions are called unless you set their
      ! function point in extras_control.
      
      
      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup
      

      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% use_other_cgrav .eqv. .true.) then 
            how_many_extra_history_columns = 5
         else
            how_many_extra_history_columns = 0
         end if

         if (s% use_other_torque .eqv. .true.) then 
            how_many_extra_history_columns = 11;
         end if

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
         
         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         if (s% use_other_cgrav .eqv. .true.) then

            !column 1, age threshold 
            names(1) = 'tfm_s_age'
            vals(1) = s% x_ctrl(1)

            !column 2, distance percentage to planet threshold
            names(2) = 'tfm_perc_dp'
            vals(2) = s% x_ctrl(2)

            !column 3, distance to planet
            names(3) = 'tfm_dp'
            vals(3) = s% x_ctrl(3)

            !column 4, mass of planet
            names(4) = 'tfm_mp'
            vals(4) = s% x_ctrl(4)

            !column 5, distance from the outtest zone of the star to the planet 
            names(5) = 'tfm_s_p_d'
            vals(5) = s% x_ctrl(3)*100000 - s% r(1)

         end if

         if (s% use_other_torque .eqv. .true.) then
            !column 1, age threshold 
            names(1) = 'tfm_s_age'
            vals(1) = 0.0

            !column 2, distance percentage to planet threshold
            names(2) = 'tfm_perc_dp'
            vals(2) = 0.0

            !column 3, distance to planet
            names(3) = 'tfm_dp'
            vals(3) = 0.0

            !column 4, mass of planet
            names(4) = 'tfm_mp'
            vals(4) = 0.0

            !column 5, distance from the outtest zone of the star to the planet 
            names(5) = 'tfm_s_p_d'
            vals(5) = 0.0


            !magnetic field torque
            names(6) = 'B'
            vals(6) = s% x_ctrl(6) 

            !photospheric escape velocity
            names(7) = 'v_esc'
            vals(7) = s% x_ctrl(7)
            
            !terminal velocity
            names(8) = 'v_inf'
            vals(8) = s% x_ctrl(8)
            
            !Wind-confinmenet parameter
            names(9) = 'eta_surf'
            vals(9) = s% x_ctrl(9)
            
            !angular momentum lost
            names(10) = 'j_dot'
            vals(10) = s% x_ctrl(10)

            names(11) = 'm_dot'
            vals(11) = s% star_mdot
         end if
         ierr = 0         
  
      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (s% use_other_cgrav .eqv. .true.) then
            how_many_extra_profile_columns = 1
         else
            how_many_extra_profile_columns = 0
         endif
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
         
         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

         if (s% use_other_cgrav .eqv. .true.) then
            !column 1
            names(1) = 'tfm_other_cgrav'
            do k = 1, nz
               vals(k,1) = s% cgrav(k)
            end do
            
            ierr = 0
         end if

         
      end subroutine data_for_extra_profile_columns
      

      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         logical :: tfm_start_planet_influence = .false.
         real :: tfm_Lnuc_div_L
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         !TFM
         !Put here the activation condition, for example, when
         !the star reaches a particular parameter
         ! activate other_energy after central hydrogen depletion
         !if (s% center_h1 < 0.01 ) s% use_other_mlt = .true.
         !if (s% center_h1 < 0.01 ) s% use_other_difussion = .true.

         !tfm_start_planet_influence = s% x_logical_ctrl(1)

         !if (tfm_start_planet_influence .and. (s% L_phot > 0d0)) then
         !   tfm_Lnuc_div_L = s% L_nuc_burn_total / s% L_phot
         !   if (tfm_Lnuc_div_L >= s% Lnuc_div_L_zams_limit) then
         !      s% use_other_mlt = .true.
         !   end if
         !end if

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
      
      ! routines for saving and restoring extra data so can do restarts
         
         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3
      
      
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
         ! call move_int or move_flg    
         num_ints = i
         
         i = 0
         ! call move_dbl       
         
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
      
