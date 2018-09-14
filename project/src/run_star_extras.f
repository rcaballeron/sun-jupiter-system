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
      rot_full_off = s% job% extras_rpar(1) !1.2
      rot_full_on = s% job% extras_rpar(2) !1.8
      
      if (s% job% extras_rpar(3) > 0.0) then
         if (s% star_mass < rot_full_off) then
            frac2 = 0
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'no rotation'
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         else if (s% star_mass >= rot_full_off .and. s% star_mass <= rot_full_on) then
            frac2 = (s% star_mass - rot_full_off) / &
            (rot_full_on - rot_full_off)
            frac2 = 0.5d0*(1 - cos(pi*frac2))
            s% job% set_near_zams_omega_div_omega_crit_steps = 10
            s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         else
            frac2 = 1.0
            s% job% set_near_zams_omega_div_omega_crit_steps = 10
            s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2 !nominally 0.4
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            write(*,*) 'new omega_div_omega_crit, fraction', s% job% new_omega_div_omega_crit, frac2
            write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         end if
      else
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) 'no rotation'
         write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      end if
      
      
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
      integer :: ierr
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
