!***********************************************************************
! Declaration of Precision of Variables

    MODULE precisions

      IMPLICIT NONE
      SAVE
   
      INTEGER, PARAMETER :: int_p  = SELECTED_INT_KIND(8)
      INTEGER, PARAMETER :: real_p = SELECTED_REAL_KIND(8)
    
    END MODULE precisions
!***********************************************************************
! constants

    MODULE constants

      USE precisions
      IMPLICIT NONE
      SAVE

      REAL(real_p),PARAMETER :: one=1.0d0,zero=0.0d0,half=0.5d0, &
                                two=2.0d0,four=4.0d0,tiny=1.0d-40
    
    END MODULE constants

!***********************************************************************
! control parameters

    MODULE control_parameters

      USE precisions
      IMPLICIT NONE
      SAVE

      INTEGER(int_p),PARAMETER :: iter_mom = 2, iter_pp = 30, &
                                  iter_global = 2000
      REAL(real_p), PARAMETER :: relax_uv = 0.8d0, relax_p = 0.2d0, &
                                 tol_inner = 1.0d-3, tol_outer = 1.0d-10, &
                                 rin_uv = 0.1d0
    
    END MODULE control_parameters

!***********************************************************************
! Storage associated with geometry

    MODULE geometry

      USE precisions
      IMPLICIT NONE
      SAVE

      INTEGER(int_p), PARAMETER :: nx_in=40,ny_in=20,nx=200,ny=40
      REAL(real_p),PARAMETER :: length_in = 0.02d0, height_in = 0.01d0
      REAL(real_p) :: length_out,height_out

      REAL(real_p) :: ar,ra,dx,dy
      REAL(real_p), DIMENSION(nx) :: xc
      REAL(real_p), DIMENSION(ny) :: yc
    
    END MODULE geometry
!***********************************************************************
! Flow related variables

    MODULE flow

      USE precisions
      USE geometry, ONLY : nx,ny
      IMPLICIT NONE
      SAVE

      LOGICAL :: comprs = .false.
      INTEGER(int_p), PARAMETER :: ivar_u = 1,ivar_v = 2, &
                                   ivar_p = 3,ivar_t = 4
      REAL(real_p), PARAMETER :: u_in = 0.2d0, &   ! u_in = 0.1 => Re = 100
                                 p_ref = 1.0d5, &
                                 cp_const = 1000.0d0, R_gas = 8314.0d0, &
                                 rho_const = 1.0d0, grav = 0.0d0, &
                                 mu_const = 2.0d-5
      REAL(real_p), DIMENSION(nx,ny) :: mu,cond,cp,rho ! properties
      REAL(real_p), DIMENSION(nx,ny) :: uu,vv,pres,pp,tem
      REAL(real_p), DIMENSION(nx,ny,2) :: uface,vface
      REAL(real_p), DIMENSION(nx,ny,4) :: rhof
      REAL(real_p), DIMENSION(nx,ny) :: ao,aw,ae,an,as
      REAL(real_p), DIMENSION(ny) :: uu_l,uu_r,vv_l,vv_r
      REAL(real_p), DIMENSION(nx) :: uu_b,uu_t,vv_b,vv_t
      REAL(real_p), DIMENSION(nx) :: rho_b,rho_t
      REAL(real_p), DIMENSION(ny) :: rho_l
      REAL(real_p), DIMENSION(ny) :: pres_l, pres_r
      REAL(real_p) :: res_u,res_v,res_p
    
    END MODULE flow
