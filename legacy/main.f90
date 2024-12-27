! Start of Main Program
! This is a program to solve a variable-density backward facing 
! step problem with a heated wall using the finite-volume method on
! a co-located mesh and the SIMPLE algorithm

    program main

    USE precisions
    USE geometry
    USE flow
    USE constants
    USE control_parameters

! Declaration of variables
    IMPLICIT NONE
    INTEGER(int_p) :: i,j,iter
    REAL(real_p) :: xx,yy
    
! Open output and residual files
    OPEN(unit=10, file="hw7.out", status = "unknown")
    OPEN(unit=20, file="hw7.rsl", status = "unknown")

! Mesh calculation
    dy = height_in/ny_in
    dx = length_in/nx_in
    ar = dx/dy
    ra = one/ar

! Initialize arrays
    uu=u_in;vv=zero;pres=zero;tem=300.0d0
    uface(:,:,:) = u_in; vface(:,:,:) = zero

! Properties setting
    mu = mu_const
    cp = cp_const
    IF(comprs)THEN
      rho = (p_ref+pres)/(R_gas*tem)
    ELSE
      rho = rho_const
      rho_l = rho_const
      rho_b = rho_const
      rho_t = rho_const
    ENDIF

! Boundary conditions
    uu_l(:) = u_in; uu_r = zero; uu_b(:) = zero; uu_t(:) = zero
    vv_l(:) = zero; vv_r = zero; vv_b(:) = zero; vv_t(:) = zero
    uu_l(1:ny_in) = zero; pres_r(:) = zero
    DO j = ny_in+1,ny
      uface(1,j,2) = uu_l(j)
    ENDDO

! Start of Outer Loop Iterations

    DO iter = 1,iter_global
      
      CALL links(ivar_u)

      CALL solve_x_mom

      CALL solve_y_mom

      CALL face_velocity   ! Face Velocity using PWIM

      CALL solve_pp

      CALL uv_correct

      CALL pres_correct
     
      WRITE(*,11) iter, res_u,res_v,res_p
      WRITE(20,11) iter, res_u,res_v,res_p

      IF(res_u < tol_outer .and. res_v < tol_outer  &
            .and. res_p < tol_outer)EXIT

    ENDDO

! Output for Tecplot

    WRITE(*,*) ""
    WRITE(10,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
    WRITE(10,*) 'ZONE I=', nx+2, ', J=', ny+2, ', DATAPACKING=POINT'
    DO i = 1,nx_in
      DO j = 1,ny_in
        uu(i,j) = zero
        vv(i,j) = zero
        pres(i,j) = zero
      ENDDO
    ENDDO
    WRITE(10,12) zero,zero,zero,zero,zero
    xx = - half*dx
    DO i = 1,nx
      xx = xx + dx
      WRITE(10,12) xx,zero,zero,zero,pres(i,1)
    ENDDO
    WRITE(10,12) xx+half*dx,zero,zero,zero,pres(nx,1)
    yy = - half*dy
    DO j = 1,ny
      yy = yy + dy
      WRITE(10,12) zero,yy,uu_l(j),vv_l(j),pres(1,j)
      xx = - half*dx
      DO i = 1,nx
        xx = xx + dx
        WRITE(10,12) xx,yy,uu(i,j),vv(i,j),pres(i,j)
      ENDDO
      WRITE(10,12) xx+half*dx,yy,uu(nx,j),vv(nx,j),pres(nx,j)
    ENDDO
    WRITE(10,12) zero,yy+half*dy,zero,zero,pres(1,ny)
    xx = - half*dx
    DO i = 1,nx
      xx = xx + dx
      WRITE(10,12) xx,yy+half*dy,zero,zero,pres(i,ny)
    ENDDO
    WRITE(10,12) xx+half*dx,yy+half*dy,zero,zero,pres(nx,ny)

    CLOSE(unit=10)
    CLOSE(unit=20)

 11 FORMAT(I6,1x,3(E14.6))
 12 FORMAT(5(1X,E14.6))

    END Program main

