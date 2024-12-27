!---------------------------------------------------------------------
    SUBROUTINE solve_x_mom

! Solves X-Momentum equation using the ADI method

    USE precisions
    USE flow
    USE constants
    USE geometry
    USE control_parameters

    IMPLICIT NONE
    INTEGER(int_p) :: i,j,iter,il
    REAL(real_p) :: ce,cw,cn,cs,de,dw,dn,ds,res,pw,pe
    REAL(real_p), DIMENSION(nx) :: ax,bx,cx,ddx,solx
    REAL(real_p), DIMENSION(nx,ny) :: gam
!---------------------------------------------------------------------

    gam(:,:) = mu(:,:)

    DO iter = 1,iter_mom

! Calculate Residual of X-mom equation

      res = zero
      DO i=2,nx-1
        DO j=ny_in+2,ny-1
          res = res + (ao(i,j)*uu(i,j)+aw(i,j)*uu(i-1,j)+ &
                ae(i,j)*uu(i+1,j)+an(i,j)*uu(i,j+1)+ &
                as(i,j)*uu(i,j-1)-half*(pres(i-1,j)-pres(i+1,j))*dy)**2
        ENDDO
      ENDDO
      DO i=nx_in+2,nx-1
        DO j=2,ny_in+1
          res = res + (ao(i,j)*uu(i,j)+aw(i,j)*uu(i-1,j)+ &
                ae(i,j)*uu(i+1,j)+an(i,j)*uu(i,j+1)+ &
                as(i,j)*uu(i,j-1)-half*(pres(i-1,j)-pres(i+1,j))*dy)**2
        ENDDO
      ENDDO
      res = SQRT(MAX(res,tiny))
      IF(iter == 1) res_u = res

! Row-wise sweep

! Bottom row: Bottom Section
      j = 1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = nx_in+1
      il = i-nx_in
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = ae(i,j)
      ds = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
             + 8.0d0*ds*ar*uu_b(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j)
      DO i = nx_in+2,nx-1
        il = i-nx_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        ax(il-1) = aw(i,j)
        ds = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
               + 8.0d0*ds*ar*uu_b(i)/3.0d0
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                 aw(i,j)*uu(i-1,j)
      ENDDO
      i = nx
      il = i-nx_in
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = aw(i,j)
      ds = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = pres_r(j)
      bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
               + 8.0d0*ds*ar*uu_b(i)/3.0d0 
      bx(il) = bx(il) - aw(i,j)*uu(i-1,j) - ao(i,j)*uu(i,j)
  
      CALL TRI(nx-nx_in,ax,ddx,cx,bx,solx)
  
      DO i = nx_in+1,nx
        uu(i,j) = uu(i,j) + solx(i-nx_in)
      ENDDO
  
!   Interior rows: Bottom Section
      DO j = 2,ny_in
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        i = nx_in+1
        il = i-nx_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        dw = gam(i,j)
        pw = pres(i,j)
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                 -as(i,j)*uu(i,j-1) + 8.0d0*dw*ra*uu_l(j)/3.0d0
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j)
        DO i = nx_in+2,nx-1
          il = i-nx_in
          ddx(il) = ao(i,j)*(one+rin_uv)
          cx(il) = ae(i,j)
          ax(il-1) = aw(i,j)
          pw = half*(pres(i,j)+pres(i-1,j))
          pe = half*(pres(i,j)+pres(i+1,j))
          bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                   -as(i,j)*uu(i,j-1)
          bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                   aw(i,j)*uu(i-1,j)
        ENDDO
        i = nx
        il = i-nx_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        ax(il-1) = aw(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = pres_r(j)
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                 -as(i,j)*uu(i,j-1)
        bx(il) = bx(il) - aw(i,j)*uu(i-1,j) - ao(i,j)*uu(i,j)
    
        CALL TRI(nx-nx_in,ax,ddx,cx,bx,solx)
  
        DO i = nx_in+1,nx
          uu(i,j) = uu(i,j) + solx(i-nx_in)
        ENDDO
  
      ENDDO

! Bottom row: Top (Inlet) Section
      j = ny_in+1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = 1
      il = i
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = ae(i,j)
      ds = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -an(i,j)*uu(i,j+1) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
             + (pw-pe)*dy &
             + 8.0d0*ds*ar*uu_b(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j)
      DO i = 2,nx_in
        il = i
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        ax(il-1) = aw(i,j)
        ds = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
               + 8.0d0*ds*ar*uu_b(i)/3.0d0
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                 aw(i,j)*uu(i-1,j)
      ENDDO
      DO i = nx_in+1,nx-1
        il = i
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        ax(il-1) = aw(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                 -as(i,j)*uu(i,j-1)
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                 aw(i,j)*uu(i-1,j)
      ENDDO
      i = nx
      il = i
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = aw(i,j)
      ds = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = pres_r(j)
      bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
               -as(i,j)*uu(i,j-1)
      bx(il) = bx(il) - aw(i,j)*uu(i-1,j) - ao(i,j)*uu(i,j)
  
      CALL TRI(nx,ax,ddx,cx,bx,solx)
  
      DO i = 1,nx
        uu(i,j) = uu(i,j) + solx(i)
      ENDDO

!   Interior rows: Top Section
      DO j = ny_in+2,ny-1
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        i = 1
        il = i
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        dw = gam(i,j)
        pw = pres(i,j)
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -an(i,j)*uu(i,j+1) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
               + (pw-pe)*dy &
               -as(i,j)*uu(i,j-1) + 8.0d0*dw*ra*uu_l(j)/3.0d0
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j)
        DO i = 2,nx-1
          il = i
          ddx(il) = ao(i,j)*(one+rin_uv)
          cx(il) = ae(i,j)
          ax(il-1) = aw(i,j)
          pw = half*(pres(i,j)+pres(i-1,j))
          pe = half*(pres(i,j)+pres(i+1,j))
          bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                   -as(i,j)*uu(i,j-1)
          bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                   aw(i,j)*uu(i-1,j)
        ENDDO
        i = nx
        il = i
        ddx(il) = ao(i,j)*(one+rin_uv)
        ax(il-1) = aw(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = pres_r(j)
        bx(il) = -an(i,j)*uu(i,j+1) + (pw-pe)*dy &
                 -as(i,j)*uu(i,j-1)
        bx(il) = bx(il) - aw(i,j)*uu(i-1,j) - ao(i,j)*uu(i,j)
    
        CALL TRI(nx,ax,ddx,cx,bx,solx)
  
        DO i = 1,nx
          uu(i,j) = uu(i,j) + solx(i)
        ENDDO
  
      ENDDO
  
!   Top row
      j = ny
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = 1
      il = i
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = ae(i,j)
      dn = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -as(i,j)*uu(i,j-1) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
             + (pw-pe)*dy &
             + 8.0d0*dn*ar*uu_t(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j)
      DO i = 2,nx-1
        il = i
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = ae(i,j)
        ax(il-1) = aw(i,j)
        dn = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -as(i,j)*uu(i,j-1) + (pw-pe)*dy &
               + 8.0d0*dn*ar*uu_t(i)/3.0d0
        bx(il) = bx(il) - ae(i,j)*uu(i+1,j) - ao(i,j)*uu(i,j) - &
                 aw(i,j)*uu(i-1,j)
      ENDDO
      i = nx
      il = i
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = aw(i,j)
      dn = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = pres_r(j)
      bx(il) = -as(i,j)*uu(i,j-1) + (pw-pe)*dy &
               + 8.0d0*dn*ar*uu_t(i)/3.0d0 
      bx(il) = bx(il) - aw(i,j)*uu(i-1,j) - ao(i,j)*uu(i,j)
  
      CALL TRI(nx,ax,ddx,cx,bx,solx)
  
      DO i = 1,nx
        uu(i,j) = uu(i,j) + solx(i)
      ENDDO
 
! Column-wise Sweep
! Inlet Section
      i = 1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = ny_in+1
      il = j-ny_in
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = an(i,j)
      ds = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -ae(i,j)*uu(i+1,j) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
             + (pw-pe)*dy &
             + 8.0d0*ds*ar*uu_b(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j)
      DO j = ny_in+2,ny-1
        il = j-ny_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ax(il-1) = as(i,j)
        dw = gam(i,j)
        pw = pres(i,j)
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
               + (pw-pe)*dy &
               + 8.0d0*dw*ra*uu_l(j)/3.0d0
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                 as(i,j)*uu(i,j-1)
      ENDDO
      j = ny
      il = j-ny_in
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = as(i,j)
      dn = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -ae(i,j)*uu(i+1,j) + rhof(i,j,2)*uu_l(j)*uu_l(j)*dy &
             + (pw-pe)*dy &
             + 8.0d0*dn*ar*uu_t(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - as(i,j)*uu(i,j-1) - ao(i,j)*uu(i,j)
  
      CALL TRI(ny-ny_in,ax,ddx,cx,bx,solx)
  
      DO j = ny_in+1,ny
        uu(i,j) = uu(i,j) + solx(j-ny_in)
      ENDDO

! Inlet Section: Interior Columns
      DO i = 2,nx_in
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        j = ny_in+1
        il = j-ny_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ds = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
               + 8.0d0*ds*ar*uu_b(i)/3.0d0 
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j)
        DO j = ny_in+2,ny-1
          il = j-ny_in
          ddx(il) = ao(i,j)*(one+rin_uv)
          cx(il) = an(i,j)
          ax(il-1) = as(i,j)
          pw = half*(pres(i,j)+pres(i-1,j))
          pe = half*(pres(i,j)+pres(i+1,j))
          bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy 
          bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                   as(i,j)*uu(i,j-1)
        ENDDO
        j = ny
        il = j-ny_in
        ddx(il) = ao(i,j)*(one+rin_uv)
        ax(il-1) = as(i,j)
        dn = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
               + 8.0d0*dn*ar*uu_t(i)/3.0d0 
        bx(il) = bx(il) - as(i,j)*uu(i,j-1) - ao(i,j)*uu(i,j)
    
        CALL TRI(ny-ny_in,ax,ddx,cx,bx,solx)
  
        DO j = ny_in+1,ny
          uu(i,j) = uu(i,j) + solx(j-ny_in)
        ENDDO
      ENDDO

! Outlet Section: First Column
      i = nx_in+1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = 1
      il = j
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = an(i,j)
      ds = gam(i,j)
      dw = gam(i,j)
      pw = pres(i,j)
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -ae(i,j)*uu(i+1,j) + (pw-pe)*dy &
             + 8.0d0*ds*ar*uu_b(i)/3.0d0 + 8.0d0*dw*ra*uu_l(j)/3.0d0
      bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j)
      DO j = 2,ny_in
        il = j
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ax(il-1) = as(i,j)
        dw = gam(i,j)
        pw = pres(i,j)
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) + (pw-pe)*dy &
               + 8.0d0*dw*ra*uu_l(j)/3.0d0
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                 as(i,j)*uu(i,j-1)
      ENDDO
      DO j = ny_in+1,ny-1
        il = j
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ax(il-1) = as(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy 
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                 as(i,j)*uu(i,j-1)
      ENDDO
      j = ny
      il = j
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = as(i,j)
      dn = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = half*(pres(i,j)+pres(i+1,j))
      bx(il) = -ae(i,j)*uu(i+1,j) -aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
             + 8.0d0*dn*ar*uu_t(i)/3.0d0 
      bx(il) = bx(il) - as(i,j)*uu(i,j-1) - ao(i,j)*uu(i,j)
  
      CALL TRI(ny,ax,ddx,cx,bx,solx)
  
      DO j = 1,ny
        uu(i,j) = uu(i,j) + solx(j)
      ENDDO

! Outlet Section: Interior Columns
      DO i = nx_in+2,nx-1
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        j = 1
        il = j
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ds = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
               + 8.0d0*ds*ar*uu_b(i)/3.0d0 
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j)
        DO j = 2,ny-1
          il = j
          ddx(il) = ao(i,j)*(one+rin_uv)
          cx(il) = an(i,j)
          ax(il-1) = as(i,j)
          pw = half*(pres(i,j)+pres(i-1,j))
          pe = half*(pres(i,j)+pres(i+1,j))
          bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy 
          bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                   as(i,j)*uu(i,j-1)
        ENDDO
        j = ny
        il = j
        ddx(il) = ao(i,j)*(one+rin_uv)
        ax(il-1) = as(i,j)
        dn = gam(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = half*(pres(i,j)+pres(i+1,j))
        bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
               + 8.0d0*dn*ar*uu_t(i)/3.0d0 
        bx(il) = bx(il) - as(i,j)*uu(i,j-1) - ao(i,j)*uu(i,j)
    
        CALL TRI(ny,ax,ddx,cx,bx,solx)
  
        DO j = 1,ny
          uu(i,j) = uu(i,j) + solx(j)
        ENDDO
      ENDDO

! Outlet Section: Rightmost Column
      i = nx
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = 1
      il = j
      ddx(il) = ao(i,j)*(one+rin_uv)
      cx(il) = an(i,j)
      ds = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = pres_r(j)
      bx(il) = - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
             + 8.0d0*ds*ar*uu_b(i)/3.0d0 
      bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j)
      DO j = 2,ny-1
        il = j
        ddx(il) = ao(i,j)*(one+rin_uv)
        cx(il) = an(i,j)
        ax(il-1) = as(i,j)
        pw = half*(pres(i,j)+pres(i-1,j))
        pe = pres_r(j)
        bx(il) = - aw(i,j)*uu(i-1,j) + (pw-pe)*dy 
        bx(il) = bx(il) - an(i,j)*uu(i,j+1) - ao(i,j)*uu(i,j) - &
                 as(i,j)*uu(i,j-1)
      ENDDO
      j = ny
      il = j
      ddx(il) = ao(i,j)*(one+rin_uv)
      ax(il-1) = as(i,j)
      dn = gam(i,j)
      pw = half*(pres(i,j)+pres(i-1,j))
      pe = pres_r(j)
      bx(il) = -ae(i,j)*uu(i+1,j) - aw(i,j)*uu(i-1,j) + (pw-pe)*dy &
             + 8.0d0*dn*ar*uu_t(i)/3.0d0 
      bx(il) = bx(il) - as(i,j)*uu(i,j-1) - ao(i,j)*uu(i,j)
  
      CALL TRI(ny,ax,ddx,cx,bx,solx)

      DO j = 1,ny
        uu(i,j) = uu(i,j) + solx(j)
      ENDDO

    ENDDO ! Inner momentum iterations
  
    END SUBROUTINE solve_x_mom
