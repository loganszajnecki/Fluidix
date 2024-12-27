!---------------------------------------------------------------------
    SUBROUTINE solve_pp

! Solves pressure correction equation using the ADI method

    USE precisions
    USE flow
    USE constants
    USE geometry
    USE control_parameters

    IMPLICIT NONE
    INTEGER(int_p) :: i,j,iter,il
    REAL(real_p) :: mdot,dy2,dx2,res
    REAL(real_p) :: aep,awp,anp,asp,aop
    REAL(real_p) , PARAMETER :: rnp = 1.0d0
    REAL(real_p), DIMENSION(nx) :: ax,bx,cx,ddx,solx
!---------------------------------------------------------------------

    dx2 = dx*dx
    dy2 = dy*dy
    pp = zero

    DO iter = 1,iter_pp

! Calculate Residual of pressure correction equation
      res = zero
      DO i=2,nx-1
        DO j=ny_in+2,ny-1
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = aep+awp+anp+asp
          res = res + (aep*pp(i+1,j)+awp*pp(i-1,j) &
                    +  anp*pp(i,j+1)+asp*pp(i,j-1)-aop*pp(i,j)-mdot)**2
        ENDDO
      ENDDO
      DO i=nx_in+2,nx-1
        DO j=2,ny_in
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = aep+awp+anp+asp
          res = res + (aep*pp(i+1,j)+awp*pp(i-1,j) &
                    +  anp*pp(i,j+1)+asp*pp(i,j-1)-aop*pp(i,j)-mdot)**2
        ENDDO
      ENDDO
      res = SQRT(MAX(res,tiny))
      
      IF(iter == 1) res_p = res

! Row-wise sweep

! Bottom row
      j = 1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = nx_in+1
      il = i-nx_in
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(aep+anp)
      ddx(il) = aop
      cx(il) = aep
      bx(il) = -mdot -anp*pp(i,j+1)
      DO i = nx_in+2,nx-1
        il = i-nx_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        aop = -(aep+awp+anp)
        ddx(il) = aop
        cx(il) = aep
        ax(il-1) = awp
        bx(il) = -mdot -anp*pp(i,j+1)
      ENDDO
      i = nx
      il = i-nx_in
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(anp+awp)+half*dy2*rhof(i,j,1)*(one/ao(i,j))
      awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!      aop = -(anp+awp)
      ddx(il) = aop
      ax(il-1) = awp
      bx(il) = -mdot -anp*pp(i,j+1)
  
      CALL TRI(nx-nx_in,ax,ddx,cx,bx,solx)
  
      DO i = nx_in+1,nx
        pp(i,j) = pp(i,j) + rnp*(solx(i-nx_in)-pp(i,j))
      ENDDO
  
! Interior Rows: Bottom Section
      DO j = 2,ny_in
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        i = nx_in+1
        il = i-nx_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+anp+asp)
        ddx(il) = aop
        cx(il) = aep
        bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
        DO i = nx_in+2,nx-1
          il = i-nx_in
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = -(aep+awp+anp+asp)
          ddx(il) = aop
          cx(il) = aep
          ax(il-1) = awp
          bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
        ENDDO
        i = nx
        il = i-nx_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(anp+asp+awp)+half*dy2*rhof(i,j,1)*(one/ao(i,j))
        awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!        aop = -(anp+asp+awp)
        ddx(il) = aop
        ax(il-1) = awp
        bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
    
        CALL TRI(nx-nx_in,ax,ddx,cx,bx,solx)

        DO i = nx_in+1,nx
          pp(i,j) = pp(i,j) + rnp*(solx(i-nx_in)-pp(i,j))
        ENDDO
      ENDDO

! Inlet Section Bottom Row
      j = ny_in+1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = 1
      il = i
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(aep+anp)
      ddx(il) = aop
      cx(il) = aep
      bx(il) = -mdot -anp*pp(i,j+1)
      DO i = 2,nx_in
        il = i
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        aop = -(aep+awp+anp)
        ddx(il) = aop
        cx(il) = aep
        ax(il-1) = awp
        bx(il) = -mdot -anp*pp(i,j+1)
      ENDDO
      DO i = nx_in+1,nx-1
        il = i
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+awp+anp+asp)
        ddx(il) = aop
        cx(il) = aep
        ax(il-1) = awp
        bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
      ENDDO
      i = nx
      il = i
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(anp+asp+awp)+half*dy2*rhof(i,j,1)*(one/ao(i,j))
      awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!      aop = -(anp+asp+awp)
      ddx(il) = aop
      ax(il-1) = awp
      bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
  
      CALL TRI(nx,ax,ddx,cx,bx,solx)
  
      DO i = 1,nx
        pp(i,j) = pp(i,j) + rnp*(solx(i)-pp(i,j))
      ENDDO

! Interior Rows: Top Section
      DO j = ny_in+2,ny-1
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        i = 1
        il = i
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+anp+asp)
        ddx(il) = aop
        cx(il) = aep
        bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
        DO i = 2,nx-1
          il = i
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = -(aep+awp+anp+asp)
          ddx(il) = aop
          cx(il) = aep
          ax(il-1) = awp
          bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
        ENDDO
        i = nx
        il = i
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(anp+asp+awp)+half*dy2*rhof(i,j,1)*(one/ao(i,j))
        awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!        aop = -(awp+anp+asp)
        ddx(il) = aop
        ax(il-1) = awp
        bx(il) = -mdot -anp*pp(i,j+1) -asp*pp(i,j-1)
    
        CALL TRI(nx,ax,ddx,cx,bx,solx)
    
        DO i = 1,nx
          pp(i,j) = pp(i,j) + rnp*(solx(i)-pp(i,j))
        ENDDO
    
      ENDDO
  
! Top Row
      j = ny
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      i = 1
      il = i
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(aep+asp)
      ddx(il) = aop
      cx(il) = aep
      bx(il) = -mdot -asp*pp(i,j-1)
      DO i = 2,nx-1
        il = i
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+awp+asp)
        ddx(il) = aop
        cx(il) = aep
        ax(il-1) = awp
        bx(il) = -mdot -asp*pp(i,j-1)
      ENDDO
      i = nx
      il = i
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(asp+awp)+half*dy2*rhof(i,j,1)*(one/ao(i,j))
      awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!      aop = -(awp+asp)
      ddx(il) = aop
      ax(il-1) = awp
      bx(il) = -mdot -asp*pp(i,j-1)
  
      CALL TRI(nx,ax,ddx,cx,bx,solx)
  
      DO i = 1,nx
        pp(i,j) = pp(i,j) + rnp*(solx(i)-pp(i,j))
      ENDDO

! Column-wise sweep

! Left (Inlet) Column
      i = 1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = ny_in+1
      il = j-ny_in
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(aep+anp)
      ddx(il) = aop
      cx(il) = anp
      bx(il) = -mdot -aep*pp(i+1,j)
      DO j = ny_in+2,ny-1
        il = j-ny_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+asp+anp)
        ddx(il) = aop
        cx(il) = anp
        ax(il-1) = asp
        bx(il) = -mdot -aep*pp(i+1,j)
      ENDDO
      j = ny
      il = j-ny_in
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(asp+aep)
      ddx(il) = aop
      ax(il-1) = asp
      bx(il) = -mdot -aep*pp(i+1,j)
  
      CALL TRI(ny-ny_in,ax,ddx,cx,bx,solx)
  
      DO j = ny_in+1,ny
        pp(i,j) = pp(i,j) + rnp*(solx(j-ny_in)-pp(i,j))
      ENDDO

! Interior (Inlet) Column
      DO i = 2,nx_in
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        j = ny_in+1
        il = j-ny_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        aop = -(aep+awp+anp)
        ddx(il) = aop
        cx(il) = anp
        bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
        DO j = ny_in+2,ny-1
          il = j-ny_in
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = -(aep+awp+asp+anp)
          ddx(il) = aop
          cx(il) = anp
          ax(il-1) = asp
          bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
        ENDDO
        j = ny
        il = j-ny_in
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(asp+aep+awp)
        ddx(il) = aop
        ax(il-1) = asp
        bx(il) = -mdot -aep*pp(i+1,j) - awp*pp(i-1,j)
    
        CALL TRI(ny-ny_in,ax,ddx,cx,bx,solx)
    
        DO j = ny_in+1,ny
          pp(i,j) = pp(i,j) + rnp*(solx(j-ny_in)-pp(i,j))
        ENDDO

      ENDDO

! Left Column : Outlet section
      i = nx_in+1
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = 1
      il = j
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(aep+anp)
      ddx(il) = aop
      cx(il) = anp
      bx(il) = -mdot -aep*pp(i+1,j)
      DO j = 2,ny_in
        il = j
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+asp+anp)
        ddx(il) = aop
        cx(il) = anp
        ax(il-1) = asp
        bx(il) = -mdot -aep*pp(i+1,j)
      ENDDO
      DO j = ny_in+1,ny-1
        il = j
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(aep+awp+asp+anp)
        ddx(il) = aop
        cx(il) = anp
        ax(il-1) = asp
        bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
      ENDDO
      j = ny
      il = j
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(asp+aep+awp)
      ddx(il) = aop
      ax(il-1) = asp
      bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
  
      CALL TRI(ny,ax,ddx,cx,bx,solx)
  
      DO j = 1,ny
        pp(i,j) = pp(i,j) + rnp*(solx(j)-pp(i,j))
      ENDDO

! Interior Columns: Outlet Section
      DO i = nx_in+2,nx-1
        ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
        j = 1
        il = j
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        aop = -(aep+awp+anp)
        ddx(il) = aop
        cx(il) = anp
        bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
        DO j = 2,ny-1
          il = j
          mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
               + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
          aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
          awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
          anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
          asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
          aop = -(aep+awp+asp+anp)
          ddx(il) = aop
          cx(il) = anp
          ax(il-1) = asp
          bx(il) = -mdot -aep*pp(i+1,j) -awp*pp(i-1,j)
        ENDDO
        j = ny
        il = j
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        aep = -half*dy2*rhof(i,j,1)*(one/ao(i,j)+one/ao(i+1,j))
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(asp+aep+awp)
        ddx(il) = aop
        ax(il-1) = asp
        bx(il) = -mdot -aep*pp(i+1,j) - awp*pp(i-1,j)
    
        CALL TRI(ny,ax,ddx,cx,bx,solx)
    
        DO j = 1,ny
          pp(i,j) = pp(i,j) + rnp*(solx(j)-pp(i,j))
        ENDDO

      ENDDO
  
! Rightmost Column: Outlet Section
      i = nx
      ax=zero;bx=zero;cx=zero;ddx=zero;solx=zero
      j = 1
      il = j
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
      aop = -(awp+anp) + half*dy2*rhof(i,j,1)*(one/ao(i,j))
      awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!      aop = -(awp+anp)
      ddx(il) = aop
      cx(il) = anp
      bx(il) = -mdot -awp*pp(i-1,j)
      DO j = 2,ny-1
        il = j
        mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
             + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
        awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
        anp = -half*dx2*rhof(i,j,3)*(one/ao(i,j)+one/ao(i,j+1))
        asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
        aop = -(awp+asp+anp) + half*dy2*rhof(i,j,1)*(one/ao(i,j))
        awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!        aop = -(awp+asp+anp)
        ddx(il) = aop
        cx(il) = anp
        ax(il-1) = asp
        bx(il) = -mdot -awp*pp(i-1,j)
      ENDDO
      j = ny
      il = j
      mdot = (rhof(i,j,1)*uface(i,j,1)-rhof(i,j,2)*uface(i,j,2))*dy &
           + (rhof(i,j,3)*vface(i,j,1)-rhof(i,j,4)*vface(i,j,2))*dx
      awp = -half*dy2*rhof(i,j,2)*(one/ao(i,j)+one/ao(i-1,j))
      asp = -half*dx2*rhof(i,j,4)*(one/ao(i,j)+one/ao(i,j-1))
      aop = -(asp+awp) + half*dy2*rhof(i,j,1)*(one/ao(i,j))
      awp = awp + half*dy2*rhof(i,j,1)*(one/ao(i,j))
!      aop = -(asp+awp)
      ddx(il) = aop
      ax(il-1) = asp
      bx(il) = -mdot - awp*pp(i-1,j)
  
      CALL TRI(ny,ax,ddx,cx,bx,solx)
  
      DO j = 1,ny
        pp(i,j) = pp(i,j) + rnp*(solx(j)-pp(i,j))
      ENDDO

    ENDDO ! Inner iterations

!    DO i=1,nx
!      DO j=1,ny
!        WRITE(15,10) i,j,uu(i,j),vv(i,j),pp(i,j)
!      ENDDO
!    ENDDO
 10 FORMAT(I6,I6,3(1x,E14.6))
  
    END SUBROUTINE solve_pp
