      SUBROUTINE uv_correct

! Corrects velocities (both cell center and face) after mass 
! conservation has been satisfied

      USE precisions
      USE control_parameters
      USE flow
      USE geometry
      USE constants
      IMPLICIT NONE
      INTEGER(int_p) :: i,j

! Correct Cell Center Velocities
! U-Velocity
      DO i = 2,nx-1
        DO j = ny_in+1,ny
          uu(i,j) = uu(i,j) + relax_uv*half*(pp(i-1,j)-pp(i+1,j))*dy/ao(i,j)
        ENDDO
      ENDDO
      i = 1
      DO j = ny_in+1,ny
        uu(i,j) = uu(i,j) + relax_uv*(pp(i,j)- &
                  half*(pp(i+1,j)+pp(i,j)))*dy/ao(i,j)
      ENDDO
      i = nx
      DO j = ny_in+1,ny
        uu(i,j) = uu(i,j) + relax_uv*(half*(pp(i-1,j)+pp(i,j)))*dy/ao(i,j)
      ENDDO
      DO i = nx_in+2,nx-1
        DO j = 1,ny_in
          uu(i,j) = uu(i,j) + relax_uv*half*(pp(i-1,j)-pp(i+1,j))*dy/ao(i,j)
        ENDDO
      ENDDO
      i = nx_in+1
      DO j = 1,ny_in
        uu(i,j) = uu(i,j) + relax_uv*(pp(i,j)- &
                  half*(pp(i+1,j)+pp(i,j)))*dy/ao(i,j)
      ENDDO
      i = nx
      DO j = 1,ny_in
        uu(i,j) = uu(i,j) + relax_uv*(half*(pp(i-1,j)+pp(i,j)))*dy/ao(i,j)
      ENDDO


! V-Velocity
      DO i = 1,nx_in
        DO j = ny_in+2,ny-1
          vv(i,j) = vv(i,j) + relax_uv*half*(pp(i,j-1)-pp(i,j+1))*dx/ao(i,j)
        ENDDO
      ENDDO
      j = ny_in+1
      DO i = 1,nx_in
        vv(i,j) = vv(i,j) + relax_uv*(pp(i,j)- &
                  half*(pp(i,j+1)+pp(i,j)))*dx/ao(i,j)
      ENDDO
      j = ny
      DO i = 1,nx_in
        vv(i,j) = vv(i,j) + relax_uv*(half*(pp(i,j)+pp(i,j-1))- &
                  pp(i,j))*dx/ao(i,j)
      ENDDO
      DO i = nx_in+1,nx
        DO j = 2,ny-1
          vv(i,j) = vv(i,j) + relax_uv*half*(pp(i,j-1)-pp(i,j+1))*dx/ao(i,j)
        ENDDO
      ENDDO
      j = 1
      DO i = nx_in+1,nx
        vv(i,j) = vv(i,j) + relax_uv*(pp(i,j)- &
                  half*(pp(i,j+1)+pp(i,j)))*dx/ao(i,j)
      ENDDO
      j = ny
      DO i = nx_in+1,nx
        vv(i,j) = vv(i,j) + relax_uv*(half*(pp(i,j)+pp(i,j-1))- &
                  pp(i,j))*dx/ao(i,j)
      ENDDO

! Correct Cell Face Velocity
! U-velocity
      DO i = 2,nx-1
        DO j = ny_in+1,ny
          uface(i,j,1) = uface(i,j,1) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i+1,j))*(pp(i,j)-pp(i+1,j))
          uface(i,j,2) = uface(i,j,2) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i-1,j))*(pp(i-1,j)-pp(i,j))
        ENDDO
      ENDDO
      i = 1
      DO j = ny_in+1,ny
        uface(i,j,1) = uface(i,j,1) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i+1,j))*(pp(i,j)-pp(i+1,j))
        uface(i,j,2) = uface(i,j,2)
      ENDDO
      i = nx
      DO j = ny_in+1,ny
        uface(i,j,1) = uu(i,j)
        uface(i,j,2) = uface(i,j,2) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i-1,j))*(pp(i-1,j)-pp(i,j))
      ENDDO
      DO i = nx_in+2,nx-1
        DO j = 1,ny_in
          uface(i,j,1) = uface(i,j,1) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i+1,j))*(pp(i,j)-pp(i+1,j))
          uface(i,j,2) = uface(i,j,2) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i-1,j))*(pp(i-1,j)-pp(i,j))
        ENDDO
      ENDDO
      i = nx_in+1
      DO j = 1,ny_in
        uface(i,j,1) = uface(i,j,1) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i+1,j))*(pp(i,j)-pp(i+1,j))
        uface(i,j,2) = uface(i,j,2)
      ENDDO
      i = nx
      DO j = 1,ny_in
        uface(i,j,1) = uu(i,j)
        uface(i,j,2) = uface(i,j,2) + relax_uv*half*dy* &
                         (one/ao(i,j)+one/ao(i-1,j))*(pp(i-1,j)-pp(i,j))
      ENDDO

! V-Velocity
      DO i = 1,nx_in
        DO j = ny_in+2,ny-1
          vface(i,j,1) = vface(i,j,1) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j+1))*(pp(i,j)-pp(i,j+1))
          vface(i,j,2) = vface(i,j,2) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j-1))*(pp(i,j-1)-pp(i,j))
        ENDDO
      ENDDO
      j = ny_in+1
      DO i = 1,nx_in
        vface(i,j,1) = vface(i,j,1) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j+1))*(pp(i,j)-pp(i,j+1))
        vface(i,j,2) = vface(i,j,2) 
      ENDDO
      j = ny
      DO i = 1,nx_in
        vface(i,j,1) = vface(i,j,1) 
        vface(i,j,2) = vface(i,j,2) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j-1))*(pp(i,j-1)-pp(i,j))
      ENDDO
      DO i = nx_in+1,nx
        DO j = 2,ny-1
          vface(i,j,1) = vface(i,j,1) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j+1))*(pp(i,j)-pp(i,j+1))
          vface(i,j,2) = vface(i,j,2) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j-1))*(pp(i,j-1)-pp(i,j))
        ENDDO
      ENDDO
      j = 1
      DO i = nx_in+1,nx
        vface(i,j,1) = vface(i,j,1) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j+1))*(pp(i,j)-pp(i,j+1))
        vface(i,j,2) = vface(i,j,2) 
      ENDDO
      j = ny
      DO i = nx_in+1,nx
        vface(i,j,1) = vface(i,j,1) 
        vface(i,j,2) = vface(i,j,2) + relax_uv*half*dx* &
                         (one/ao(i,j)+one/ao(i,j-1))*(pp(i,j-1)-pp(i,j))
      ENDDO

      END SUBROUTINE uv_correct
