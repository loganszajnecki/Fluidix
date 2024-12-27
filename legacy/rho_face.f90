        SUBROUTINE rho_face

        USE precisions
        USE flow, ONLY: rho,rhof,rho_l,rho_t,rho_b
        USE constants, ONLY: half
        USE geometry

        IMPLICIT NONE
        INTEGER(int_p) :: i,j

! Inlet Section
        DO i = 2,nx_in
          DO j = ny_in+2,ny-1
            rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
            rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
            rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
            rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
          ENDDO
        ENDDO

        i = 1
        DO j = ny_in+2,ny-1
          rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
          rhof(i,j,2) = rho_l(j)
          rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
          rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
        ENDDO

        j = ny_in+1
        DO i = 2,nx_in
          rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
          rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
          rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
          rhof(i,j,4) = rho_b(i)
        ENDDO

        j = ny
        DO i = 2,nx-1
          rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
          rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
          rhof(i,j,3) = rho_t(i)
          rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
        ENDDO

        i = 1
        j = ny_in+1
        rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
        rhof(i,j,2) = rho_l(j)
        rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
        rhof(i,j,4) = rho_b(i)

        i = 1
        j = ny
        rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
        rhof(i,j,2) = rho_l(j)
        rhof(i,j,3) = rho_t(i)
        rhof(i,j,4) = half*(rho(i,j-1)+rho(i,j))

! Outlet Section
        DO i = nx_in+2,nx-1
          DO j = 2,ny_in
            rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
            rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
            rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
            rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
          ENDDO
        ENDDO

        DO i = nx_in+1,nx-1
          DO j = ny_in+1,ny-1
            rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
            rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
            rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
            rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
          ENDDO
        ENDDO

        i = nx_in+1
        DO j = 2,ny_in
          rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
          rhof(i,j,2) = rho_l(j)
          rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
          rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
        ENDDO

        j = 1
        DO i = nx_in+2,nx-1
          rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
          rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
          rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
          rhof(i,j,4) = rho_b(i)
        ENDDO

        i = nx
        DO j = 2,ny-1
          rhof(i,j,1) = rho(i,j)
          rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
          rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
          rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))
        ENDDO

        i = nx_in+1
        j = 1
        rhof(i,j,1) = half*(rho(i+1,j)+rho(i,j))
        rhof(i,j,2) = rho_l(j)
        rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
        rhof(i,j,4) = rho_b(i)

        i = nx
        j = 1
        rhof(i,j,1) = rho(i,j)
        rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
        rhof(i,j,3) = half*(rho(i,j)+rho(i,j+1))
        rhof(i,j,4) = rho_b(i)

        i = nx
        j = ny 
        rhof(i,j,1) = rho(i,j)
        rhof(i,j,2) = half*(rho(i-1,j)+rho(i,j))
        rhof(i,j,3) = rho_t(i)
        rhof(i,j,4) = half*(rho(i,j)+rho(i,j-1))

        END SUBROUTINE rho_face
