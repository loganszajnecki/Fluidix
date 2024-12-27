!---------------------------------------------------------------------
     SUBROUTINE face_velocity

! Computes velocity at face using PWIM

       USE precisions
       USE flow
       USE constants
       USE geometry

       IMPLICIT NONE
       INTEGER(int_p) :: i,j
       REAL(real_p) :: vol
       REAL(real_p) :: dpdxO,dpdxE,dpdxW,dpdx_e,dpdx_w
       REAL(real_p) :: dpdyO,dpdyN,dpdyS,dpdy_n,dpdy_s
!---------------------------------------------------------------------

       vol = dx*dy

! X component
! Top part
       DO i = 3,nx-2
         DO j = ny_in+1,ny
           dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
           dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
           dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
           dpdx_e = (pres(i+1,j)-pres(i,j))/dx
           dpdx_w = (pres(i,j)-pres(i-1,j))/dx
           uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                        + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                        - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
           uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                        + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                        - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
         ENDDO
       ENDDO
       i=2
       DO j = ny_in+1,ny
         dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
         dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
         dpdxW = (half*(pres(i,j)+pres(i-1,j))-pres(i-1,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO
       i=1
       DO j = ny_in+1,ny
         dpdxO = (half*(pres(i+1,j)+pres(i,j))-pres(i,j))/dx
         dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = uu_l(j)
       ENDDO
       i = nx-1
       DO j = ny_in+1,ny
         dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
         dpdxE = (pres_r(j)-half*(pres(i+1,j)+pres(i,j)))/dx
         dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO
       i = nx
       DO j = ny_in+1,ny
         dpdxO = (pres_r(j)-half*(pres(i-1,j)+pres(i,j)))/dx
         dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = uu(i,j)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO

! Bottom part
       DO i = nx_in+3,nx-2
         DO j = 1,ny_in
           dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
           dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
           dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
           dpdx_e = (pres(i+1,j)-pres(i,j))/dx
           dpdx_w = (pres(i,j)-pres(i-1,j))/dx
           uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                        + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                        - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
           uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                        + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                        - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
         ENDDO
       ENDDO
       i=nx_in+2
       DO j = 1,ny_in
         dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
         dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
         dpdxW = (half*(pres(i,j)+pres(i-1,j))-pres(i-1,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO
       i=nx_in+1
       DO j = 1,ny_in
         dpdxO = (half*(pres(i+1,j)+pres(i,j))-pres(i,j))/dx
         dpdxE = half*(pres(i+2,j)-pres(i,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = uu_l(j)
       ENDDO
       i = nx-1
       DO j = 1,ny_in
         dpdxO = half*(pres(i+1,j)-pres(i-1,j))/dx
         dpdxE = (pres_r(j)-half*(pres(i+1,j)+pres(i,j)))/dx
         dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
         dpdx_e = (pres(i+1,j)-pres(i,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = half*(uu(i+1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxE/ao(i+1,j) &
                      - (one/ao(i,j)+one/ao(i+1,j))*dpdx_e)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO
       i = nx
       DO j = 1,ny_in
         dpdxO = (pres_r(j)-half*(pres(i-1,j)+pres(i,j)))/dx
         dpdxW = half*(pres(i,j)-pres(i-2,j))/dx
         dpdx_w = (pres(i,j)-pres(i-1,j))/dx
         uface(i,j,1) = uu(i,j)
         uface(i,j,2) = half*(uu(i-1,j)+uu(i,j)) &
                      + half*vol*(dpdxO/ao(i,j)+dpdxW/ao(i-1,j) &
                      - (one/ao(i,j)+one/ao(i-1,j))*dpdx_w)
       ENDDO

! Y component
! Inlet part
       DO i = 1,nx_in
         DO j = ny_in+3,ny-2
           dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
           dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
           dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
           dpdy_n = (pres(i,j+1)-pres(i,j))/dy
           dpdy_s = (pres(i,j)-pres(i,j-1))/dy
           vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                        + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                        - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
           vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                        + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                        - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
         ENDDO
       ENDDO
       j = ny_in+2
       DO i = 1,nx_in
         dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
         dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
         dpdyS = (half*(pres(i,j)+pres(i,j-1))-pres(i,j-1))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO
       j = ny_in+1
       DO i = 1,nx_in
         dpdyO = (half*(pres(i,j+1)+pres(i,j))-pres(i,j))/dy
         dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = vv_b(i)
       ENDDO
       j = ny-1
       DO i = 1,nx_in
         dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
         dpdyN = (pres(i,j+1)-half*(pres(i,j+1)+pres(i,j)))/dy
         dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO
       j = ny
       DO i = 1,nx_in
         dpdyO = (pres(i,j)-half*(pres(i,j)+pres(i,j-1)))/dy
         dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = vv_t(i)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO

! Outlet part
       DO i = nx_in+1,nx
         DO j = 3,ny-2
           dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
           dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
           dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
           dpdy_n = (pres(i,j+1)-pres(i,j))/dy
           dpdy_s = (pres(i,j)-pres(i,j-1))/dy
           vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                        + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                        - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
           vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                        + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                        - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
         ENDDO
       ENDDO
       j = 2
       DO i = nx_in+1,nx
         dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
         dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
         dpdyS = (half*(pres(i,j)+pres(i,j-1))-pres(i,j-1))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO
       j = 1
       DO i = nx_in+1,nx
         dpdyO = (half*(pres(i,j+1)+pres(i,j))-pres(i,j))/dy
         dpdyN = half*(pres(i,j+2)-pres(i,j))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = vv_b(i)
       ENDDO
       j = ny-1
       DO i = nx_in+1,nx
         dpdyO = half*(pres(i,j+1)-pres(i,j-1))/dy
         dpdyN = (pres(i,j+1)-half*(pres(i,j+1)+pres(i,j)))/dy
         dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
         dpdy_n = (pres(i,j+1)-pres(i,j))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = half*(vv(i,j+1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyN/ao(i,j+1) &
                      - (one/ao(i,j)+one/ao(i,j+1))*dpdy_n)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO
       j = ny
       DO i = nx_in+1,nx
         dpdyO = (pres(i,j)-half*(pres(i,j)+pres(i,j-1)))/dy
         dpdyS = half*(pres(i,j)-pres(i,j-2))/dy
         dpdy_s = (pres(i,j)-pres(i,j-1))/dy
         vface(i,j,1) = vv_t(i)
         vface(i,j,2) = half*(vv(i,j-1)+vv(i,j)) &
                      + half*vol*(dpdyO/ao(i,j)+dpdyS/ao(i,j-1) &
                      - (one/ao(i,j)+one/ao(i,j-1))*dpdy_s)
       ENDDO

!      WRITE(25,*)"Printing From Face Velocity"
!      i = nx_in+1
!      DO j = 1,ny_in
!        WRITE(25,*)i,j,uface(i,j,1),uface(i,j,2)
!      ENDDO
!      DO j = 1,ny_in
!        WRITE(25,*)i,j,vface(i,j,1),vface(i,j,2)
!      ENDDO

     END SUBROUTINE face_velocity
