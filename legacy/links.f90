!---------------------------------------------------------------------
    SUBROUTINE links(var_id)

! Calculates links for a generic advection-diffusion equation for a 
! co-located mesh finite-volume method

    USE precisions
    USE flow
    USE constants
    USE geometry

    IMPLICIT NONE
    INTEGER(int_p), INTENT(in) :: var_id 
    INTEGER(int_p) :: i,j
    REAL(real_p) :: ce,cw,cn,cs,de,dw,dn,ds,re,rw,rn,rs,ue,uw,un,us, &
                    ve,vw,vn,vs
    REAL(real_p), DIMENSION(nx,ny,4) :: rcp
    REAL(real_p), DIMENSION(nx,ny) :: gam
!---------------------------------------------------------------------

    CALL rho_face  ! Calculate density at faces

    IF(var_id == ivar_u .OR. var_id == ivar_v)THEN
      gam(:,:) = mu(:,:)
      rcp(:,:,:) = rhof(:,:,:)
    ELSE
      gam(:,:) = cond(:,:)
      rcp(:,:,:) = rhof(:,:,:)*cp_const
    ENDIF

! Interior cells
! Inlet Section
    DO j = ny_in+2,ny-1
      DO i = 2,nx_in
        re = rcp(i,j,1)
        rw = rcp(i,j,2)
        rn = rcp(i,j,3)
        rs = rcp(i,j,4)
        ue = uface(i,j,1)
        uw = uface(i,j,2)
        vn = vface(i,j,1)
        vs = vface(i,j,2)
        ce = re*ue
        cw = rw*uw
        cn = rn*vn
        cs = rs*vs
        de = half*(gam(i,j)+gam(i+1,j))
        dw = half*(gam(i,j)+gam(i-1,j))
        dn = half*(gam(i,j)+gam(i,j+1))
        ds = half*(gam(i,j)+gam(i,j-1))
        ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
              + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
              + (de+dw)*ra + (dn+ds)*ar
        aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
        ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
        as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
        an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
      ENDDO
    ENDDO
! Outlet Section, Bottom Part
    DO j = 2,ny_in
      DO i = nx_in+2,nx-1
        re = rcp(i,j,1)
        rw = rcp(i,j,2)
        rn = rcp(i,j,3)
        rs = rcp(i,j,4)
        ue = uface(i,j,1)
        uw = uface(i,j,2)
        vn = vface(i,j,1)
        vs = vface(i,j,2)
        ce = re*ue
        cw = rw*uw
        cn = rn*vn
        cs = rs*vs
        de = half*(gam(i,j)+gam(i+1,j))
        dw = half*(gam(i,j)+gam(i-1,j))
        dn = half*(gam(i,j)+gam(i,j+1))
        ds = half*(gam(i,j)+gam(i,j-1))
        ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
              + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
              + (de+dw)*ra + (dn+ds)*ar
        aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
        ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
        as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
        an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
      ENDDO
    ENDDO
! Outlet Section, Top Part
    DO j = ny_in+1,ny-1
      DO i = nx_in+1,nx-1
        re = rcp(i,j,1)
        rw = rcp(i,j,2)
        rn = rcp(i,j,3)
        rs = rcp(i,j,4)
        ue = uface(i,j,1)
        uw = uface(i,j,2)
        vn = vface(i,j,1)
        vs = vface(i,j,2)
        ce = re*ue
        cw = rw*uw
        cn = rn*vn
        cs = rs*vs
        de = half*(gam(i,j)+gam(i+1,j))
        dw = half*(gam(i,j)+gam(i-1,j))
        dn = half*(gam(i,j)+gam(i,j+1))
        ds = half*(gam(i,j)+gam(i,j-1))
        ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
              + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
              + (de+dw)*ra + (dn+ds)*ar
        aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
        ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
        as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
        an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
      ENDDO
    ENDDO

! Inlet Boundary
    i = 1
    DO j = ny_in+2,ny-1
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = half*(gam(i,j)+gam(i+1,j))
      dw = gam(i,j)
      dn = half*(gam(i,j)+gam(i,j+1))
      ds = half*(gam(i,j)+gam(i,j-1))
      ao(i,j) = half*(ABS(ce)+ce)*dy+ &
            + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
            + (de+3.0d0*dw)*ra + (dn+ds)*ar
      aw(i,j) = zero
      ae(i,j) = -half*(ABS(ce)-ce)*dy - (de+dw/3.0d0)*ra
      as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
      an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
    ENDDO

! Vertical Left Wall Boundary
    i = nx_in+1
    DO j = 2,ny_in
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = half*(gam(i,j)+gam(i+1,j))
      dw = gam(i,j)
      dn = half*(gam(i,j)+gam(i,j+1))
      ds = half*(gam(i,j)+gam(i,j-1))
      ao(i,j) = half*(ABS(ce)+ce)*dy+ &
            + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
            + (de+3.0d0*dw)*ra + (dn+ds)*ar
      aw(i,j) = zero
      ae(i,j) = -half*(ABS(ce)-ce)*dy - (de+dw/3.0d0)*ra
      as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
      an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
    ENDDO

! Bottom Wall Boundary
! Inlet Section
    j = ny_in+1 
    DO i = 2,nx_in
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = half*(gam(i,j)+gam(i+1,j))
      dw = half*(gam(i,j)+gam(i-1,j))
      dn = half*(gam(i,j)+gam(i,j+1))
      ds = gam(i,j)
      ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
            + half*(ABS(cn)+cn)*dx+ &
            + (de+dw)*ra + (dn+3.0d0*ds)*ar
      aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
      ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
      as(i,j) = zero
      an(i,j) = -half*(ABS(cn)-cn)*dx - (dn+ds/3.0d0)*ar
    ENDDO
! Outlet Section
    j = 1
    DO i = nx_in+2,nx-1
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = half*(gam(i,j)+gam(i+1,j))
      dw = half*(gam(i,j)+gam(i-1,j))
      dn = half*(gam(i,j)+gam(i,j+1))
      ds = gam(i,j)
      ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
            + half*(ABS(cn)+cn)*dx+ &
            + (de+dw)*ra + (dn+3.0d0*ds)*ar
      aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
      ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
      as(i,j) = zero
      an(i,j) = -half*(ABS(cn)-cn)*dx - (dn+ds/3.0d0)*ar
    ENDDO

! Top Wall
    j = ny
    DO i = 2,nx-1
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = half*(gam(i,j)+gam(i+1,j))
      dw = half*(gam(i,j)+gam(i-1,j))
      dn = gam(i,j)
      ds = half*(gam(i,j)+gam(i,j-1))
      ao(i,j) = half*(ABS(ce)+ce)*dy+half*(ABS(cw)-cw)*dy &
            + half*(ABS(cs)-cs)*dx &
            + (de+dw)*ra + (3.0d0*dn+ds)*ar
      aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
      ae(i,j) = -half*(ABS(ce)-ce)*dy - de*ra
      as(i,j) = -half*(ABS(cs)+cs)*dx - (ds+dn/3.0d0)*ar
      an(i,j) = zero
    ENDDO

! Outlet Boundary
    i = nx
    DO j = 2,ny-1
      re = rcp(i,j,1)
      rw = rcp(i,j,2)
      rn = rcp(i,j,3)
      rs = rcp(i,j,4)
      ue = uface(i,j,1)
      uw = uface(i,j,2)
      vn = vface(i,j,1)
      vs = vface(i,j,2)
      ce = re*ue
      cw = rw*uw
      cn = rn*vn
      cs = rs*vs
      de = zero
      dw = half*(gam(i,j)+gam(i-1,j))
      dn = half*(gam(i,j)+gam(i,j+1))
      ds = half*(gam(i,j)+gam(i,j-1))
      ao(i,j) = half*(ABS(cw)-cw)*dy + ce*dy &
            + half*(ABS(cn)+cn)*dx+half*(ABS(cs)-cs)*dx &
            + (de+dw)*ra + (dn+ds)*ar
      aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
      ae(i,j) = zero
      as(i,j) = -half*(ABS(cs)+cs)*dx - ds*ar
      an(i,j) = -half*(ABS(cn)-cn)*dx - dn*ar
    ENDDO

! Inlet Bottom
    i = 1
    j = ny_in+1
    re = rcp(i,j,1)
    rw = rcp(i,j,2)
    rn = rcp(i,j,3)
    rs = rcp(i,j,4)
    ue = uface(i,j,1)
    uw = uface(i,j,2)
    vn = vface(i,j,1)
    vs = vface(i,j,2)
    ce = re*ue
    cw = rw*uw
    cn = rn*vn
    cs = rs*vs
    de = half*(gam(i,j)+gam(i+1,j))
    dw = gam(i,j)
    dn = half*(gam(i,j)+gam(i,j+1))
    ds = gam(i,j)
    ao(i,j) = half*(ABS(ce)+ce)*dy+ &
             + half*(ABS(cn)+cn)*dx+ &
             + (de+3.0d0*dw)*ra + (dn+3.0d0*ds)*ar
    aw(i,j) = zero
    ae(i,j) = -half*(ABS(ce)-ce)*dy - (de+dw/3.0d0)*ra
    as(i,j) = zero
    an(i,j) = -half*(ABS(cn)-cn)*dx - (dn+ds/3.0d0)*ar

! Inlet Top
    i = 1
    j = ny
    re = rcp(i,j,1)
    rw = rcp(i,j,2)
    rn = rcp(i,j,3)
    rs = rcp(i,j,4)
    ue = uface(i,j,1)
    uw = uface(i,j,2)
    vn = vface(i,j,1)
    vs = vface(i,j,2)
    ce = re*ue
    cw = rw*uw
    cn = rn*vn
    cs = rs*vs
    de = half*(gam(i,j)+gam(i+1,j))
    dw = gam(i,j)
    dn = gam(i,j)
    ds = half*(gam(i,j)+gam(i,j-1))
    ao(i,j) = half*(ABS(ce)+ce)*dy+ &
             + half*(ABS(cs)-cs)*dx &
             + (de+3.0d0*dw)*ra + (3.0d0*dn+ds)*ar
    aw(i,j) = zero
    ae(i,j) = -half*(ABS(ce)-ce)*dy - (de+dw/3.0d0)*ra
    as(i,j) = -half*(ABS(cs)+cs)*dx - (dn/3.0d0+ds)*ar
    an(i,j) = zero

! Vertical Wall bottom corner
    i = nx_in+1
    j = 1
    re = rcp(i,j,1)
    rw = rcp(i,j,2)
    rn = rcp(i,j,3)
    rs = rcp(i,j,4)
    ue = uface(i,j,1)
    uw = uface(i,j,2)
    vn = vface(i,j,1)
    vs = vface(i,j,2)
    ce = re*ue
    cw = rw*uw
    cn = rn*vn
    cs = rs*vs
    de = half*(gam(i,j)+gam(i+1,j))
    dw = gam(i,j)
    dn = half*(gam(i,j)+gam(i,j+1))
    ds = gam(i,j)
    ao(i,j) = half*(ABS(ce)+ce)*dy+ &
             + half*(ABS(cn)+cn)*dx+ &
             + (de+3.0d0*dw)*ra + (dn+3.0d0*ds)*ar
    aw(i,j) = zero
    ae(i,j) = -half*(ABS(ce)-ce)*dy - (de+dw/3.0d0)*ra
    as(i,j) = zero
    an(i,j) = -half*(ABS(cn)-cn)*dx - (dn+ds/3.0d0)*ar

! Outlet Bottom
    i = nx
    j = 1
    re = rcp(i,j,1)
    rw = rcp(i,j,2)
    rn = rcp(i,j,3)
    rs = rcp(i,j,4)
    ue = uface(i,j,1)
    uw = uface(i,j,2)
    vn = vface(i,j,1)
    vs = vface(i,j,2)
    ce = re*ue
    cw = rw*uw
    cn = rn*vn
    cs = rs*vs
    de = zero
    dw = half*(gam(i,j)+gam(i-1,j))
    dn = half*(gam(i,j)+gam(i,j+1))
    ds = gam(i,j)
    ao(i,j) = half*(ABS(cw)-cw)*dy + ce*dy &
             + half*(ABS(cn)+cn)*dx &
             + (de+dw)*ra + (dn+3.0d0*ds)*ar
    aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
    ae(i,j) = zero
    as(i,j) = zero
    an(i,j) = -half*(ABS(cn)-cn)*dx - (dn+ds/3.0d0)*ar

! Outlet Top
    i = nx
    j = ny
    re = rcp(i,j,1)
    rw = rcp(i,j,2)
    rn = rcp(i,j,3)
    rs = rcp(i,j,4)
    ue = uface(i,j,1)
    uw = uface(i,j,2)
    vn = vface(i,j,1)
    vs = vface(i,j,2)
    ce = re*ue
    cw = rw*uw
    cn = rn*vn
    cs = rs*vs
    de = zero
    dw = half*(gam(i,j)+gam(i-1,j))
    dn = gam(i,j)
    ds = half*(gam(i,j)+gam(i,j-1))
    ao(i,j) = half*(ABS(cw)-cw)*dy + ce*dy &
             + half*(ABS(cs)-cs)*dx &
             + (de+dw)*ra + (3.0d0*dn+ds)*ar
    aw(i,j) = -half*(ABS(cw)+cw)*dy - dw*ra
    ae(i,j) = zero
    as(i,j) = -half*(ABS(cs)+cs)*dx - (dn/3.0d0+ds)*ar
    an(i,j) = zero

    END subroutine links
