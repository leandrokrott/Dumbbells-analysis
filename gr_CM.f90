program gr_bordin
  implicit none
  integer, parameter :: snaps = 500        !number of snapshots
  integer, parameter :: Nhis = 250
  real(8), parameter :: time = 2000.0
  real(8), parameter :: rho_0 = 0.025
  integer, parameter :: binx = 2500
  integer, parameter :: biny = 1200
  REAL*8, ALLOCATABLE  :: xc(:,:), yc(:,:), zc(:,:), xcc(:,:), ycc(:,:), zcc(:,:)
  REAL*8, ALLOCATABLE  :: xcm(:,:), ycm(:,:), zcm(:,:)
  REAL*8, ALLOCATABLE  :: xcf(:,:), ycf(:,:), zcf(:,:)
  REAL*8, ALLOCATABLE  :: xp(:,:), yp(:,:), zp(:,:),dr(:,:)
  REAL*8, ALLOCATABLE  :: x(:,:), y(:,:), z(:,:), vx(:,:), vy(:,:), vz(:,:)
  real(8) :: x0,y0, z0, msd, msdvec(snaps), timevec(snaps)
  real(8) :: Lado, zmean, LadoH
  real(8) ::  bin_bx(Binx), bin_rx(Binx), bin_drx(Binx), bin_by(Biny), bin_ry(Biny), bin_dry(Biny), bin_V(Binx),fx, fy !for density histogram
  real(8) ::  histo(Binx), histo_measure(Binx) !for density histogram
  real(8) :: sx, sy, sxx, sxy, a, b, D
  real(8) :: zini, zfina,dz,zmax,dumb6,dumb7, dumb8
  real(8) ddd, TEMP, rough, lambda
  integer ncc(snaps), ncf(snaps),aux
  integer dumb2,  dumb4, nn, N,stepAvg,N2
  integer i,j,k, nc, np, nncc, nncf, zcount,tipo,indice
  real(8) :: avg(Nhis),avg2(Nhis),r(Nhis),avgm(Nhis)
  real(8) :: dx, dy, dzp  ! variáveis locais para o loop

  DOUBLE PRECISION::rr,xr,yr,zr,r2,rrm,xrm,yrm,zrm,r2m

  CHARACTER, allocatable :: type(:,:)
  character(LEN=1024) :: flag, dumb1, dumb3, dumb5, dumb9
  character ::  aaaa
  CHARACTER (LEN=6)    :: ESTADO


  open(10, file= 'dump.lammpstrj')
  open(20, file= 'dump_cm.lammpstrj') ! Arquivo para salvar os centros de massa

  READ(*,*) LadoH
  Lado = 2.D0*LadoH

  N = 2000
  N2 = 1000

  allocate (xcm(snaps,N2), ycm(snaps,N2), zcm(snaps,N2))
  allocate (xc(snaps,N), yc(snaps,N), zc(snaps,N), xcc(snaps,N), ycc(snaps,N), zcc(snaps,N))
  allocate (xcf(snaps,N), ycf(snaps,N), zcf(snaps,N),type(snaps,N))
  allocate (xp(snaps,N), yp(snaps,N), zp(snaps,N), dr(snaps,N))
  allocate (x(snaps,N), y(snaps,N), z(snaps,N), vx(snaps,N), vy(snaps,N), vz(snaps,N))

  dz = 1.0
  zcount = 1

  zfina = zini+dz

  do j=1,snaps
    read(10,*) dumb1
    read(10,*) dumb2
    read(10,*) dumb3
    read(10,*) dumb4
    read(10,*) dumb5
    read(10,*) dumb6
    read(10,*) dumb7
    read(10,*) dumb8
    read(10,*) dumb9
    k = 0
    do i= 1, N
       read(10,*)  indice, tipo, x(j,i), y(j,i), z(j,i), vx(j,i), vy(j,i), vz(j,i)
       ! Preenchendo xc, yc, zc para cada partícula
       xc(j,k+1) = x(j,i)
       yc(j,k+1) = y(j,i)
       zc(j,k+1) = z(j,i)
       k = k+1
    enddo

    ! Escrever cabeçalho do arquivo de centros de massa
    write(20,'(A)') 'ITEM: TIMESTEP'
    write(20,*) dumb2
    write(20,'(A)') 'ITEM: NUMBER OF ATOMS'
    write(20,*) N2
    write(20,'(A)') 'ITEM: BOX BOUNDS pp pp pp'
    write(20,*) -Lado/2.0, Lado/2.0
    write(20,*) -Lado/2.0, Lado/2.0
    write(20,*) -Lado/2.0, Lado/2.0
    write(20,'(A)') 'ITEM: ATOMS id type x y z'

    k = 0
    do i = 1, N, 2
      k = k+1
      ! Corrigir deslocamento periódico em x
      dx = x(j,i+1) - x(j,i)
      dx = dx - Lado * dble(nint(dx / Lado))
      xcm(j,k) = x(j,i) + dx / 2.0d0
      if (xcm(j,k) < -Lado/2.0d0) xcm(j,k) = xcm(j,k) + Lado
      if (xcm(j,k) >  Lado/2.0d0) xcm(j,k) = xcm(j,k) - Lado

      ! Corrigir deslocamento periódico em y
      dy = y(j,i+1) - y(j,i)
      dy = dy - Lado * dble(nint(dy / Lado))
      ycm(j,k) = y(j,i) + dy / 2.0d0
      if (ycm(j,k) < -Lado/2.0d0) ycm(j,k) = ycm(j,k) + Lado
      if (ycm(j,k) >  Lado/2.0d0) ycm(j,k) = ycm(j,k) - Lado

      ! Corrigir deslocamento periódico em z
      dzp = z(j,i+1) - z(j,i)
      dzp = dzp - Lado * dble(nint(dzp / Lado))
      zcm(j,k) = z(j,i) + dzp / 2.0d0
      if (zcm(j,k) < -Lado/2.0d0) zcm(j,k) = zcm(j,k) + Lado
      if (zcm(j,k) >  Lado/2.0d0) zcm(j,k) = zcm(j,k) - Lado

      ! Escrever posições dos centros de massa no arquivo
      write(20,'(I10,1X,I1,3(1X,F14.6))') k, 1, xcm(j,k), ycm(j,k), zcm(j,k)
    enddo
    ncc(j) = k

  enddo

  close(20) ! Fechar arquivo de centros de massa

  call RADIAL_XYZ(N,N2,ncc,snaps,xc,yc,zc,r,xcm,ycm,zcm,avg,avg2,avgm,Lado,flag,Nhis)

end program gr_bordin

SUBROUTINE RADIAL_XYZ(N,N2,ncc,snaps,xc,yc,zc,r,xcm,ycm,zcm,avg,avg2,avgm,Lado,flag,Nhis)
  IMPLICIT NONE

  INTEGER::Nhis,i,j,k,snaps,ig,igm,N1,sumN1,frames, N, N2
  DOUBLE PRECISION::rr,rrm,delg,lbox,pi,xr,yr,zr,r2,xrm,yrm,zrm,r2m,vb,nid,rho,Lado,ddd
  DOUBLE PRECISION,DIMENSION(snaps,N)::xc,yc,zc
  DOUBLE PRECISION,DIMENSION(snaps,N2)::xcm,ycm,zcm
  INTEGER, DIMENSION(snaps,N):: type
  integer ncc(snaps)
  real(8) :: avg(Nhis),avg2(Nhis),avgm(Nhis),r(Nhis), gr(snaps,Nhis), gr2(snaps,Nhis),grm(snaps,Nhis)
  character(len=1024) flag

  delg = dble(Lado/(2.d0*nhis))
  pi=4*ATAN(1.)

  rho = N/((Lado*Lado*Lado))
  gr = 0.0d0
  gr2 = 0.0d0
  grm = 0.0d0
  avg(:)=0.d0
  avg2(:)=0.d0
  avgm(:)=0.d0
! !
  DO k=1,snaps
     DO i=1,N-1
        DO j=i+1,N
!
           xr=xc(k,i)-xc(k,j)
           yr=yc(k,i)-yc(k,j)
           zr=zc(k,i)-zc(k,j)
!
           xr=xr-Lado*(NINT(xr/Lado))
           yr=yr-Lado*(NINT(yr/Lado))
           zr=zr-Lado*(NINT(zr/Lado))
           r2=xr*xr+yr*yr+zr*zr
           rr=SQRT(r2)
           IF(rr<=Lado/2.)THEN
              if ((-1.0)**i == 1.0) then !se i é par cacula sempre
                 ig=ceiling(rr/delg)
                 gr(k,ig)=gr(k,ig)+2.
                 gr2(k,ig)=gr2(k,ig)+2.
              else if ((-1.0)**i == -1.0 .and. j /= i+1) then !se i é impar e j diferente de i+1
                 ig=ceiling(rr/delg)
                 gr(k,ig)=gr(k,ig)+2.
              endif
           END IF
        END DO
     END DO

        DO i=1,N2-1
        DO j=i+1,N2
!
           xrm=xcm(k,i)-xcm(k,j)
           yrm=ycm(k,i)-ycm(k,j)
           zrm=zcm(k,i)-zcm(k,j)
!
           xrm=xrm-Lado*(NINT(xrm/Lado))
           yrm=yrm-Lado*(NINT(yrm/Lado))
           zrm=zrm-Lado*(NINT(zrm/Lado))
           r2m=xrm*xrm+yrm*yrm+zrm*zrm
           rrm=SQRT(r2m)
           IF(rrm<=Lado/2.)THEN
              igm=ceiling(rrm/delg)
              grm(k,igm)=grm(k,igm)+2.
           END IF
        END DO
     END DO
END DO

    DO j=1,Nhis
     DO i=1,snaps
        r(j)=delg*(j+0.5)
        vb=((j+1)**3.-j**3.)*delg**3.0d0
               nid=(4./3.)*pi*vb*rho
        gr(i,j)=gr(i,j)/(nid)
        gr2(i,j)=gr2(i,j)/(nid/2.0)
        grm(i,j)=grm(i,j)/(nid/2.0)
     END DO
  END DO
  !
  DO i=1,Nhis
     DO j=1,snaps
        avg(i)=avg(i)+gr(j,i)
        avg2(i)=avg(i)+gr2(j,i)
        avgm(i)=avgm(i)+grm(j,i)
     END DO
  END DO
!
  OPEN(unit=2,file='rdf.dat',action="write")
  OPEN(unit=222,file='rdf2.dat',action="write")
  OPEN(unit=22,file='rdfcm.dat',action="write")

  DO i=1,Nhis
     WRITE(2,'(2(f17.10,1X))')r(i),avg(i)/snaps/N
     WRITE(222,'(2(f17.10,1X))')r(i),avg2(i)/snaps/N
     WRITE(22,'(2(f17.10,1X))')r(i),avgm(i)/snaps/N2
  END DO


END SUBROUTINE RADIAL_XYZ
