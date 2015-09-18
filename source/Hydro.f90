!
! Hydro module
!
module Hydro
   !
   implicit none
   !
   public  :: Hydro_init, Hydro_solve
   !
   private :: hydro1D, fill_guardcells, minmod, &
              interface_flux, cfl_timestep
   !
   integer :: ieke,ke1,iah,jah,ieu,iev,jeh,jev,jehgg,jahgg,ieh
   integer :: na,na2,ne,iaa,iahgg,iau,iav,iea,iehgg,nj,nj2
   integer :: neighbor(4)
   integer :: kfl850,kfl500,kfl300
   real    :: r,rddrm1,wls,eddphi,edadphi,dt2,edg,edfakinf
   real    :: alcnva,dtdeh,ed2dt,eddlam,emrdrd,rdrd,rerd,vbcfl,vbmxv,wcp,wcpr,wlk
   logical :: laistep
   !
   ! CB:AUSGEWAEHLTE PARAMETRISIERUNGSKONSTANTEN
   REAL    ::          B1    , B2W   , B2E  , B3    , B4W  , B4E   ,   &
                       UC1   , UC2   , UCL  , RHDE  ,                    &
                       AKS2  , AKS4  , AKT
   !   These fields are initiated by progec4
   !
   REAL    ::  OM850M(1), OM500M(1), OM300M(1)
   !
   !
   !
contains
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_init
      !
      implicit none
      !
      write(*,*) '----- Hydro_init ------------------'
      write(*,*) 'Hydro is not doing anything yet,   '
      write(*,*) 'because i have to admit, i was lazy'
      write(*,*) '----- Hydro_init done -------------'
      !
   end subroutine Hydro_init
   !
   !----------------------------------------------------------------------------------------------
   !
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine Hydro_solve(dt)
      !
      use Grid, only: ndim
      !use Database
      use mo_ec4org
      !
      implicit none
      !
      real, intent(inout) :: dt
      !
      if(ndim>1) then
         write(*,*) '------------------------'
         write(*,*) 'ERROR in Hydro_solve:'
         write(*,*) 'only 1D possilbe!'
         write(*,*) '------------------------'
         stop
      endif
      !
      call fill_guardcells
      !
      call progexp(                                                        &
        !
        ! WARNING! Here are 3 dummy argument names that
        ! differ from actual parameter names
        !
        !     |        |        |
        !     v        v        v
          OM850M(1), OM500M(1), OM300M(1),                                 &
        !     ^        ^        ^                                            
        !     |        |        |                                            
          AK    , BK    , AKH    , BKH    , DAK    , DBK  , A1T   ,        &
          A2T   , VVFH  , FIB    , FC     , ACPHIR , CPHI , PS    ,        &
          TGL   , TGW   , TGI    , QDBL   , QDBW   , QDBI , BFLHSL,        &
          BFLHSW, BFLHSI, BFLQDSL, BFLQDSW, BFLQDSI, TG   , QDB   ,        &
          BFLHS , BFLQDS, BFLUS  , BFLVS  , U      , V    , T     ,        &
          QD    , QW    , FI     , TMKVMH , TMCHL  , TMCHW,                &
          TMCHI , TMCM  , TMCH   , SOTHDT , TTK    , QDTK , UVTK  ,        &
          TTS   , QDTS  , QWTS   , INFRL  , INFRW  , INFRI, QI    ,        &
          QITS  , PINT  , DWDT   , ETAS )                                  
      !
   end subroutine Hydro_solve
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine hydro1D(dt)
      !
      use Grid
      use Database
      !
      implicit none
      !
      real, intent(in) :: dt
      !
      !
   end subroutine hydro1D
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine fill_guardcells
      !
      use Grid
      use Database
      !
      implicit none
      !
      !
   end subroutine fill_guardcells
   !
   !----------------------------------------------------------------------------------------------
   !
   real function interface_flux(q1,q2,q3,q4,v_face,dt)
      !
      use Grid
      !
      implicit none
      !
      real :: q1,q2,q3,q4
      real :: v_face,dt
      real :: r,phi,flux,theta
      !
      theta = sign(1.e0,v_face)
      !
      if(abs(q3-q2).gt.0.e0) then
         if(v_face.ge.0e0) then
            r = (q2-q1)/(q3-q2)
         else 
            r = (q4-q3)/(q3-q2)
         endif
      else
         r = 0.d0
      endif
      !
      select case(fl)
         !
         case('donor-cell')
            phi = 0.d0
         case('Lax-Wendroff')
            phi = 1.d0
         case('Beam-Warming')
            phi = r
         case('Fromm')
            phi = 0.5e0*(1.e0+r)
         case('minmod')
            phi = minmod(1.e0,r)
         case('superbee')
            phi = max(0.e0,min(1.e0,2.e0*r),min(2.e0,r))
         case default
            phi = 0.e0
      end select
      !
      flux = 0.5d0*v_face*((1.e0+theta)*q2+(1.e0-theta)*q3) +  &
             0.5d0*abs(v_face)*(1.e0-abs(v_face*dt/dx))*phi*(q3-q2)
             !
      interface_flux = flux
      !
   end function interface_flux
   !
   !----------------------------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------------------------
   !
   real function phi(fl,r)
   !
   implicit none
   !
   real,         intent(in) :: r
   character(*), intent(in) :: fl
   !
   reaL :: limiter
   !
   select case(fl)
     !
     case('donor-cell')
       limiter = 0.e0
     case('Lax-Wendroff')
       limiter = 1.e0
     case('Beam-Warming')
       limiter = r
     case('Fromm')
       limiter = 0.5e0*(1.e0+r)
     case('minmod')
       limiter = minmod(1.e0,r)
     case('superbee')
       limiter = max(0.e0,min(1.e0,2.e0*r),min(2.e0,r))
     case('hyperbee')
       limiter = 1.0!hyperbee(r)
     case('MC')
       limiter = max(0.e0,min(0.5e0*(1.e0+r),2.e0,2.e0*r))
     case('van Leer')
       limiter = (r+abs(r))/(1.e0+abs(r))
     case('van Albada 1')
       limiter = (r*r+r)/(r*r+1.e0)
     case('van Albada 2')
       limiter = (2.e0*r)/(r*r+1.e0)
     case default
       limiter = 0.e0
   end select
   !
   phi = limiter
   !
   end function phi
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine compute_slope(a,l,n,nvar,r)
     !
     implicit none
     !
     integer, intent(in) :: n,nvar
     real, dimension(nvar,n), intent(in)  :: a,l
     real, dimension(nvar,n), intent(out) :: r
     !
     integer :: i,ivar
     !
     do ivar=1,nvar
       do i=2,n-1
         if(abs(a(ivar,i)) > 0.0) then
            if(l(ivar,i) > 0.0) then
              r(ivar,i) = a(ivar,i-1)/a(ivar,i)
            else
              r(ivar,i) = a(ivar,i+1)/a(ivar,i)
            endif
         else
            r(ivar,i) = 0.0
         endif
       enddo
       r(ivar,1) = r(ivar,2)
       r(ivar,n) = r(ivar,n-1)
     enddo
     !
   end subroutine compute_slope
   !
   !----------------------------------------------------------------------------------------------
   !
   subroutine roe_average(q,qa,w,n)
      !
      ! This subroutine computes Roe's average qa
      ! of a quantity a using the weights w.
      !
      implicit none
      !
      integer,               intent(in)  :: n
      real,    dimension(n), intent(in)  :: q,w
      real,    dimension(n), intent(out) :: qa
      !
      integer :: i
      !
      do i=2,n-1
        qa(i) = (q(i-1)*w(i-1)+q(i)*w(i))/(w(i-1)+w(i))
      enddo
      !
      qa(1) = qa(2)
      qa(n) = qa(n-1)
      !
   end subroutine roe_average
   !
   !----------------------------------------------------------------------------------------------
   !
   !real function slope_limiter(slope,q1,q2,q3,q4,v)
   !   !
   !   implicit none
   !   !
   !   character(80)    :: slope
   !   real :: q1,q2,q3,q4
   !   real :: v
   !   !
   !   select case(slope)
   !      !
   !   end select
   !   !
   !   slope_limiter = 0.d0
   !   !
   !end function slope_limiter
   real function cfl_timestep(n,nvar,l,dx)
   !
   use RuntimeParameters
   
   implicit none
   !
   integer, intent(in) :: n,nvar
   real, intent(in), dimension(nvar,n) :: l
   real, intent(in)    :: dx
   !
   integer :: i,ivar
   real :: dt,lmax,lmin
   real, dimension(n) :: dti
   !
   lmax = 0.0
   lmin = 1.e20
   !
   do i=1,n
     lmax = maxval(l(:,i))
     lmax = minval(l(:,i))
     if(lmax-lmin > 0.0) then
       dti(i) = dx/(lmax-lmin)
     else
       dti(i) = dtmax
     endif
   enddo 
   !
   dt = 1.e-3!cfl * minval(dti)
   !
   if(dt < dtmin) then
      write(*,*) 'WARNING: cfl timestep is less than minimum timestep'
      write(*,*) 'using dtmin'
      dt = dtmin
   endif
   !
   cfl_timestep = dt
   !
   end function cfl_timestep
   !
   !----------------------------------------------------------------------------------------------
   !
   real function minmod(a,b)
      !
      implicit none
      !
      real :: a,b,c
      !
      if(a*b.gt.0.d0) then
         if (abs(a).lt.abs(b)) then 
             c = a
         else
             c = b
         endif
      else
         c = 0
      endif
      !
      minmod = c 
      !
   end function minmod
   !
   !----------------------------------------------------------------------------------------------
   !
end module Hydro
