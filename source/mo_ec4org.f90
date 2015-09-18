!
MODULE mo_ec4org
!
! These USE statements will replace the INCLUDE statments
! USE Paraorg_Module, ONLY: IE,JE,KE,IEJE,NOZ,KE1,
!&                          MOIE,MOJE,IEJEKE
! USE Param_Module,   ONLY: NOZ
!
! This module should contain all the data that is passed
! through the numerous subroutines used in EC4ORG and should be
! stored in the "long-term" memory.
! In the end, this module replaces the approach of using
! a large number of local variables in EC4ORG. Each subroutine
! called in EC4ORG then simply uses this module and one
! can get rid of the gigantic subroutine headers.
!
! TODO: Some fields are declared 4-dimensional like in ec4org but
! some are declared as 2-dimensional like in progec4? What
! should we do about this?
!
! Here is what "European Standards For Writing and Documenting Exchangeable Fortran 90 Code" says:
! "Implicitly changing the shape of an array when passing it into a subroutine. Although actually
!  forbidden in the standard it was very common practice in FORTRAN 77 to pass 'n' dimensional
!  arrays into a subroutine where they would, say, be treated as a 1 dimensional array. This
!  practice, though banned in Fortran 90, is still possible with external routines for which no
!  Interface block has been supplied. This only works because of assumptions made about how
!  the data is stored: it is therefore unlikely to work on a massively parallel computer. Hence the
!  practice is banned"
!
IMPLICIT NONE
!*--EC4ORG_MODULE21
!
!
! A metadata type for the the field type, see below.
! Contains all important information concerning a field.
!
TYPE field_info
  !
  SEQUENCE
  !
  ! Content of the field
  !
  CHARACTER(len= 64) :: name          ! variable name
  CHARACTER(len= 64) :: units         ! units
  CHARACTER(len=128) :: longname      ! long name
  !
  ! Memory buffer information
  !
  INTEGER            :: gdim(4)       ! global dimensions of variable
  INTEGER            :: dim (4)       ! local dimensions (lon,lat)
  INTEGER            :: dima(4)       ! allocated dimensions (nproma,ngpblks)
  LOGICAL            :: lreg          ! true for dim==dima
  INTEGER            :: ndim          ! rank of variable
  INTEGER            :: klev          ! number of vertical levels
  INTEGER            :: iklev         ! vertical level index
  LOGICAL            :: alloc         ! pointer must be dealloc. by destruct
  INTEGER            :: repr          ! representation (gridpoint,spectral..)
  !
  ! GRIB output information
  !
  INTEGER            :: gribtable     ! gribcode table number
  INTEGER            :: gribcode      ! gribcode number
  INTEGER            :: gribbits      ! number of bits used for GRIB encoding
  INTEGER            :: levelindx     ! HYBRID, ABOVESUR, SURFACE, BELOWSUR,
                                      ! TILES, SOILLEV, ROOTZONES, CANOPY
  !
  ! NETCDF output information
  !
  INTEGER            :: tracidx       ! tracer index
  INTEGER            :: IO_var_indx(4)! index to dimension table
  INTEGER            :: IO_var_id     ! NETCDF id for internal use (restart)
  INTEGER            :: IO_var_stid   ! NETCDF id for internal use (stream)
  CHARACTER(len=128) :: IO_name       ! name specifier for NETCDF file
  CHARACTER(len= 32) :: IO_unit       ! unit specifier for NETCDF file
  !
  ! Output information
  !
  INTEGER            :: gridID
  INTEGER            :: zaxisID
  !
END TYPE field_info
!
!
! The field type. It contains a pointer to the actual data which 
! is still stored in its original location and variable so it is
! still useable for the rest of REMO.
! The field type will be used in the output module.
!
TYPE field
    REAL, POINTER      :: ptr (:,:,:,:)      ! pointer to the field data
    TYPE (field_info)  :: info               ! meta data for this entry
END TYPE field
!
!INCLUDE "parorg.h"
!INCLUDE "param.h"
!
!-------------------------------------------------------------------------------------------------
!
! Here, we define the fields in a slightly different manner than they use to be
! in ec4org.f.
!
! 3D Fields + Time
! (ie,je,ke,3)
!                          \ Field Position
!                           \    | Layer
!                            \   | / Time
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::                                                   &
                                       t, qd, qw, u, v, fi, zso4all, zso4nat, zozact, sothdt,      &  
                                       uvtk, tke, qi
! (ie,je,ke,2)
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET :: ur, vr, qdr, tr, qwr
REAL, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: tmkvmh
!
!-----------------
! 3D Fields
! (ie,je,ke)
!                           \ Field Position
!                            \ | Layer
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET   :: vervel, ozonpl, so4all, so4nat, aclc, aclcac,     &
                                                 rprac 
! (ie,ke,je)
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET   :: ttk, qdtk, tts, qdts, qwts, qits
! (ieje,ke1)
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET   :: emter, trsol
! 3D Fields
! (iejeke)
!                           \ Field Position
!                            \ 
!REAL, ALLOCATABLE, DIMENSION(:), TARGET   :: 
!-----------------
! 2D Fields + Time
! (ieje,3)
!                           \ Field Position
!                            \ | Time
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET  :: ps, rmy, qdb, qdbl, qdbw, qdbi, ts, tb, tg, tgl, tgw,&
                                       tgi, tsech, tslech, tswech, tsiech, wsech, sn, wl, td, tdcl,& 
                                       td3, td4, td5, tsn, wi3, wi4, wi5, wi, wicl
! 2D Fields + Time
! (ie,je,3)
!                           \ Field Position
!                            \   | Time
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET  :: alphabound 
! (ieje,2)
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET  :: tswechr, tsiechr, seaicer, qdblr, sicedr, psr
!-----------------
! 2D Fields
! (ieje)
!
!                            |Field Position
REAL, ALLOCATABLE, DIMENSION(:), TARGET  ::                                                        &
                                       seaice, bla, fib, phi, rla, coslat, sinlat, coslon, sinlon, &
                                       tmcm, tmch, tmchl, tmchw, tmchi, glac, siced, teff, tsurf,  &
                                       tsmax, tsmin, tslin, qres, dsnac, snmel, runoff, drain,     &
                                       varor, srfl, thfl, qhfl, xhfl, rsfc, ssfc, rsfl, ssfl, ahfl,&
                                       ahfs, ahfsl, ahfsw, ahfsi, ahfice, dhft, dhfqw, dhfqs,      &
                                       topmax, az0, az0l, az0w, az0i, aprc, aprl, aprs, evap,      &
                                       evapm, evapl, evapw, evapi, vgrat, forest, ghpbl, beta,     &
                                       wminlok, wmaxlok, cape, albech, albedo, alsol, alsow, alsoi,&
                                       dew2, wsmx, vlt, fao, rgcgn, tlambda, dlambda, porvol, fcap,&
                                       temp2, t2max, t2min, ustar3, ustr, ustrl, ustrw, ustri,     &
                                       u10, vdis, vstr, vstrl, vstrw, vstri, v10, wind10, wimax,   &
                                       vbm10m, ws1, ws2, ws3, ws4, ws5, dzr, dzs, fksat, fmpot,    &
                                       bclapp, vpor, etrans, ebsoil, esnow, eskin, eres, aclcov,   &
                                       alwcvi, qvi, sclf0, sclfs, sraf0, srafs, tclf0, tclfs,      &
                                       traf0, trafs, aclcv, srads, sradsu, srad0, srad0u, trads,   &
                                       tradsu, trad0, ustrgw, vdisgw, vstrgw, fc, bflhs, bflhsl,   &
                                       bflhsw, bflhsi, bflqds, bflqdsl, bflqdsw, bflqdsi, bflus,   &
                                       bflvs, zt2max, zt2min, qdboxs, qwboxs, ekboxs, fhboxs,      &
                                       fiboxs, qivi, zwimax, qiboxs
REAL, ALLOCATABLE, DIMENSION(:), TARGET     :: emtef, trsof
! (ieje*4)
REAL, ALLOCATABLE, DIMENSION(:), TARGET     :: var
! (ieje,12)
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET   :: zalb, zvgr, zvlt
! (ieje)
INTEGER, ALLOCATABLE, DIMENSION(:), TARGET  :: infrl, infrw, infri
LOGICAL, ALLOCATABLE, DIMENSION(:), TARGET  :: loland, losea, loice, loglac, laland
!        
!-----------------
! One value per layer
! (ke)
!                            |Layer
REAL, ALLOCATABLE, DIMENSION(:), TARGET     :: ak, bk, akh, bkh, dak, dbk, a1t, a2t, vvfh, cevapcu,&
                                               sinue, sicq
! (ke1)
REAL, ALLOCATABLE, DIMENSION(:), TARGET     :: cvdaes, cvdael, cvdaeu, cvdaed
! (ke,ke)
REAL, ALLOCATABLE, DIMENSION(:,:), TARGET   :: sistm, sigam, sitau, sivmt, sivmti
!-----------------
! One value per layer + Time
!
!                           \ Layer 
!                            \ | Time
REAL, ALLOCATABLE, DIMENSION(:,:)   :: gcphi, gacphir, acphir, cphi  
!-----------------
! Fixed
REAL,    ALLOCATABLE, DIMENSION(:)  :: trigsi, rz1i, rz2i, trigsj, rz1j, rz2j
INTEGER, ALLOCATABLE, DIMENSION(:)  :: ifaxi,  ifaxj    
!
!-----------------
!
! ADDITIONAL FIELDS FOR NON-HYDROSTATIC
! PINT          - PRESSURE
! DWDT          - CONVERSION BETWEEN HYDROSTATIC AND NON-HYDROSTATIC
! ETAS          - DIAGNOSTIC VERTICAL VELOCITY
! W             - VERTICAL VELOCITY
!
REAL, ALLOCATABLE, DIMENSION(:,:,:,:), TARGET ::  PINT, DWDT, W
REAL, ALLOCATABLE, DIMENSION(:,:,:)  , TARGET ::  ETAS
!
!--------------------------------------------------------------------------------------------------
!
CONTAINS
!
!--------------------------------------------------------------------------------------------------
! This subroutine inits the EC4ORG database
! The routine allocates all fields
! and connect them with metadata information.
!--------------------------------------------------------------------------------------------------
!
SUBROUTINE EC4ORG_INIT
  !
  IMPLICIT NONE
  !
  TYPE(field) :: test_field
  !
  test_field%ptr => t
  !
  !
END SUBROUTINE EC4ORG_INIT
!
!--------------------------------------------------------------------------------------------------
! This subroutine allocates all fields used in ec4org
!--------------------------------------------------------------------------------------------------
SUBROUTINE EC4ORG_ALLOCATE(ie,je,ke,ke1,noz,moje,moie)
implicit none
integer, intent(in) :: ie,je,ke,ke1,noz,moje,moie
integer :: ieje,ieke,iejeke
!
ieje = ie*je
ieke = ie*ke
iejeke = ie*je*ke
!
ALLOCATE(t(ie,je,ke,3), qd(ie,je,ke,3), qw(ie,je,ke,3), u(ie,je,ke,3), v(ie,je,ke,3),                   &
         fi(ie,je,ke,2), vervel(ie,je,ke,1), ps(ieje,3), seaice(ieje), bla(ieje),                      &
         ozonpl(ie,je,noz,1), so4all(ie,je,ke,1), so4nat(ie,je,ke,1), zso4all(ie,je,ke,14),                  &
         zso4nat(ie,je,ke,14), zozact(ie,je,noz,14), alphabound(ie,je,1,3),                            &
         fib(ieje), ak(ke1), bk(ke1), akh(ke), bkh(ke), dak(ke), dbk(ke), phi(ieje),               &
         rla(ieje), coslat(ieje))
ALLOCATE(gcphi(moje,2), gacphir(moje,2))
ALLOCATE(sinlat(ieje), coslon(ieje), sinlon(ieje), rmy(ieje,3), a1t(ke1), a2t(ke1),                &
         acphir(je,2), cphi(je,2), qdb(ieje,3), qdbl(ieje,3), qdbw(ieje,3), qdbi(ieje,3),          &
         ts(ieje,3))
ALLOCATE(tb(ieje,3), tg(ieje,3), tgl(ieje,3), tgw(ieje,3), tgi(ieje,3), sothdt(ie,ke,je,2),         &
         uvtk(ie,ke,je,2), tmkvmh(ie*(ke-1),je,2), ttk(ie,ke,je,1), qdtk(ie,ke,je,1), tts(ie,ke,je,1),       &
         qdts(ie,ke,je,1), qwts(ie,ke,je,1))
ALLOCATE(tmcm(ieje), tmch(ieje), tmchl(ieje), tmchw(ieje), tmchi(ieje), glac(ieje),                &
         siced(ieje), teff(ieje), tsech(ieje,3), tslech(ieje,3), tswech(ieje,3),                   &
         tsiech(ieje,3), wsech(ieje,3), sn(ieje,3))
ALLOCATE(wl(ieje,3), td(ieje,3), tdcl(ieje,3), td3(ieje,3), td4(ieje,3), td5(ieje,3),              &
         tsn(ieje,3), tsurf(ieje), tsmax(ieje), tsmin(ieje))
ALLOCATE(tslin(ieje), qres(ieje), dsnac(ieje), snmel(ieje), runoff(ieje), drain(ieje),             &
         varor(ieje), srfl(ieje), thfl(ieje), qhfl(ieje), xhfl(ieje), rsfc(ieje), ssfc(ieje))
ALLOCATE(rsfl(ieje), ssfl(ieje), ahfl(ieje), ahfs(ieje), ahfsl(ieje), ahfsw(ieje),                 &
         ahfsi(ieje), ahfice(ieje), dhft(ieje), dhfqw(ieje), dhfqs(ieje), topmax(ieje),            &
         az0(ieje), az0l(ieje), az0w(ieje), az0i(ieje))
ALLOCATE(aprc(ieje), aprl(ieje), aprs(ieje), evap(ieje), evapm(ieje), evapl(ieje),                 &
         evapw(ieje), evapi(ieje), aclc(ie,je,ke,1), aclcac(ie,je,ke,1), vgrat(ieje), forest(ieje),      &
         ghpbl(ieje), beta(ieje), wminlok(ieje), wmaxlok(ieje), cape(ieje))
ALLOCATE(albech(ieje), albedo(ieje), alsol(ieje), alsow(ieje), alsoi(ieje), tke(ie,je,ke,3),        &
         dew2(ieje), wsmx(ieje), vlt(ieje), fao(ieje), rgcgn(ieje), tlambda(ieje),                 &
         dlambda(ieje), porvol(ieje), fcap(ieje), wi3(ieje,3), wi4(ieje,3), wi5(ieje,3),           &
         wi(ieje,3), wicl(ieje,3), temp2(ieje))
ALLOCATE(t2max(ieje), t2min(ieje), ustar3(ieje), ustr(ieje), ustrl(ieje), ustrw(ieje),             &
         ustri(ieje), u10(ieje), vdis(ieje), vstr(ieje), vstrl(ieje), vstrw(ieje),                 &
         vstri(ieje), v10(ieje), wind10(ieje), wimax(ieje), vbm10m(ieje))
ALLOCATE(ws1(ieje), ws2(ieje), ws3(ieje), ws4(ieje), ws5(ieje), dzr(ieje), dzs(ieje),              &
         fksat(ieje), fmpot(ieje), bclapp(ieje), vpor(ieje), etrans(ieje), ebsoil(ieje),           &
         esnow(ieje), eskin(ieje), eres(ieje))
ALLOCATE(aclcov(ieje), alwcvi(ieje), qvi(ieje), emter(ie,je,ke1,1), trsol(ie,je,ke1,1),                  &
         emtef(ieje*2), trsof(ieje*2), sclf0(ieje), sclfs(ieje), sraf0(ieje))
ALLOCATE(srafs(ieje), tclf0(ieje), tclfs(ieje), traf0(ieje), trafs(ieje), aclcv(ieje),             &
         srads(ieje), sradsu(ieje), srad0(ieje), srad0u(ieje))
ALLOCATE(trads(ieje), tradsu(ieje), trad0(ieje), ustrgw(ieje), vdisgw(ieje), vstrgw(ieje),         &
         var(ieje*4), fc(ieje), bflhs(ieje), bflhsl(ieje), bflhsw(ieje), bflhsi(ieje))
ALLOCATE(bflqds(ieje), bflqdsl(ieje), bflqdsw(ieje), bflqdsi(ieje), bflus(ieje), bflvs(ieje),      &
         vvfh(ke), cevapcu(ke), cvdaes(ke1), cvdael(ke1), cvdaeu(ke1), cvdaed(ke1),                & 
         sistm(ke,ke))
ALLOCATE(sigam(ke,ke), sitau(ke,ke), sinue(ke), sivmt(ke,ke), sivmti(ke,ke), sicq(ke),             &
         trigsi(5*(moie-1)/2), rz1i(6), rz2i(6), ifaxi(11))
ALLOCATE(trigsj(5*(moje-1)/2), rz1j(6), rz2j(6), ifaxj(11), ur(ie,je,ke,2), vr(ie,je,ke,2),          &
         tr(ie,je,ke,2), qdr(ie,je,ke,2), psr(ie,je,1,2), qwr(ie,je,ke,2), zt2max(ieje), zt2min(ieje),   &
         zwimax(ieje), zalb(ieje,12), zvgr(ieje,12), zvlt(ieje,12))
ALLOCATE(tswechr(ie,je,1,2), tsiechr(ie,je,1,2), seaicer(ie,je,1,2), qdblr(ie,je,1,2), sicedr(ie,je,1,2))
ALLOCATE(qdboxs(ieje), qwboxs(ieje), ekboxs(ieje), fhboxs(ieje), fiboxs(ieje))
ALLOCATE(qiboxs(ieje), qi(ie,je,ke,3), qivi(ieje), qits(ie,ke,je,1))
ALLOCATE(loland(ieje), losea(ieje), loice(ieje), loglac(ieje), laland(ieje))
ALLOCATE(infrl(ieje), infrw(ieje), infri(ieje))
ALLOCATE(rprac(ie,je,ke,1))
ALLOCATE(PINT(IE,JE,KE1,3), DWDT(IE,JE,KE ,3), ETAS(IE,JE,KE1), W(IE,JE,KE1,3))
!
END SUBROUTINE EC4ORG_ALLOCATE
!
!
!
!--------------------------------------------------------------------------------------------------
! This subroutine deallocates all fields in mo_ec4org
!--------------------------------------------------------------------------------------------------
SUBROUTINE EC4ORG_DEALLOCATE
  !
  IMPLICIT NONE
  !
END SUBROUTINE EC4ORG_DEALLOCATE
!
END MODULE mo_ec4org
!
! AKH   (K) = 0.5*(AK(K) + AK(K+1))
! BKH   (K) = 0.5*(BK(K) + BK(K+1))
! DAK   (K) =    - AK(K) + AK(K+1)
! DBK   (K) =    - BK(K) + BK(K+1)
! VERTIKALKOORDINATEN-PARAMETER (AK, BK)
