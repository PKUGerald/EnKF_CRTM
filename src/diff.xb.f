21,30c21,22
< real, dimension(ix, jx, kx+1)           :: ph, phb
< real, dimension(ix, jx, kx  )           :: pt, qv, qc, qr, pb
< real, dimension(ix, jx      )           :: mu, mub, t2m, th2m, q2m, u10m, v10m, psfc
< real, dimension(ix+1, jx, kx)           :: u
< real, dimension(ix, jx+1, kx)           :: v
< real, dimension(2, 2, kx+1)             :: p
< real, dimension(kx)                     :: pres, ptt, qvt, ht
< real                                    :: mu1, mub1, long, grid_u, grid_v
< integer                                 :: i1, j1, k1, i, j, k, m, ii, jj, kk, obs_ii,obs_jj
< integer                                 :: i_ph, i_phb, i_mu, i_mub, i_pt, i_qv, i_qc, i_qr, i_var, i_u, i_v
---
> real, dimension(ix, jx      )           :: t2m, th2m, q2m, u10m, v10m, psfc
> integer                                 :: i1, j1, i, j, m, ii, jj, obs_ii,obs_jj
33,34c25
< real                                    :: tsfc, u10, v10, t2, q2, th2, ps
< real, dimension(ix, jx      )           :: rough
---
> real                                    :: ps, u10, v10, t2, q2, th2
45,57d35
< call roughness_from_landuse ( 'USGS', times, ix, jx, lu_index, rough ) 
< 
< !- calculate q, pressure profile on (obs_ii, obs_jj)
< i_qv = 0 
< i_qc = 0 
< i_qr = 0 
< i_mu = 0
< i_mub = 0
< i_ph = 0
< i_phb = 0
< i_pt = 0
< i_u = 0
< i_v = 0
65,74d42
<    if ( enkfvar(m) == 'QVAPOR    ' ) i_qv=m
<    if ( enkfvar(m) == 'QCLOUD    ' ) i_qc=m
<    if ( enkfvar(m) == 'QRAIN     ' ) i_qr=m
<    if ( enkfvar(m) == 'MU        ' ) i_mu=m
<    if ( enkfvar(m) == 'MUB       ' ) i_mub=m
<    if ( enkfvar(m) == 'PH        ' ) i_ph=m
<    if ( enkfvar(m) == 'PHB       ' ) i_phb=m
<    if ( enkfvar(m) == 'T         ' ) i_pt=m
<    if ( enkfvar(m) == 'U         ' ) i_u=m
<    if ( enkfvar(m) == 'V         ' ) i_v=m
82,86d49
< if(i_qv>0) qv (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qv )
< if(i_qc>0) qc (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qc )
< if(i_qr>0) qr (i1:i1+1, j1:j1+1, 1:kx) = xa( 1:2, 1:2, 1:kx, i_qr )
< if(i_mu>0)  mu (i1:i1+1, j1:j1+1)   = xa( 1:2, 1:2, 1, i_mu )
< if(i_mub>0) mub (i1:i1+1, j1:j1+1)  = xa( 1:2, 1:2, 1, i_mub )
93,97d55
< if ( i_qv == 0 ) call get_variable3d(inputfile, 'QVAPOR    ', ix, jx, kx,   1, qv )
< if ( i_qc == 0 ) call get_variable3d(inputfile, 'QCLOUD    ', ix, jx, kx,   1, qc )
< if ( i_qr == 0 ) call get_variable3d(inputfile, 'QRAIN     ', ix, jx, kx,   1, qr )
< if ( i_mu == 0 ) call get_variable2d(inputfile, 'MU        ', ix, jx, 1,    mu  )
< if ( i_mub== 0 ) call get_variable2d(inputfile, 'MUB       ', ix, jx, 1,    mub )
105,169d62
< !qv (i1:i1+1, j1:j1+1, 1:kx) = qv (i1:i1+1, j1:j1+1, 1:kx) + qc (i1:i1+1, j1:j1+1, 1:kx) + qr (i1:i1+1, j1:j1+1, 1:kx)
< !
< !!...... get total qv at obs' position from horizontal interpolation
< !qvt(1:kx) = dym*(dx*qv(i1+1,j1,1:kx) + dxm*qv(i1,j1,1:kx)) + dy*(dx*qv(i1+1,j1+1,1:kx) + dxm*qv(i1,j1+1,1:kx))
< !!...... get mu,mub at obs' position from horizontal interpolation
< !mu1 = dym*(dx*mu(i1+1,j1  ) + dxm*mu(i1,j1  )) + dy*(dx*mu(i1+1,j1+1) + dxm*mu(i1,j1+1))
< !mub1 = dym*(dx*mub(i1+1,j1  ) + dxm*mub(i1,j1  )) + dy*(dx*mub(i1+1,j1+1) + dxm*mub(i1,j1+1))
< !!...... calculate pressure profile from qv
< !call cal_press_from_q( kx, znu, znw, qvt, mu1, mub1, p_top, pres )
< !
< !!- calculate t (not theta) and height profile
< !if(i_ph>0) ph (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_ph )
< !if(i_phb>0) phb (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_phb )
< !if(i_pt>0) pt (i1:i1+1, j1:j1+1, 1:3) = xa( 1:2,1:2,1:3,i_pt )
< !if ( i_ph == 0 ) call get_variable3d(inputfile, 'PH        ', ix, jx, kx+1, 1, ph )
< !if ( i_phb== 0 ) call get_variable3d(inputfile, 'PHB       ', ix, jx, kx+1, 1, phb)
< !if ( i_pt == 0 ) call get_variable3d(inputfile, 'T         ', ix, jx, kx,   1, pt )
< !
< !!...... get geopotential height profile around obs%position, then horizontal interpolated to obs's position
< !p(1:2,1:2,1:3) = (ph(i1:i1+1, j1:j1+1, 1:3) + phb(i1:i1+1, j1:j1+1, 1:3))/g
< !p(1,1,1:3) = dym*(dx*p(2,1,1:3) + dxm*p(1,1,1:3)) + dy*(dx*p(2,2,1:3) + dxm*p(1,2,1:3))
< !ht(1:2) = 0.5*(p(1,1,1:2)+p(1,1,2:3))
< !!...... get ptt(theta)  at obs' position from horizontal interpolation
< !ptt(1) = dym*(dx*pt(i1+1,j1,  1) + dxm*pt(i1,j1,  1)) + dy*(dx*pt(i1+1,j1+1,1) + dxm*pt(i1,j1+1,1))
< !ptt(1) = theta_to_temp(ptt(1)+to, pres(1))
< ! 
< !!-------------------------------------------------
< !!- calculate surface p and ground t 
< !   psfc = pres(1) + (1. - znu(1))*mub1
< !   tsfc = ptt(1) + 0.0065*(ht(1)-p(1,1,1))
< !!
< !!-------------------------------------------------
< !!...... get model u and v
< !if(i_u>0) u(i1:i1+2,j1:j1+1,1)=xa(1:3,1:2,1,i_u)
< !if(i_v>0) v(i1:i1+1,j1:j1+2,1)=xa(1:2,1:3,1,i_v)
< !if ( i_u == 0 ) call get_variable3d(inputfile, 'U         ', ix+1, jx, kx, 1, u )
< !if ( i_v == 0 ) call get_variable3d(inputfile, 'V         ', ix, jx+1, kx, 1, v )
< !
< !!...... horizontal interp for U
< !if ( obs_ii-i1 == 0.500 ) then
< !   grid_u = u(i1+1,j1,1)*(j1+1-obs_jj) + u(i1+1,j1+1,1)*(obs_jj-j1)
< !else if ( obs_ii-i1 > 0.500 ) then
< !   grid_u = (j1+1-obs_jj)*( u(i1+1,j1  ,1)*(i1+1.5-obs_ii)+u(i1+2,j1  ,1)*(obs_ii-i1-0.5) ) + &
< !             (obs_jj-j1)  *( u(i1+1,j1+1,1)*(i1+1.5-obs_ii)+u(i1+2,j1+1,1)*(obs_ii-i1-0.5) )
< !else if ( obs_ii-i1 < 0.500 ) then
< !   grid_u = (j1+1-obs_jj)*( u(i1,j1  ,1)*(i1+0.5-obs_ii)+u(i1+1,j1  ,1)*(obs_ii-i1+0.5) ) + &
< !             (obs_jj-j1)  *( u(i1,j1+1,1)*(i1+0.5-obs_ii)+u(i1+1,j1+1,1)*(obs_ii-i1+0.5) )
< !endif
< !
< !!...... horizontal interp for V
< !if ( obs_jj-j1 == 0.500 ) then
< !   grid_v = v(i1,j1+1,1)*(i1+1-obs_ii) + v(i1+1,j1+1,1)*(obs_ii-i1)
< !else if ( obs_jj-j1 > 0.500 ) then
< !   grid_v = (i1+1-obs_ii)*( v(i1  ,j1+1,1)*(j1+1.5-obs_jj) + v(i1  ,j1+2,1)*(obs_jj-j1-0.5) ) + &
< !             (obs_ii-i1)  *( v(i1+1,j1+1,1)*(j1+1.5-obs_jj) + v(i1+1,j1+2,1)*(obs_jj-j1-0.5) )
< !else if ( obs_jj-j1 < 0.500 ) then
< !   grid_v = (i1+1-obs_ii)*( v(i1  ,j1,1)*(j1+0.5-obs_jj) + v(i1  ,j1+1,1)*(obs_jj-j1+0.5) ) + &
< !             (obs_ii-i1)  *( v(i1+1,j1,1)*(j1+0.5-obs_jj) + v(i1+1,j1+1,1)*(obs_jj-j1+0.5) )
< !endif
< 
< !-------------------------------------------------
< !- calculate 10m wind, 2m t and q
< !call sfc_wtq( psfc, tsfc, pres(1), ptt(1), qvt(1), grid_u, grid_v,          &
< !              pres(2), ptt(2), qvt(2), ht(1), rough(i1,j1), xland(i1,j1),   &
< !              u10, v10, t2, q2 )  
177,187d69
< 
< !-------------------------------------------------
< !- Correct surface pressure
< !call da_sfc_pre ( psfcm, psfc, t2, q2, p(1,1,1), obs%sta(iob,1), obs%sta(iob,3), obs%sta(iob,4)/1000.)
< !
< !if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1  ,1), u(i1+1,j1  ,1), u(i1+2,j1  ,1)  
< !if( print_detail > 100 )write(*,'(3x,a,5f)')'model u =', u(i1,j1+1,1), u(i1+1,j1+1,1), u(i1+2,j1+1,1)  
< !if( print_detail > 100 )write(*,'(3x,a,5f)')'grid_u, grid_v, u10, v10 =',grid_u, grid_v, u10, v10
< !if( print_detail > 100 )write(*,'(3x,a,4f)')'station elevation, p, t, q =', obs%sta(iob,1:4)
< !if( print_detail > 100 )write(*,'(3x,a,5f)')'t2, q2, model_terrain =', t2, q2, p(1,1,1)
< !if( print_detail > 100 )write(*,'(3x,a,5f)')'psfc and corrected psfc =', psfc, psfcm
192,193d73
< !     xb = psfcm
< !     xb = psfc
592c472
< 
---
>     
594c474
<   i_ph = 0 ; i_phb = 0 ; i_pb = 0 ; i_qr = 0 ; i_psfc = 0;    
---
>   i_ph = 0 ; i_phb = 0 ; i_pb = 0 ; i_qr = 0 ; i_psfc = 0;
727c607
<      if( i_t <1  ) call get_variable3d(inputfile, 'T_2       ', ix, jx, kx, 1, t ) 
---
>      if( i_t <1  ) call get_variable3d(inputfile, 'T_2       ', ix, jx, kx, 1, t )
739c619
< 
---
>    
1065,1608d944
< 
< 
< 
< 
< 
< 
< !=======================================================================================
< subroutine xb_to_radiance(inputfile,proj,ix,jx,kx,xlong,xlat,landmask,iob_radmin,iob_radmax,xb_tb)
< 
< !---------------------
< ! radiance subroutine calculates brightness temperature for satellite channels
< !---------------------
< 
<   USE constants
<   USE netcdf
<   USE mpi_module
<   USE CRTM_Module
<   use namelist_define
<   use obs_define
<   use wrf_tools
< 
< 
<   implicit none
< 
<   integer, intent(in)                      :: ix, jx, kx
<   integer, intent(out)                     :: iob_radmin,iob_radmax
<   character(len=10), intent(in)            :: inputfile
<   type(proj_info), intent(in)              :: proj                   ! 1st guestmap info
<   real, dimension(obs%num), intent(out)    :: xb_tb
<   real, dimension(ix, jx ), intent(in)     :: xlong
<   real, dimension(ix, jx ), intent(in)     :: xlat
<   real, dimension(ix, jx ), intent(in)     :: landmask
<   integer                                  :: iob,irad
<   real                                     :: obs_ii, obs_jj, dx,dxm,dy,dym
< 
<   CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'ctrm'
<   CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/'
<   CHARACTER(*), PARAMETER :: CRTM_code ='/work/03154/tg824524/tools/EnKF_crtm/code/CRTM/crtm_wrf/'
<   REAL, PARAMETER :: P1000MB=100000.D0
<   REAL, PARAMETER :: R_D=287.D0
<   REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0
<   REAL, PARAMETER :: Re=6378000.0
<   !====================
<   !setup for GOES-ABI
<    REAL, PARAMETER :: sat_h=35780000.0
<    REAL, PARAMETER :: sat_lon=-75.0/180.0*3.14159
<    INTEGER, parameter :: n_ch=3        !for GOES-ABI
<   !====================
< !  INTEGER, intent(in) :: ix = ix  !total number of the x-grid
< !  INTEGER, parameter, intent(in) :: jx = jx  !total number of the y-grid
< !  INTEGER, parameter, intent(in) :: kx = kx        !level range
<   ! Profile dimensions...
<   INTEGER, PARAMETER :: N_PROFILES  = 1 
< !  INTEGER, PARAMETER :: N_LAYERS    = kx
<   INTEGER, PARAMETER :: N_ABSORBERS = 2 
< !  INTEGER, PARAMETER :: N_CLOUDS    = kx*5
<   INTEGER, PARAMETER :: N_AEROSOLS  = 0
<   INTEGER, PARAMETER :: N_SENSORS = 1
<   REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis
< 
<   ! Variables
<   CHARACTER(256) :: Message
<   CHARACTER(256) :: Version
<   CHARACTER(256) :: Sensor_Id
<   CHARACTER(256) :: FILE_NAME
<   CHARACTER(256) :: obstype
<   CHARACTER(12)  :: sat_id
<   INTEGER :: Error_Status
<   INTEGER :: Allocate_Status
<   INTEGER :: n_Channels
<   INTEGER :: l, m, irec, yend, ystart, nyi
<   integer :: ncid,ncrcode
<   character(LEN=16) :: var_name
<   character(LEN=3)  :: file_ens
<   integer :: x,y,tt,v,z,n,reci,ens,n_ec,num_radgrid
<   INTEGER :: ncl,icl,k1,k2
<   integer :: lat_radiance(ix*jx)  ! latitude
<   integer :: lon_radiance(ix*jx) ! longitude
<   real :: lat(ix,jx)   ! in radian
<   real :: lon(ix,jx)   ! in radian
<   real :: p(ix,jx,kx)
<   real :: pb(ix,jx,kx)
<   real :: pres(ix,jx,kx)
<   real :: ph(ix,jx,kx+1)
<   real :: phb(ix,jx,kx+1)
<   real :: delz(kx)
<   real :: t(ix,jx,kx)
<   real :: tk(ix,jx,kx)
<   real :: qvapor(ix,jx,kx)
<   real :: qcloud(ix,jx,kx)
<   real :: qrain(ix,jx,kx)
<   real :: qice(ix,jx,kx)
<   real :: qsnow(ix,jx,kx)
<   real :: qgraup(ix,jx,kx)
<   real :: psfc(ix,jx)
<   real :: hgt(ix,jx)
<   real :: tsk(ix,jx)
< !  real :: landmask(ix,jx)
<   real :: Tbsend(ix,jx,n_ch)
<   real :: Tb(ix,jx,n_ch)
< 
<   ! ============================================================================
<   ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
<   !
<   TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
<   TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
<   TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
<   TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
<   TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
<   TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
<   ! ============================================================================
< 
<   ! ============================================================================
<   ! 1.5. **** make a loop to get the number of satellite-radiance-iob ****
<   !
<   num_radgrid = 0
<   check_cycle:do iob=1,obs%num
<     obstype = obs%type(iob)
<     if ( obstype(1:8) == 'Radiance' ) then
<      if(num_radgrid == 0) then
<       num_radgrid = num_radgrid + 1
<       lon_radiance(num_radgrid) = int(obs%position(iob,1))
<       lat_radiance(num_radgrid) = int(obs%position(iob,2))
<       iob_radmin = iob
<       iob_radmax = iob
<      else
<       iob_radmax = iob
<       do irad = 1,num_radgrid
<       if((lon_radiance(irad)==int(obs%position(iob,1))).and.(lat_radiance(irad)==int(obs%position(iob,2))))cycle check_cycle
<       enddo
<       num_radgrid = num_radgrid + 1
<       lon_radiance(num_radgrid) = int(obs%position(iob,1))
<       lat_radiance(num_radgrid) = int(obs%position(iob,2))
<      endif
<     endif
<   enddo check_cycle
< 
<   ! ============================================================================
<   ! --------------
<   CALL CRTM_Version( Version )
<   !if(my_proc_id==0)  write(*,*) "CRTM ver.",TRIM(Version) 
<   ! Get sensor id from user
<   ! -----------------------
<   !It assumes that all the Radiance data is same sattelite as the first data.
<   Sensor_Id = trim(adjustl(obs%sat(iob_radmin)))
<   
<   ! ============================================================================
<   ! 2. **** INITIALIZE THE CRTM ****
<   !
<   ! 2a. This initializes the CRTM for the sensors
<   !     predefined in the example SENSOR_ID parameter.
<   !     NOTE: The coefficient data file path is hard-
<   !           wired for this example.
<   ! --------------------------------------------------
<   !if(my_proc_id==0) WRITE( *,'(/5x,"Initializing the CRTM...")' )
<   Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
<                             ChannelInfo  , &  ! Output
<                             IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
<                             IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
<                             File_Path='coefficients/')
<   IF ( Error_Status /= SUCCESS ) THEN
<     Message = 'Error initializing CRTM'
<     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<     STOP
<   END IF
< 
<   ! 2b. Determine the total number of channels
<   !     for which the CRTM was initialized
<   ! ------------------------------------------
<   ! Specify channel 14 for GOES-R ABI
<   if (Sensor_Id == 'abi_gr' ) then
<     Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset =(/8,9,10/) )
<     IF ( Error_Status /= SUCCESS ) THEN
<       Message = 'Error initializing CRTM'
<       CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<       STOP
<     END IF
<   endif 
<   n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
<   ! ============================================================================
< 
< 
< 
<   ! ============================================================================
<   ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
<   !
<   ! 3a. Allocate the ARRAYS
<   ! -----------------------
<   ! Note that only those structure arrays with a channel
<   ! dimension are allocated here because we've parameterized
<   ! the number of profiles in the N_PROFILES parameter.
<   !
<   ! Users can make the 
<   ! then the INPUT arrays (Atm, Sfc) will also have to be allocated.
<   ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
<   IF ( Allocate_Status /= 0 ) THEN
<     Message = 'Error allocating structure arrays'
<     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<     STOP
<   END IF
< 
<   ! 3b. Allocate the STRUCTURES
<   ! ---------------------------
<   ! The input FORWARD structure
<   CALL CRTM_Atmosphere_Create( Atm, kx, N_ABSORBERS, kx*5, N_AEROSOLS)
<   IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
<     Message = 'Error allocating CRTM Atmosphere structures'
<     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<     STOP
<   END IF
<   ! ============================================================================
< 
<   ! ============================================================================
<   ! 4. **** ASSIGN INPUT DATA ****
<   !
<   ! Fill the Atm structure array.
<   ! NOTE: This is an example program for illustrative purposes only.
<   !       Typically, one would not assign the data as shown below,
<   !       but rather read it from file
<   !
<   ! 4a1. Loading Atmosphere and Surface input
<   ! --------------------------------
< !  call get_variable2d(inputfile,'XLAT',ix,jx,1,xlat)
< !  call get_variable2d(inputfile,'XLONG',ix,jx,1,xlong)
<   call get_variable3d(inputfile,'P',ix,jx,kx,1,p)
<   call get_variable3d(inputfile,'PB',ix,jx,kx,1,pb)
<   call get_variable3d(inputfile,'PH_2',ix,jx,kx+1,1,ph)
<   call get_variable3d(inputfile,'PHB',ix,jx,kx+1,1,phb)
<   call get_variable3d(inputfile,'T_2',ix,jx,kx,1,t)
<   call get_variable3d(inputfile,'QVAPOR',ix,jx,kx,1,qvapor)
<   call get_variable3d(inputfile,'QCLOUD',ix,jx,kx,1,qcloud)
<   call get_variable3d(inputfile,'QRAIN',ix,jx,kx,1,qrain)
<   call get_variable3d(inputfile,'QICE',ix,jx,kx,1,qice)
<   call get_variable3d(inputfile,'QSNOW',ix,jx,kx,1,qsnow)
<   call get_variable3d(inputfile,'QGRAUP',ix,jx,kx,1,qgraup)
<   call get_variable2d(inputfile,'PSFC',ix,jx,1,psfc)
<   call get_variable2d(inputfile,'TSK',ix,jx,1,tsk)
<   call get_variable2d(inputfile,'HGT',ix,jx,1,hgt)
< !  call get_variable2d(inputfile,'LANDMASK',ix,jx,1,landmask)
<   lat = xlat/180.0*3.14159
<   lon = xlong/180.0*3.14159
<   pres = P + PB
<   tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )
<   where(qvapor.lt.0.0) qvapor=1.0e-8
<   where(qcloud.lt.0.0) qcloud=0.0
<   where(qice.lt.0.0) qice=0.0
<   where(qrain.lt.0.0) qrain=0.0
<   where(qsnow.lt.0.0) qsnow=0.0
<   where(qgraup.lt.0.0) qgraup=0.0
< 
<   ! 4a2. Parallerization with grids
<   ! --------------------------------
<   !--- preparation for the x,y-loop
<   if(mod(num_radgrid,nprocs).eq.0) then
<      nyi=num_radgrid/nprocs
<   else
<      nyi=num_radgrid/nprocs+1
<   endif
<   ystart=my_proc_id*nyi+1
<   yend=min(num_radgrid,(my_proc_id+1)*nyi)
< 
<   do iob = ystart, yend
<      obs_ii=lon_radiance(iob)
<      obs_jj=lat_radiance(iob)
<      x = int( obs_ii )
<      y = int( obs_jj )
< 
<   ! 4a3. Converting WRF data for CRTM structure
<   ! --------------------------------
<   !--- converting the data to CRTM structure
< 
<   !*******************************************************************************
<   ! satellite information
<   !*******************************************************************************
< 
<   sat_dis=sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
<   SCAN_ANGLE=180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))**2))
<   ZENITH_ANGLE=SCAN_ANGLE+180.0/3.14159*acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
< 
<   !*******************************************************************************
<   ! load WRF data into CRTM structures
<   !*******************************************************************************
<   !--- calcurating delz
<    do z=1,kx
<     if(z.eq.1) then
<      delz(z) = (PH(x,y,z+1) + PHB(x,y,z+1)) / 9.806 - hgt(x,y)
<     else
<      delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/2/9.806
<     endif
<    enddo
<   !---Atmospheric Profile
<    atm(1)%Climatology         = TROPICAL
<    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
<    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)
<    atm(1)%Level_Pressure(0) = (pres(x,y,kx)*3.0/2.0 - pres(x,y,kx-1)/2.0)/100.0  ! convert from Pa to hPA
< !   atm(1)%Level_Pressure(0) = 0.05
<    do z=kx,1,-1
<      if(z.eq.1) then
<        atm(1)%Level_Pressure(kx-z+1) = psfc(x,y)/100.0  ! convert from Pa tohPA
<      else
<        atm(1)%Level_Pressure(kx-z+1) = ((pres(x,y,z-1)+pres(x,y,z))/2.0)/100.0  ! convert from Pa to hPA
<      endif
<      atm(1)%Pressure(kx-z+1)       = pres(x,y,z) / 100.0
<      atm(1)%Temperature(kx-z+1)    = tk(x,y,z)
<      atm(1)%Absorber(kx-z+1,1)     = qvapor(x,y,z)*1000.0
<    enddo
<    atm(1)%Absorber(:,2) = 5.0E-02 
<    ! when # of vertical layer is 60
<    ! (/1.26E+00, 5.55E-01, 3.24E-01, 1.07E-01, 7.03E-02, 5.87E-02, 6.15E-02,6.43E-02, 6.99E-02, 7.17E-02,&
<    !   7.27E-02, 7.35E-02, 7.38E-02, 7.41E-02, 7.42E-02, 7.41E-02, 7.35E-02,7.31E-02, 7.27E-02, 7.27E-02,&
<    !   7.27E-02, 7.26E-02, 7.17E-02, 7.05E-02, 6.80E-02, 6.73E-02, 6.73E-02,6.76E-02, 6.72E-02, 6.62E-02,&
<    !   6.51E-02, 6.45E-02, 6.44E-02, 6.46E-02, 6.48E-02, 6.49E-02, 6.46E-02,6.42E-02, 6.38E-02, 6.38E-02,&
<    !   6.42E-02, 6.48E-02, 6.56E-02, 6.64E-02, 6.64E-02, 6.72E-02, 6.84E-02,6.84E-02, 6.84E-02, 6.94E-02,&
<    !   6.94E-02, 6.72E-02, 6.72E-02, 6.72E-02, 6.05E-02, 6.05E-02, 6.05E-02,4.12E-02, 4.12E-02, 4.12E-02/)
<   !---Cloud Profile
<   do z=1,kx*5
<    atm(1)%Cloud(z)%Type = 0
<    atm(1)%Cloud(z)%Effective_Radius = 0.0
<    atm(1)%Cloud(z)%Water_Content = 0.0
<   enddo
<    ncl = 0
<    icl = 0
<    !--calculating # of clouds (cloud and rain)
<    do z=kx,1,-1
<      if(qcloud(x,y,z).gt.0.0) then
<        ncl = ncl + 1
<      endif
<      if(qrain(x,y,z).gt.0.0) then
<        ncl = ncl + 1
<      endif
<      if(qice(x,y,z).gt.0.0) then
<        ncl = ncl + 1
<      endif
<      if(qsnow(x,y,z).gt.0.0) then
<        ncl = ncl + 1
<      endif
<      if(qgraup(x,y,z).gt.0.0) then
<        ncl = ncl + 1
<      endif
<    enddo
<    !--Data for cloud
<    atm(1)%n_Clouds         = ncl
<    IF ( atm(1)%n_Clouds > 0 ) THEN
<    do z=kx,1,-1
<      if(qcloud(x,y,z).gt.0.0) then
<        icl = icl + 1
<        k1 = kx-z+1
<        k2 = kx-z+1
<        atm(1)%Cloud(icl)%Type = WATER_CLOUD
<        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
<        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
<            qcloud(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
<      endif
<    enddo
<    do z=kx,1,-1
<      if(qrain(x,y,z).gt.0.0) then
<        icl = icl + 1
<        k1 = kx-z+1
<        k2 = kx-z+1
<        atm(1)%Cloud(icl)%Type = RAIN_CLOUD
<        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
<        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
<            qrain(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
<      endif
<    enddo
<    do z=kx,1,-1
<      if(qice(x,y,z).gt.0.0) then
<        icl = icl + 1
<        k1 = kx-z+1
<        k2 = kx-z+1
<        atm(1)%Cloud(icl)%Type = ICE_CLOUD
<        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
<        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
<            qice(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
<      endif
<    enddo
<    do z=kx,1,-1
<      if(qsnow(x,y,z).gt.0.0) then
<        icl = icl + 1
<        k1 = kx-z+1
<        k2 = kx-z+1
<        atm(1)%Cloud(icl)%Type = SNOW_CLOUD
<        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
<        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
<            qsnow(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
<      endif
<    enddo
<    do z=kx,1,-1
<      if(qgraup(x,y,z).gt.0.0) then
<        icl = icl + 1
<        k1 = kx-z+1
<        k2 = kx-z+1
<        atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
<        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
<        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
<            qgraup(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
<      endif
<    enddo
<    ENDIF
< 
<   !---Surface data
<    if(landmask(x,y).eq.1.0) then
<     sfc(1)%Water_Coverage = 0.0_fp
<     sfc(1)%Land_Coverage = 1.0_fp
<     sfc(1)%Land_Temperature = tsk(x,y)
<     sfc(1)%Soil_Temperature = tsk(x,y)
<    else
<     sfc(1)%Water_Coverage = 1.0_fp
<     sfc(1)%Land_Coverage = 0.0_fp
<     sfc(1)%Water_Type = 1  ! Sea water
<     sfc(1)%Water_Temperature = tsk(x,y)
<    endif
< 
< 
<   ! 4b. GeometryInfo input
<   ! ----------------------
<   ! All profiles are given the same value
<   !  The Sensor_Scan_Angle is optional.
<   CALL CRTM_Geometry_SetValue( Geometry, &
<                                Sensor_Zenith_Angle = ZENITH_ANGLE, &
<                                Sensor_Scan_Angle   = SCAN_ANGLE )
< 
< 
<   ! 4c. Use the SOI radiative transfer algorithm
<   ! --------------------------------------------
<   Options%RT_Algorithm_ID = RT_SOI
<   ! ============================================================================
< 
<   ! ============================================================================
<   ! 5. **** CALL THE CRTM FORWARD MODEL ****
<   !
<   Error_Status = CRTM_Forward( Atm        , &
<                                Sfc        , &
<                                Geometry   , &
<                                ChannelInfo, &
<                                RTSolution , &
<                                Options = Options )
<   IF ( Error_Status /= SUCCESS ) THEN
<     Message = 'Error in CRTM Forward Model'
<     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<     STOP
<   END IF
<   ! ============================================================================
< 
< 
< 
<   ! ============================================================================
<   ! 6. **** Collecting output ****
<   !
<   ! User should read the user guide or the source code of the routine
<   ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
<   ! select the needed variables for outputs.  These variables are contained
<   ! in the structure RTSolution.
<   !
<   !DO m = 1, N_PROFILES
<   !  WRITE( *,'(//7x,"Profile ",i0," output for ",a )') n, TRIM(Sensor_Id)
<   !  DO l = 1, n_Channels
<   !    WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
<   !    CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
<   !  END DO
<   !END DO
< 
<   !---for file output, edited 2014.9.26
<   do l = 1, n_Channels
<       Tbsend(x,y,l) = real(RTSolution(l,1)%Brightness_Temperature)
<   enddo
<   !WRITE(*,'(7x,"Profile (",i0,", ",i0,") finished Tb =  ",f6.2)')x,y,Tbsend(x,y,2)
< 
<   !--- end of iob(x,y)-loop
<   enddo
< 
<   CALL MPI_Allreduce(Tbsend,Tb,ix*jx*n_ch,MPI_REAL,MPI_SUM,comm,ierr)
< !  if(x==24 .and. y==184) write(*,*) my_proc_id
< 
<   ! ============================================================================
< 
<   ! ============================================================================
<   !6.5  **** writing the output ****
<   !
<   if(my_proc_id==0) then
<   do iob = iob_radmin, iob_radmax
<      obs_ii=obs%position(iob,1)
<      obs_jj=obs%position(iob,2)
<      x = int( obs_ii )
<      y = int( obs_jj )
<      if (Sensor_Id == 'abi_gr' ) then
<          if (obs%ch(iob) .eq. 8) xb_tb(iob) = Tb(x,y,1) !6.19um
<          if (obs%ch(iob) .eq. 9) xb_tb(iob) = Tb(x,y,2) !6.95um
<          if (obs%ch(iob) .eq. 10) xb_tb(iob) = Tb(x,y,3) !7.34um
<          if (obs%ch(iob) .eq. 14) write(*,*)'change channel setting for ch14' !xb_tb(iob) = Tb(x,y,4) !11.2um
<      elseif (Sensor_Id == 'imgr_g13' ) then
<          if (obs%ch(iob) .eq. 3) xb_tb(iob) = Tb(x,y,2) !6.19um
<          if (obs%ch(iob) .eq. 4) xb_tb(iob) = Tb(x,y,3) !11.2um
<      endif
<   enddo
<     !--initializing the Tbsend fields for Bcast
<     Tbsend = 0.0
<   endif
<   if(my_proc_id==0) &
<    WRITE(*,'(a10," Tb=",f6.2,"~",f6.2)')inputfile,minval(xb_tb(iob_radmin:iob_radmax)),maxval(xb_tb(iob_radmin:iob_radmax))
< 
<   ! ============================================================================
<   !  **** initializing all Tb and Tbsend fields ****
<   !
<   Tb = 0.0
<   CALL MPI_BCAST(Tbsend,ix*jx*n_ch,MPI_REAL,0,comm,ierr)
< 
<   ! ============================================================================
<   ! 7. **** DESTROY THE CRTM ****
<   !
< !  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
<   Error_Status = CRTM_Destroy( ChannelInfo )
<   IF ( Error_Status /= SUCCESS ) THEN
<     Message = 'Error destroying CRTM'
<     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
<     STOP
<   END IF
<   ! ============================================================================
< 
<   !  call parallel_finish()
< 
<   ! ============================================================================
<    !---for debug by Minamide
<    !write(*,*) 'lpres',atm(1)%Level_Pressure
<    !write(*,*) 'Pres',atm(1)%Pressure
<    !write(*,*) 'Temp', atm(1)%Temperature
<    !write(*,*) 'H2O', atm(1)%Absorber(:,1)
<    !write(*,*) 'delz',delz
<    !write(*,*) 'hgt',hgt(x,y)
<    !write(*,*) 'ph',ph(x,y,:)
<    !write(*,*) 'phb',phb(x,y,:)
<    !write(*,*) 'qcloud',qcloud(x,y,:)
<    !write(*,*) 'qice',qice(x,y,:)
<    !write(*,*) 'qsnow',qsnow(x,y,:)
<    !write(*,*) 'qrain',qrain(x,y,:)
<    !write(*,*) 'qgraup',qgraup(x,y,:)
<    !do z=1,ncl
<    !write(*,*)
<    !'cloud',atm(1)%Cloud(z)%Type,minval(atm(1)%Cloud(z)%Water_Content),'~',maxval(atm(1)%Cloud(z)%Water_Content)
<    !enddo
<   ! ============================================================================
< 
< end subroutine xb_to_radiance
< 
