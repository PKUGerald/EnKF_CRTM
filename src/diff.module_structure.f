111c111
<    integer       :: vroi_metar          ! vertical radius of influence for sfcshp   
---
>    integer       :: vroi_metar          ! vertical radius of influence for sfcshp
158,165d157
< !-- use_radiance
<    logical       :: use_radiance        ! .true. : assimilated radiance 
<    integer       :: datathin_radiance   ! 0=all data, 2=1/2 data, 10=1/10 radiance
<    integer       :: hroi_radiance       ! horizontal radius of influence for radiance    
<    integer       :: vroi_radiance       ! vertical radius of influence for radiance  
< 
< 
< 
179c171
<    namelist /metar_obs      / use_metar , datathin_metar , hroi_metar , vroi_metar 
---
>    namelist /metar_obs      / use_metar , datathin_metar , hroi_metar , vroi_metar
187d178
<    namelist /radiance / use_radiance, datathin_radiance, hroi_radiance, vroi_radiance
292,299d282
<    type radiance_data_type
<         integer                                    :: num
<         character(len=12),allocatable,dimension(:) :: platform
<         real, allocatable,dimension(:)             :: lat, lon,ii, jj, tb, err
<         integer, allocatable,dimension(:)          :: ch, hroi, hroi_d
<    end type radiance_data_type
< 
< 
305d287
<         type ( radiance_data_type      )         :: radiance
320,324c302
<        integer, allocatable, dimension(:,:)      :: roi         !! (ob_num,3) : 1=horizontal, 2=vertical, 3=h.(non-Q in radiance)
<        character(len=12), allocatable, dimension(:) :: sat      !! Name of the Satellite
<        integer, allocatable, dimension(:)        :: ch          !! channel of the satellite
< 
< 
---
>        integer, allocatable, dimension(:,:)      :: roi         !! (ob_num,2) : 1=horizontal, 2=vertical
