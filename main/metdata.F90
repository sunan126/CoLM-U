#include <define.h>

MODULE METDAT

! Initial author: Hua Yuan, /04/2014/

   use precision
   use netcdf
   use timemanager
   use user_specified_forcing
   use omp_lib

   implicit none

   save
   private

 ! local variables
   character(len=256) :: fmetdat                 ! directory of forcing data
   character(len=256) :: fmetnam                 ! filename of forcing data
   character(len=256) :: fname(NVAR)             ! file name of forcing data
   integer  :: deltim                            ! model time step length
   integer  :: fid(NVAR)                         ! file id of nc file
   integer  :: vid(NVAR)                         ! variable id of nc file
   integer  :: land(nlands)                      ! coordinate infor
#ifndef  USE_ERA5LAND_DATA
   real(r4) :: dlats(nlats)                      ! latitude values
   real(r4) :: dlons(nlons)                      ! longitude values
   real(r8) :: rlats(nlats)                      ! latitude values
   real(r8) :: rlons(nlons)                      ! longitude values
#endif
   real(r4), allocatable :: latxy(:,:)           ! latitude values in 2d
   real(r4), allocatable :: lonxy(:,:)           ! longitude values in 2d
   real(r4), allocatable :: metdata(:,:)         ! forcing data
   real(r4), allocatable :: metdatatmp(:,:)      ! forcing data
   real(r4), allocatable :: metdata1d(:)         ! forcing data
   real(r8), allocatable :: avgcos(:,:)          ! time-average of cos(zenith)
   type(timestamp) :: tstamp_LB(NVAR)            ! time stamp of low boundary data
   type(timestamp) :: tstamp_UB(NVAR)            ! time stamp of up boundary data
#ifdef USE_ERA5LAND_DATA
   real(r4), allocatable :: metdata_int(:,:)     ! forcing data 
   real(r4) :: dlats_int(nlats+1)                ! latitude values
   real(r4) :: dlons_int(nlons)                ! latitude values
   real(r4) :: dlats(nlats/5)                      ! latitude values
   real(r4) :: dlons(nlons/5)                      ! longitude values
   real(r8) :: rlats(nlats/5)                      ! latitude values
   real(r8) :: rlons(nlons/5)                      ! longitude values
#endif

 ! external functions
   real(r8), external :: orb_coszen              ! cosine of solar zenith angle

 ! define interface
   INTERFACE metreadLBUB
      MODULE procedure metreadLBUB
   END INTERFACE

   INTERFACE metreadpoint
      MODULE procedure metreadpoint
   END INTERFACE

 ! public
   public NVAR, tstamp_LB, tstamp_UB, tintalgo, avgcos, dlats, dlons, rlats, rlons, vname
   public metinit, metreadLBUB, metreadpoint, metpreprocess, metfinal, sanity

CONTAINS

 ! initialization of forcing data
 ! ------------------------------------------------------------
   SUBROUTINE metinit(fmetdatin, fmetnamin, deltatime)

      implicit none
      character(len=256), intent(in) :: fmetdatin  ! forcing data directory
      character(len=256), intent(in) :: fmetnamin  ! forcing data filename
      real(r8),           intent(in) :: deltatime  ! model time step

    ! get value of fmetdat and deltim
      fmetdat = fmetdatin
      fmetnam = fmetnamin
      deltim  = int(deltatime)

    ! set initial values
      fname(:) = 'NULL'
      fid(:)   = -1
      vid(:)   = -1
      tstamp_LB(:) = timestamp(-1, -1, -1)
      tstamp_UB(:) = timestamp(-1, -1, -1)

    ! allocate memory for forcing data
      allocate(latxy(nlons,nlats))       ! latitudes for each point
      allocate(lonxy(nlons,nlats))       ! longitude for each point
#ifdef USE_ERA5LAND_DATA
      allocate(metdata_int(nlons,nlats+1))   ! forcing data
      allocate(metdata(nlons/5,nlats/5))     ! forcing data
      allocate(metdatatmp(nlons/5,nlats/5))  ! forcing data
      allocate(avgcos(nlons/5,nlats/5))      ! time-average of cos(zenith)
#else
      allocate(metdata(nlons,nlats))     ! forcing data
      allocate(metdatatmp(nlons,nlats))  ! forcing data
      allocate(avgcos(nlons,nlats))      ! time-average of cos(zenith)
#endif
      allocate(metdata1d(nlands))        ! forcing data

   END SUBROUTINE metinit

 ! finalization of forcing data
 ! ------------------------------------------------------------
   SUBROUTINE metfinal

      implicit none
      integer :: i

    ! close POINT data file
#ifdef USE_POINT_DATA
      if (fid(1) > 0) then
         ! close(fid(1))
         if (trim(fprefix)=='NC') then
            call sanity( nf90_close(fid(1)) )
            fid(1) = -1
         else
            fid(1) = -1
         endif
      end if
#endif

    ! close files
      do i = 1, NVAR
         if (fid(i) > 0) then
            call sanity( nf90_close(fid(i)) )
            fid(i) = -1
            vid(i) = -1
         end if
      end do

    ! free memory
      deallocate(latxy)
      deallocate(lonxy)
      deallocate(metdata)
      deallocate(metdatatmp)
      deallocate(metdata1d)
#ifdef USE_ERA5LAND_DATA
      deallocate(metdata_int)
#endif
      deallocate(avgcos)

   END SUBROUTINE metfinal

 ! ------------------------------------------------------------
 ! FUNCTION:
 !    metreadLBUB
 !
 ! PURPOSE:
 !    read lower and upper boundary forcing data
 !
 ! NOTE:
 !    major interface of this module
 ! ------------------------------------------------------------
   SUBROUTINE metreadLBUB(idate, lat_i, lon_i, lat_n, lon_n, forcn_LB, forcn_UB)

      implicit none
      integer,  intent(in) :: idate(3)
      integer,  intent(in) :: lat_i
      integer,  intent(in) :: lon_i
      integer,  intent(in) :: lat_n
      integer,  intent(in) :: lon_n
      real(r8), intent(inout) :: forcn_LB(:,:,:)
      real(r8), intent(inout) :: forcn_UB(:,:,:)

      integer         :: i, year, month, time_i
      type(timestamp) :: mtstamp

      mtstamp = idate

      do i = 1, NVAR

         if (trim(vname(i)) == 'NULL') cycle     ! no data, cycle

       ! lower and upper boundary data already exist, cycle
         if ( .NOT.(tstamp_LB(i)=='NULL') .AND. .NOT.(tstamp_UB(i)=='NULL') .AND. &
                     tstamp_LB(i)<=mtstamp .AND. mtstamp<=tstamp_UB(i) ) then
            cycle
         end if

       ! set lower boundary time stamp and get data
         if (tstamp_LB(i) == 'NULL') then
            call setstampLB(mtstamp, i, year, month, time_i)
            call openmetfile(year, month, i)
            call readvar(i, time_i)
            forcn_LB(:,:,i) = metdata(lon_i:(lon_i+lon_n-1), lat_i:(lat_i+lat_n-1))
         end if

       ! set upper boundary time stamp and get data
         if (tstamp_UB(i) == 'NULL' .OR. tstamp_UB(i) < mtstamp) then
            if ( .NOT. (tstamp_UB(i) == 'NULL') ) then
               forcn_LB(:,:,i) = forcn_UB(:,:,i)
            end if
            call setstampUB(i, year, month, time_i)
          ! when reaching the end of forcing data, always reuse the last time step data
            if (year <= endyr) then
               call openmetfile(year, month, i)
               call readvar(i, time_i)
               forcn_UB(:,:,i) = metdata(lon_i:(lon_i+lon_n-1), lat_i:(lat_i+lat_n-1))
            else
               print *, 'NOTE: reaching the end of forcing data, always reuse the last time step data!'
            end if
            if (i == 7) then  ! calculate time average coszen, for shortwave radiation
               call calavgcos(lat_i, lon_i, lat_n, lon_n)
            end if
         end if

         ! modified by wanyi for fixed year
         if (mtstamp%day == 1 .AND. mtstamp%sec == 0) then
            call setstampLB(mtstamp, i, year, month, time_i)
            call openmetfile(year, month, i)
            call readvar(i, time_i)
            forcn_LB(:,:,i) = metdata(lon_i:(lon_i+lon_n-1), lat_i:(lat_i+lat_n-1))
            ! set upper boundary
            tstamp_UB(i) = tstamp_LB(i) + dtime(i)
            call openmetfile(year, month, i)
            call readvar(i, time_i)
            forcn_UB(:,:,i) = metdata(lon_i:(lon_i+lon_n-1), lat_i:(lat_i+lat_n-1))
            if (i == 7) then  ! calculate time average coszen, for shortwave radiation
               call calavgcos(lat_i, lon_i, lat_n, lon_n)
            end if
         end if

      end do

   END SUBROUTINE metreadLBUB

 ! open nc file and set file id, close 'old' file
 ! ------------------------------------------------------------
   SUBROUTINE openmetfile(year, month, var_i)

      implicit none

      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: var_i

      character(256) :: filename

      filename = trim(fmetdat)//trim(getfilename(year, month, var_i))
      if (fname(var_i) .NE. filename) then

       ! close the old file first
         if (fname(var_i) .NE. 'NULL' .AND. fid(var_i) > 0) then
            call sanity( nf90_close(fid(var_i)) )
         end if

       ! open the new file and set the fid/vid
         fname(var_i) = filename
         if (fname(var_i) .NE. 'NULL' .AND. vname(var_i) .NE. 'NULL') then
            call sanity( nf90_open(path=trim(fname(var_i)), mode=nf90_nowrite, ncid=fid(var_i)) )
            call sanity( nf90_inq_varid(fid(var_i), trim(vname(var_i)), vid(var_i)) )
         end if
      end if

   END SUBROUTINE openmetfile

 ! read forcing data and coordinates information
 ! ------------------------------------------------------------
   SUBROUTINE readvar(var_i, time_i)

      implicit none
      integer, intent(in) :: var_i
      integer, intent(in) :: time_i
      integer             :: latid
      integer             :: lonid
      integer             :: landid
      integer             :: i, lati, loni, j
      real(r8)            :: pi
      logical, save       :: firstcall = .true.

      if (fid(var_i) > 0) then

       ! set coordinates information
         if (firstcall) then

            call sanity( nf90_inq_varid(fid(var_i), latname, latid) )
            call sanity( nf90_inq_varid(fid(var_i), lonname, lonid) )
            if (dim2d) then
               call sanity( nf90_get_var(fid(var_i), latid, latxy) )
               call sanity( nf90_get_var(fid(var_i), lonid, lonxy) )
               dlats = latxy(1,:)
               dlons = lonxy(:,1)
            else
#ifdef USE_ERA5LAND_DATA
                call sanity( nf90_get_var(fid(var_i), latid, dlats_int) )
                call sanity( nf90_get_var(fid(var_i), lonid, dlons_int) )
                DO j = 1,nlats,5
                   dlats((j+4)/5) = sum(dlats_int(j:(j+4)))/5.0
                END do
                DO j = 1,nlons,5
                   dlons((j+4)/5) = sum(dlons_int(j:(j+4)))/5.0
                ENDDO
#else
                call sanity( nf90_get_var(fid(var_i), latid, dlats) )
                call sanity( nf90_get_var(fid(var_i), lonid, dlons) )
#endif
               
            end if
#ifdef USE_ERA5LAND_DATA
            if (latrev) dlats = dlats(360:1:-1)
#else
            if (latrev) dlats = dlats(nlats:1:-1) ! reverse in latitude for axis
#endif
            if (lonadj) dlons = dlons - 180.      ! adjust in longitude for axis

            pi = 4._r8 * atan(1._r8)  ! get PI
            rlons = dlons/180.*pi     ! degree to radian
            rlats = dlats/180.*pi     ! ..

            if (.NOT. data2d) then
               call sanity( nf90_inq_varid(fid(var_i), "land", landid) )
               call sanity( nf90_get_var(fid(var_i), landid, land) )
            end if

            firstcall = .false.
         end if

       ! read forcing data
         if (data2d) then
            if (hightdim) then
               call sanity( nf90_get_var(fid(var_i), vid(var_i), metdata(:,:), &
                                         start=(/1,1,1,time_i/), count=(/nlons,nlats,1,1/)) )
            else
#ifdef USE_ERA5LAND_DATA
               call sanity( nf90_get_var(fid(var_i), vid(var_i), metdata_int(:,:), &
                                         start=(/1,1,time_i/), count=(/nlons,nlats+1,1/)) )
               ! do j = 1,nlats
               !    metdata(:,j) = (metdata_int(:,j) + metdata_int(:,j+9))/10.0
               ! END DO
               IF (var_i==1) THEN
                  ! >-100
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<0) metdata_int(j:(j+4),i:(i+4))=280. 
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==2) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<0) metdata_int(j:(j+4),i:(i+4))=0.01
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==3) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<0) metdata_int(j:(j+4),i:(i+4))=96000.
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==4) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<0) metdata_int(j:(j+4),i:(i+4))=0
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==5) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<-300) metdata_int(j:(j+4),i:(i+4))=2
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==6) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<-300) metdata_int(j:(j+4),i:(i+4))=2
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==7) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<0) metdata_int(j:(j+4),i:(i+4))=0
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ELSE IF (var_i==8) THEN
                  DO i=1,nlats,5
                     DO j=1,nlons,5
                        WHERE(metdata_int(j:(j+4),i:(i+4))<-300) metdata_int(j:(j+4),i:(i+4))=300
                        metdata((j+4)/5,(i+4)/5) = sum(metdata_int(j:(j+4),i:(i+4)))/25
                     ENDDO
                  ENDDO
               ENDIF
#else
               call sanity( nf90_get_var(fid(var_i), vid(var_i), metdata(:,:), &
                                         start=(/1,1,time_i/), count=(/nlons,nlats,1/)) )
#endif
            end if
         else
            call sanity( nf90_get_var(fid(var_i), vid(var_i), metdata1d(:), &
                                      start=(/1,time_i/), count=(/nlands,1/)) )
            metdata(:,:) = 0.
            do i = 1, nlands
               lati = land(i)/nlons + 1
               loni = mod(land(i), nlons)
               metdata(loni, lati) = metdata1d(i)
            end do
         end if

       ! data adjustment
         if (latrev) metdata(:,:) = metdata(:,nlats:1:-1) ! reverse in latitude for data
         if (lonadj) then
            metdatatmp(:,:) = metdata(:,:)
#ifdef USE_ERA5LAND_DATA
            metdata(1:nlons/5/2,:) = metdatatmp((nlons/5/2+1):nlons/5,:) ! adjust in longitude for data
            metdata((nlons/5/2+1):nlons/5,:) = metdatatmp(1:nlons/5/2,:) ! 0.5E~359.5E -> 179.5W~179.5E
#else
            metdata(1:nlons/2,:) = metdatatmp((nlons/2+1):nlons,:) ! adjust in longitude for data
            metdata((nlons/2+1):nlons,:) = metdatatmp(1:nlons/2,:) ! 0.5E~359.5E -> 179.5W~179.5E
#endif
         end if
      end if

   END SUBROUTINE readvar

 ! ------------------------------------------------------------
 ! FUNCTION:
 !    setstampLB
 !
 ! PURPOSE:
 !    set the lower boundary time stamp and record information
 !
 ! NOTE:
 !    KEY function of this module
 !
 ! - for time stamp, set it regularly as the model time step.
 ! - for record information, account for:
 !    o year alternation
 !    o month alternation
 !    o leap year
 !    o required dada just beyond the first record
 ! ------------------------------------------------------------
   SUBROUTINE setstampLB(mtstamp, var_i, year, month, time_i)

      implicit none
      type(timestamp), intent(in)  :: mtstamp
      integer,         intent(in)  :: var_i
      integer,         intent(out) :: year
      integer,         intent(out) :: month
      integer,         intent(out) :: time_i

      integer :: i, mday, day, sec
      integer :: months(0:12)

      year = mtstamp%year
      day  = mtstamp%day
      sec  = mtstamp%sec

      tstamp_LB(var_i)%year = year
      tstamp_LB(var_i)%day  = day

    ! in the case of one year one file
      if ( trim(groupby) == 'year' ) then

       ! calculate the intitial second
         sec    = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)-0.01) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(day-1)
         tstamp_LB(var_i)%sec = sec

       ! set time stamp (ststamp_LB)
         if (sec <= 0) then
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            if (tstamp_LB(var_i)%day == 0) then
               tstamp_LB(var_i)%year = year - 1
               if ( isleapyear(tstamp_LB(var_i)%year) ) then
                  tstamp_LB(var_i)%day = 366
               else
                  tstamp_LB(var_i)%day = 365
               end if
            end if
         end if

       ! set record info (year, time_i)
         if ( sec<0 .OR. (sec==0 .AND. offset(var_i).NE.0) ) then

          ! if the required dada just behind the first record
          ! -> set to the first record
            if ( year==startyr .AND. month==startmo .AND. day==1 ) then
               sec = offset(var_i)

          ! else, set to one record backward
            else
               sec = 86400 + sec
               day = day - 1
               if (day == 0) then
                  year = year - 1
                  if ( isleapyear(year) .AND. leapyear) then
                     day = 366
                  else
                     day = 365
                  end if
               end if
            end if
         end if ! end if (sec <= 0)

       ! in case of leapyear with a non-leayyear calendar
       ! use the data 1 day before after FEB 28th (Julian day 59).
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. day>59 ) then
            day = day - 1
         end if

       ! get record time index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

    ! in the case of one month one file
      if ( trim(groupby) == 'month' ) then

         if ( isleapyear(year) ) then
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         else
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         end if

       ! calculate initial month and day values
         call julian2monthday(year, day, month, mday)

       ! calculate initial second value
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)-0.01) *1. / dtime(var_i) ) + 1
         sec    = (time_i-1)*dtime(var_i) + offset(var_i) - 86400*(mday-1)
         tstamp_LB(var_i)%sec  = sec

       ! set time stamp (ststamp_LB)
         if (sec <= 0) then
            tstamp_LB(var_i)%sec = 86400 + sec
            tstamp_LB(var_i)%day = day - 1
            if (tstamp_LB(var_i)%day == 0) then
               tstamp_LB(var_i)%year = year - 1
               if ( isleapyear(tstamp_LB(var_i)%year) ) then
                  tstamp_LB(var_i)%day = 366
               else
                  tstamp_LB(var_i)%day = 365
               end if
            end if
         end if

       ! set record info (year, month, time_i)
         if ( sec<0 .OR. (sec==0 .AND. offset(var_i).NE.0) ) then

          ! if just behind the first record -> set to first record
            if ( year==startyr .AND. month==startmo .AND. mday==1 ) then
               sec = offset(var_i)

          ! set to one record backward
            else
               sec = 86400 + sec
               mday = mday - 1
               if (mday == 0) then
                  month = month - 1
                ! bug found by Zhu Siguang & Zhang Xiangxiang, 05/19/2014
                ! move the below line in the 'else' statement
                  !mday = months(month) - months(month-1)
                  if (month == 0) then
                     month = 12
                     year = year - 1
                     mday = 31
                  else
                     mday = months(month) - months(month-1)
                  end if
               end if
            end if
         end if

       ! in case of leapyear with a non-leayyear calendar
       ! use the data 1 day before, i.e., FEB 28th.
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
            mday = 28
         end if

       ! get record time index
         sec = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

      if (time_i <= 0) then
         write(6, *) "got the wrong time record of forcing! stop!"; stop
      end if

      return

   END SUBROUTINE setstampLB

 ! ------------------------------------------------------------
 ! FUNCTION:
 !    setstampUB
 !
 ! PURPOSE:
 !    set the upper boundary time stamp and record information
 !
 ! NOTE:
 !    KEY function of this module
 ! ------------------------------------------------------------
   SUBROUTINE setstampUB(var_i, year, month, time_i)

      implicit none
      integer,         intent(in)  :: var_i
      integer,         intent(out) :: year
      integer,         intent(out) :: month
      integer,         intent(out) :: time_i

      integer :: mday, day, sec
      integer :: months(0:12)

    ! calculate the time stamp
      if ( tstamp_UB(var_i) == 'NULL' ) then
         tstamp_UB(var_i) = tstamp_LB(var_i) + dtime(var_i)
      else
         tstamp_LB(var_i) = tstamp_UB(var_i)
         tstamp_UB(var_i) = tstamp_UB(var_i) + dtime(var_i)
      end if

    ! calcualte initial year, day, and second values
      year = tstamp_UB(var_i)%year
      day  = tstamp_UB(var_i)%day
      sec  = tstamp_UB(var_i)%sec

      if ( trim(groupby) == 'year' ) then

       ! adjust year value
         if ( sec==86400 .AND. offset(var_i).EQ.0 ) then
            sec = 0
            day = day + 1
            if( isleapyear(year) .AND. day==367) then
               year = year + 1; day = 1
            end if
            if( .NOT. isleapyear(year) .AND. day==366) then
               year = year + 1; day = 1
            end if
         end if

       ! in case of leapyear with a non-leayyear calendar
       ! use the data 1 day before after FEB 28th (Julian day 59).
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. day>59 ) then
            day = day - 1
         end if

       ! set record index
         sec = 86400*(day-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

      if ( trim(groupby) == 'month' ) then

         if ( isleapyear(year) ) then
            months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
         else
            months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
         end if

       ! calculate initial month and day values
         call julian2monthday(year, day, month, mday)

       ! record in the next day, adjust year, month and second values
         if ( sec==86400 .AND. offset(var_i).EQ.0 ) then
            sec  = 0
            mday = mday + 1
            if ( mday > (months(month)-months(month-1)) ) then
               mday = 1
             ! bug found by Zhu Siguang, 05/25/2014
             ! move the below line in the 'else' statement
               !month = month + 1
               if (month == 12) then
                  month = 1
                  year = year + 1
               else
                  month = month + 1
               end if
            end if
         end if

       ! in case of leapyear with a non-leayyear calendar
       ! for day 29th Feb, use the data 1 day before, i.e., 28th FEB.
         if ( .NOT. leapyear .AND. isleapyear(year) .AND. month==2 .AND. mday==29 ) then
            mday = 28
         end if

       ! set record index
         sec    = 86400*(mday-1) + sec
         time_i = floor( (sec-offset(var_i)) *1. / dtime(var_i) ) + 1
      end if

      if (time_i < 0) then
         write(6, *) "got the wrong time record of forcing! stop!"; stop
      end if

      return

   END SUBROUTINE setstampUB

 ! calculate time average coszen value bwteeen [LB, UB]
 ! ------------------------------------------------------------
   SUBROUTINE calavgcos(lat_i, lon_i, lat_n, lon_n)

      implicit none
      integer, intent(in) :: lat_i, lon_i, lat_n, lon_n

      integer  :: i, j
      real(r8) :: calday, cosz
      type(timestamp) :: tstamp

      tstamp = tstamp_LB(7)
      avgcos(:,:) = 0._r8

      do while (tstamp < tstamp_UB(7))

         tstamp = tstamp + deltim
         calday = calendarday(tstamp)

! added by yuan, 05/23/2022
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,cosz)
#endif
         do i = lat_i, lat_i+lat_n-1
            do j = lon_i, lon_i+lon_n-1
               cosz = orb_coszen(calday, rlons(j), rlats(i))
               cosz = max(0.001, cosz)
               avgcos(j, i) = avgcos(j, i) + cosz*real(deltim)
            end do
         end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      end do
      avgcos = avgcos/real(dtime(7))

   END SUBROUTINE calavgcos

 ! read point forcing data
 ! ------------------------------------------------------------
   SUBROUTINE metreadpoint(idate,forcn,s_year,s_month,s_day,s_seconds,deltim,s_julian)

      implicit none

      integer , intent(in)    :: s_year, s_month, s_day, s_seconds, s_julian
      integer , intent(in)    :: idate(3)
      real(r8), intent(in)    :: deltim
      real(r8), intent(inout) :: forcn(:,:,:)
      
      ! var ids
      INTEGER  :: ncid, varid, swid, lwid, prid, snid, taid, qaid, wvid, &
                  wuid, psid
      INTEGER  :: days, time_i, time_s, i
      INTEGER  :: months(0:12)
      REAL(r8) :: rain(1), snow(1)
      REAL(r8) :: rtime
      REAL(r8) :: metadata(1)
      LOGICAL  :: firstread

      save time_i, firstread

      character(256) :: filename

      filename = trim(fmetdat)//trim(fmetnam)

      rtime = deltim !dtime(1)
#ifdef USE_POINT_DATA      
      IF(trim(fprefix)=='NC') THEN
         IF(time_i <= 0) THEN
            IF (s_year==startyr .and. s_month==startmo .and. s_day==startday .and. s_seconds==startsec) THEN
               IF (idate(1) == s_year) THEN
                  IF (idate(2) == s_julian) THEN
                     time_i = (idate(3)-s_seconds)/rtime+1
                  ELSE
                     time_i = ((86400-s_seconds)/rtime+1)+(idate(2)-1-s_julian)*(86400/rtime)+((idate(3)/rtime))
                  ENDIF
               ELSE
                  IF (isleapyear(s_year)) THEN
                     months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
                  ELSE
                     months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
                  ENDIF

                  time_i = (idate(2)-1)*(86400/rtime)+((idate(3)/rtime+1))+ &
                           ((86400-s_seconds)/rtime)+(months(12)-months(s_month-1)-s_day)*(86400/rtime)
                  
                  DO i = s_year+1, idate(1)
                     IF (isleapyear(i) .and. (idate(1)-i)>0) THEN
                        time_i = time_i+(86400/rtime)*366
                     ENDIF

                     IF (.NOT. isleapyear(i) .and. (idate(1)-i)>0) THEN
                        time_i = time_i+(86400/rtime)*365
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               ! restart case
               print*, idate(:)
               IF (s_year==startyr) THEN
                  time_i = (idate(2)-1)*(86400/rtime)+((idate(3)/rtime+1)) - startsec/rtime
               ELSE
                  IF (isleapyear(startyr)) THEN
                     months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
                     time_i = (86400/rtime)*366 - ((months(startmo-1)+startday-1)*(86400/rtime)+startsec/rtime)
                  ELSE
                     months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
                     time_i = (86400/rtime)*365 - ((months(startmo-1)+startday-1)*(86400/rtime)+startsec/rtime)
                  ENDIF
                ! print*, time_i, ((months(startmo)+startday-1)*(86400/rtime)+startsec/rtime)
                  time_i = time_i + (idate(2)-1)*(86400/rtime)+((idate(3)/rtime+1))

                  DO i = startyr+1, idate(1)

                     IF (isleapyear(i) .and. (idate(1)-i)>0) THEN
                        time_i = time_i+(86400/rtime)*366
                     ENDIF

                     IF (.NOT. isleapyear(i) .and. (idate(1)-i)>0) THEN
                        time_i = time_i+(86400/rtime)*365
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
            firstread = .True.
         ELSE
            time_i = time_i + 1
            firstread = .False.
         ENDIF

         ! print*, time_i

         IF (firstread) THEN
            CALL sanity( nf90_open(filename, nf90_nowrite, fid(1)) )
            firstread = .False.
         ENDIF

         DO i = 1, NVAR

            IF (vname(i)=='Rainf') THEN
               CALL sanity( nf90_inq_varid(fid(1), 'Rainf', prid) )
               CALL sanity( nf90_inq_varid(fid(1), 'Snowf', snid) )
               CALL sanity( nf90_get_var  (fid(1), prid, rain(:), start=(/1,1,time_i/), count=(/1,1,1/)) )
               CALL sanity( nf90_get_var  (fid(1), snid, snow(:), start=(/1,1,time_i/), count=(/1,1,1/)) )
               
               forcn(1,1,4) = rain(1) + snow(1)
            ELSE
               CALL sanity( nf90_inq_varid(fid(1), vname(i), varid) )
               CALL sanity( nf90_get_var  (fid(1), varid   , metadata(:), start=(/1,1,time_i/), count=(/1,1,1/)) )

               forcn(1,1,i) = metadata(1)
            ENDIF
         ENDDO
         ! close in metfina     
         ! CALL sanity( nf90_close(ncid) )
      ELSE
         IF(fid(1) ==  -1) THEN
            open(unit=11, file=filename, form='formatted', status='old', action='read')
            fid(1) = 11
         ENDIF

         read (fid(1), 10) forcn(1,1,7), forcn(1,1,8), forcn(1,1,4), forcn(1,1,1), &
                        forcn(1,1,5), forcn(1,1,6), forcn(1,1,3), forcn(1,1,2)

10       format (2f7.1, e14.3, 3f10.3, f10.1, e12.3)
      ENDIF
#endif
   END SUBROUTINE metreadpoint


 ! nc file open check
 ! ------------------------------------------------------------
   SUBROUTINE sanity(ret)

      implicit none
      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6, *) trim(nf90_strerror(ret)); stop
      end if

   END SUBROUTINE sanity

END MODULE METDAT
