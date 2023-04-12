#include <define.h>

MODULE UrbanAnthropogenic

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  USE timemanager, only: julian2monthday, isleapyear
  IMPLICIT NONE
  SAVE
  PRIVATE :: timeweek, gmt2local
  PUBLIC :: LUCY 

CONTAINS

  ! Using LUCY model to calculate vehicle heat and metabolic heat
  ! Allen et al., 2011

  Subroutine LUCY(idate,deltim,patchlonr,fix_holiday,week_holiday,hum_prof, &
               wdh_prof,weh_prof,popcell,vehicle,Fahe,vehc,meta)


   IMPLICIT NONE

   ! input vars
   INTEGER , INTENT(in) :: &
      idate(3)             ! calendar (year, julian day, seconds)
   REAL(r8), INTENT(in) :: &
      fix_holiday(365) , & ! Fixed public holidays, holiday(0) or workday(1)
      week_holiday(7)

   REAL(r8), INTENT(in) :: &
      deltim      , &
      patchlonr   , &
      ! f_fac       , &    ! factor multiplying total vehicles to get number of vehicles on road
      hum_prof(24), & ! Diurnal metabolic heat profile
      wdh_prof(24), &
      weh_prof(24), &
      popcell     , & ! grid population
      vehicle(3)

   REAL(r8) :: &
      vehc_prof(24,2), &
      carscell , & ! number of cars
      frescell , & ! number of freights
      mbkscell     ! nukber of motobikes
   
   REAL(r8), INTENT(out) :: &
      Fahe, &
      vehc, &
      meta

   REAL(r8) :: & 
      ahf_flx , &  ! metabolic heat + vehicle heat
      vehc_tot, &  
      vehc_flx     ! vehicle heat

   REAL(r8) ::  &
      londeg   , &
      car_sp   , &
      traf_frac, & ! vehicle heat profile of time step
      meta_prof, & ! metabolic heat profile of time step
      meta_flx , &
      carflx   , & ! heat from car
      motflx   , & ! heat from motorbike
      freflx       ! heat from freight
      !vehc_tot    ! total vehicle

   ! local vars
   INTEGER :: &
         ldate(3), &
         iweek   , &
         ihour   , &
         day     , &
         month   , &
         day_inx , &
         H2City  , &
         k       , &
         n       , &
         EC      , &
         EF      , &
         EM

   car_sp = 50
   EC = 3975  ! Sailor and Lu, 2004, J/m
   EM = 3975
   EF = 3975

#ifndef USE_POINT_DATA
   londeg = patchlonr*180/PI
   CALL gmt2local(idate, londeg, ldate)
#endif
   ! emmission factor of cars/mobike/freight in different speed
   !IF (car_sp == 16) THEN
   !   EC = 41.21
   !   EF = 162.13
   !   EM = 22.01
   !ELSE IF (car_sp == 24) THEN
   !   EC = 35.4
   !   EF = 141.21
   !   EM = 18.805
   !ELSE IF (car_sp == 32) THEN
   !   EC = 29.59
   !   EF = 120.29
   !   EM = 15.6
   !ELSE IF (car_sp == 40) THEN
   !   EC = 27.755
   !   EF = 114.355
   !   EM = 14.38
   !ELSE IF (car_sp == 48) THEN
   !   EC = 25.92
   !   EF = 108.42
   !   EM = 13.16
   !ELSE IF (car_sp == 56) THEN
   !   EC = 25.33
   !   EF = 106.685
   !   EM = 14.815
   !ELSE IF (car_sp == 64) THEN
   !   EC = 24.74
   !   EF = 104.95
   !   EM = 16.47
   !ELSE
      ! user difine
      ! EC =
      ! EF =
      ! EM =
   !ENDIF

   vehc_prof(:,1) = wdh_prof
   vehc_prof(:,2) = weh_prof
   !
   CALL julian2monthday(ldate(1), ldate(2), month, day)
   CALL timeweek(ldate(1), month, day, iweek)

   ! which hour of the day
   ! print*, 'lon is',londeg
   ! print*, 'gmt time', idate(:)
   ! print*, 'local time',ldate(:)
   ! print*, month, day, iweek
   
   ihour = CEILING(ldate(3)*1./3600)
   ! print*, ihour, ldate(3), idate(3), londeg

   ! public holiday or weekendday
   ! print*, day, iweek
   IF (day==366)  day=365
   IF (fix_holiday(day)==0 .or. week_holiday(iweek)==0) THEN
   !   H2city  = 1
      day_inx = 1
   ELSE
   !   H2city  = 0.8
      day_inx = 2
   ENDIF

   traf_frac = vehc_prof(ihour,day_inx)
   meta_prof = hum_prof (ihour)

   carscell = vehicle(1) !*popcell
   mbkscell = vehicle(2) !*popcell
   frescell = vehicle(3) !*popcell

   ! metabolic heat
   meta_flx= popcell*meta_prof/1e6
   ! Cars
   IF (carscell > 0) THEN
      carflx = carscell*popcell/1000!(carscell*24*f_fac)*(popcell/1000)
      carflx = carflx*traf_frac &
               *EC*(car_sp*1000)/1e6
      carflx = carflx/3600
   ENDIF
   ! Motorbikes
   IF (mbkscell > 0) THEN
      motflx = mbkscell*popcell/1000
      motflx = motflx*traf_frac &
               *EM*(car_sp*1000)/1e6
      motflx   = motflx/3600
   ENDIF
   ! Freight
   IF (frescell > 0)THEN
      freflx = frescell*popcell/1000
      freflx = freflx*traf_frac &
               *EF*(car_sp*1000)/1e6
      freflx = freflx/3600
   ENDIF
   !
   vehc = carflx + motflx + freflx
   meta = meta_flx
   Fahe = meta_flx + vehc !vehc_flx
   
   ! print*, Fahe
  
  END Subroutine LUCY

  SUBROUTINE gmt2local(idate, long, ldate)
    
    IMPLICIT NONE

    INTEGER , intent(in ) :: idate(3)
    REAL(r8), intent(in ) :: long
    INTEGER , intent(out) :: ldate(3)

    INTEGER  :: maxday
    REAL(r8) :: tdiff

    tdiff = long/15.*3600

    ldate(3) = idate(3) + tdiff

    IF (ldate(3) < 0) THEN

      ldate(3) = 86400 + ldate(3)
      ldate(2) = idate(2) - 1

      IF (ldate(2) < 1) THEN
         ldate(1) = idate(1) - 1
         IF ( isleapyear(ldate(1)) ) THEN
            ldate(2) = 366
         ELSE
            ldate(2) = 365
         ENDIF
      ENDIF

    ELSE IF (ldate(3) > 86400) THEN

      ldate(3) = ldate(3) - 86400
      ldate(2) = idate(2) + 1

      IF ( isleapyear(ldate(1)) ) THEN
         maxday = 366
      ELSE
         maxday = 365
      ENDIF

      IF(ldate(2) > maxday) THEN
         ldate(1) = idate(1) + 1
         ldate(2) = 1
      ENDIF
    ELSE
      ldate(2) = idate(2)
      ldate(1) = idate(1)
    ENDIF
   END SUBROUTINE gmt2local


  SUBROUTINE timeweek(year, month, day, iweek)

    IMPLICIT NONE

    INTEGER, intent(in ) :: year, month
    INTEGER, intent(out) :: iweek, day

    INTEGER :: myear, mmonth
    INTEGER :: judy1, judy2, judy3
    INTEGER :: yy, mm, dd, y12, y34
    INTEGER :: A, B, C, D, i

    INTEGER, save :: DayOfMonth(12) = [31,28,31,30,31,30,31,31,30,31,30,31]

    judy1 = mod(year, 400)
    judy2 = mod(year, 100)
    judy3 = mod(year, 4  )

    IF (judy2 == 0) THEN
      IF (judy1 == 0) THEN
         DayOfMonth(2) = 29
      ELSE
         DayOfMonth(2) = 28
      ENDIF
    ELSE
      IF (judy3 == 0) THEN
         DayOfMonth(2) = 29
      ELSE
         DayOfMonth(2) = 28
      ENDIF
    ENDIF

!#ifdef USE_POINT_DATA
!    IF (leap) THEN
!       DayOfMonth(2) = 29
!    ELSE
!       DayOfMonth(2) = 28
!    ENDIF
!#endif

    IF (month==1 .or. month==2) THEN
      mmonth = month + 12
      myear  = year  - 1
    ENDIF

    y12 = myear/100
    y34 = myear - y12*100

    A = int(y34/4.)
    B = int(y12/4.)
    C = y12*2
    D = int(26*(mmonth+1)/10.)

    iweek = abs(mod((y34+A+B-C+D+day-1), 7))

    DO i=1, month-1
       day = day + DayOfMonth(i)
    ENDDO

    IF (iweek == 0) THEN
       iweek = 7
    ENDIF

  END SUBROUTINE timeweek

END MODULE UrbanAnthropogenic
