#include <define.h>

MODULE UrbanAnthropogenic

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  USE timemanager, only: julian2monthday
  IMPLICIT NONE
  SAVE
  PRIVATE :: timeweek
  PUBLIC :: LUCY 

CONTAINS

  ! Using LUCY model to calculate vehicle heat and metabolic heat
  ! Allen et al., 2011

  Subroutine LUCY(idate,deltim,fix_holiday,week_holiday,f_fac,car_sp,hum_prof, &
               wdh_prof,weh_prof,popcell,vehicle,Fahe)  !vehc_tot,ahf_flx,vehc_flx)


   IMPLICIT NONE

   ! input vars
   INTEGER , INTENT(in) :: &
      idate(3)             ! calendar (year, julian day, seconds)
   REAL(r8), INTENT(in) :: &
      fix_holiday(365) , & ! Fixed public holidays, holiday(0) or workday(1)
      week_holiday(7)

   REAL(r8), INTENT(in) :: &
      deltim   , &
      f_fac    , &    ! factor multiplying total vehicles to get number of vehicles on road
      car_sp   , &    ! average speed of vehicles
      hum_prof(24), & ! Diurnal metabolic heat profile
      wdh_prof(24), &
      weh_prof(24), &
      popcell  , & ! grid population
      vehicle(3)
   REAL(r8) :: &
      vehc_prof(24,2), &
      carscell , & ! number of cars
      frescell , & ! number of freights
      mbkscell     ! nukber of motobikes
   
   REAL(r8), INTENT(out) :: Fahe

   REAL(r8) :: & 
      ahf_flx , &  ! metabolic heat + vehicle heat
      vehc_tot, &  
      vehc_flx     ! vehicle heat

   REAL(r8) ::  &
      traf_frac, & ! vehicle heat profile of time step
      meta_prof, & ! metabolic heat profile of time step
      meta_flx , &
      carflx   , & ! heat from car
      motflx   , & ! heat from motorbike
      freflx       ! heat from freight
      !vehc_tot    ! total vehicle

   ! local vars
   INTEGER :: &
         hdate(3), &
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

   ! emmission factor of cars/mobike/freight in different speed
   IF (car_sp == 16) THEN
      EC = 41.21
      EF = 162.13
      EM = 22.01
   ELSE IF (car_sp == 24) THEN
      EC = 35.4
      EF = 141.21
      EM = 18.805
   ELSE IF (car_sp == 32) THEN
      EC = 29.59
      EF = 120.29
      EM = 15.6
   ELSE IF (car_sp == 40) THEN
      EC = 27.755
      EF = 114.355
      EM = 14.38
   ELSE IF (car_sp == 48) THEN
      EC = 25.92
      EF = 108.42
      EM = 13.16
   ELSE IF (car_sp == 56) THEN
      EC = 25.33
      EF = 106.685
      EM = 14.815
   ELSE IF (car_sp == 64) THEN
      EC = 24.74
      EF = 104.95
      EM = 16.47
   ELSE
      ! user difine
      ! EC =
      ! EF =
      ! EM =
   ENDIF

   vehc_prof(:,1) = wdh_prof
   vehc_prof(:,2) = weh_prof
   !
   CALL julian2monthday(idate(1), idate(2), month, day)
   CALL timeweek(idate(1), month, day, iweek)

   ! which hour of the day
   IF (mod(deltim,3600.) == 0) THEN
      ihour = int(deltim/3600)
   ELSE
      ihour = int(deltim/3600) + 1
   ENDIF

   IF (ihour == 25) ihour = 1

   ! public holiday or weekendday
   IF (fix_holiday(day)==0 .or. week_holiday(iweek)==0) THEN
      H2city  = 1
      day_inx = 1
   ELSE
      H2city  = 0.8
      day_inx = 2
   ENDIF

   k = deltim/3600
   n = mod(deltim,3600.)
   
   IF (n == 0) THEN
      n         = 1
      traf_frac = vehc_prof(ihour,day_inx)
      meta_prof = hum_prof (ihour)
   ELSE
      n         = int(n/deltim)
      traf_frac = vehc_prof(ihour,day_inx) + &
                  (vehc_prof(ihour+1,day_inx)-vehc_prof(ihour,day_inx))*k*n
      meta_prof = hum_prof (ihour) + &
                  (hum_prof(ihour+1)-hum_prof(ihour))*k*n
   ENDIF

   carscell = vehicle(1) !*popcell
   mbkscell = vehicle(2) !*popcell
   frescell = vehicle(3) !*popcell

   ! metabolic heat
   meta_flx= popcell*meta_prof/1e6
   ! Cars
   carflx = (carscell*24*f_fac)*(popcell/1000)
   carflx = carflx*traf_frac &
            *EC*(car_sp*k*1000)/1e6
   carflx = carflx/3600
   ! Motorbikes
   motflx = (mbkscell*24*f_fac)*(popcell/1000)
   motflx = motflx*traf_frac &
            *EM*(car_sp*k*1000)/1e6
   motflx   = motflx/3600
   ! Freight
   freflx = (frescell*24*f_fac)*(popcell/1000)
   freflx = freflx*traf_frac &
            *EF*(car_sp*k*1000)/1e6
   freflx = freflx/3600

   !
   vehc_tot = carflx + motflx + freflx
   vehc_flx = vehc_tot*H2city
   Fahe = meta_flx + vehc_flx
  
  END Subroutine LUCY

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

    ! IF (leap) THEN
    !    DayOfMonth(2) = 29
    ! ELSE
    !    DayOfMonth(2) = 28
    ! ENDIF

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

    DO i=1, month
      day = day + DayOfMonth(i)
    ENDDO

    IF (iweek == 0) THEN
    iweek = 7
    ENDIF

  END SUBROUTINE timeweek

END MODULE UrbanAnthropogenic
