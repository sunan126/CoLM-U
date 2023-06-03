#include <define.h>

MODULE UrbanAnthropogenic

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  USE timemanager, only: julian2monthday, isleapyear
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
      deltim      , &
      f_fac       , &    ! factor multiplying total vehicles to get number of vehicles on road
      car_sp      , &    ! average speed of vehicles
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
         n       , &
         EC      , &
         EF      , &
         EM
   REAL(r8) :: k

   ! emmission factor
   EC = 3975        ! Sailor & Lu, 2004; unit, J/m
   EM = 3975
   EF = 3975

   vehc_prof(:,1) = wdh_prof
   vehc_prof(:,2) = weh_prof
   !
   CALL julian2monthday(idate(1), idate(2), month, day)
   CALL timeweek(idate(1), month, day, iweek)

   ! which hour of the day
   !IF (mod(deltim,3600.) == 0) THEN
      ihour = floor((idate(3)/3600)*1.) + 1
   !ELSE
   !   ihour = int(deltim/3600) + 1
   !ENDIF
   !day = idate(2)
   IF (day == 366) day = 365
   IF (ihour == 25) ihour = 1

   ! public holiday or weekendday
   IF (fix_holiday(day)==0 .or. week_holiday(iweek)==0) THEN
      !H2city  = 1
      day_inx = 1
   ELSE
      !H2city  = 0.8
      day_inx = 2
   ENDIF

   k = deltim/3600
   n = mod(deltim,3600.)
   print*, k 
   !IF (n == 0) THEN
   !   n         = 1
      traf_frac = vehc_prof(ihour,day_inx)
      meta_prof = hum_prof (ihour)
   !ELSE
   !   n         = int(n/deltim)
   !   traf_frac = vehc_prof(ihour,day_inx) + &
   !               (vehc_prof(ihour+1,day_inx)-vehc_prof(ihour,day_inx))*k*n
   !   meta_prof = hum_prof (ihour) + &
   !               (hum_prof(ihour+1)-hum_prof(ihour))*k*n
   !ENDIF

   carscell = 669!vehicle(1) !*popcell
   mbkscell = vehicle(2) !*popcell
   frescell = vehicle(3) !*popcell

   ! metabolic heat
   meta_flx= popcell*meta_prof/1e6*k
   ! Cars
   carflx = carscell*(popcell/1000)
   carflx = carflx*traf_frac &
            *EC*(car_sp*1000)/1e6*k
   carflx = carflx/3600
   print*, carflx
   ! Motorbikes
   motflx = mbkscell*(popcell/1000)
   motflx = motflx*traf_frac &
            *EM*(car_sp*1000)/1e6*k
   motflx   = motflx/3600
   ! Freight
   freflx = frescell*(popcell/1000)
   freflx = freflx*traf_frac &
            *EF*(car_sp*1000)/1e6*k
   freflx = freflx/3600

   !
   vehc_tot = carflx !+ motflx + freflx
   print*, vehc_tot
   Fahe = meta_flx + vehc_tot
  
  END Subroutine LUCY

  SUBROUTINE timeweek(year, month, day, iweek)

    IMPLICIT NONE

    INTEGER, intent(in ) :: year, month
    INTEGER, intent(out) :: iweek, day

    INTEGER :: myear, mmonth
    INTEGER :: judy1, judy2, judy3
    INTEGER :: yy, mm, dd, y12, y34
    INTEGER :: A, B, C, D, i

    INTEGER, save :: DayOfMonth(0:12)

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

    IF ( isleapyear(year) ) THEN
         DayOfMonth = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      ELSE
         DayOfMonth = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      ENDIF

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

    day = day + DayOfMonth(month-1)

    IF (iweek == 0) THEN
    iweek = 7
    ENDIF

  END SUBROUTINE timeweek

END MODULE UrbanAnthropogenic
