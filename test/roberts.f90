!>
!  Test code for the Jacchia-Roberts model.
!  Uses the older INPE code for the Jacchia-Roberts model
!  as a reference for testing the new implementation.
!
!  Source:
!  * Original [dead link]: http://www.dem.inpe.br/~val/atmod/default.html
!    The Orbital Dynamics group of INPE (Brazilian National Institute for Space Research)
!  * Archive: https://github.com/jacobwilliams/INPE-atmosphere-models

module inpe_roberts_module

   use space_weather_module,       only: sw_data_type, flux_data_type
   use jacchia_roberts_module,     only: prepare_flux_data
   use jacchia_roberts_kinds,      only: dp

   implicit none

   !> Standard MJD = MJD-1950 + 33282  (JD-2400000.5 = (JD-2433282.5) + 33282)
   real(dp), parameter :: MJD_1950_OFFSET = 33282.0_dp

   real(dp),parameter :: avog = 6.02217e+26_dp
   real(dp),dimension(6),parameter :: wm = [4.0026_dp  , &
                                            31.9988_dp , &
                                            28.0134_dp , &
                                            39.948_dp  , &
                                            15.9994_dp , &
                                            1.00797_dp]

   ! Global state (ok, since we are only using this for testing/validation)
   type(sw_data_type), save :: sw_global
   logical,            save :: sw_initialized = .false.

   private

   public :: soflud_init
   public :: soflud
   public :: rdymos_cssi
   public :: rsdamo

   contains

!----
!
SUBROUTINE rdymos(Sa,Su,Rjud,Dafr,Gsti,Te,Ad,Wmol,Rhod)
!
!------
!
!
! PURPOSE:
!
!     THE SUBROUTINE RDYMOS GIVES THE  TEMPERATURE,  DENSITY
!     AND MOLECULAR WEIGHT OF ATMOSPHERE, USING THE  ROBERTS
!     VERSION (2) OF THE JACCHIA 70 DYNAMIC MODEL, WITH  THE
!     SOLAR  FLUX DATA FILE.
!
! INPUTS:
!
!     SA(1) RIGHT ASCENSION OF THE  POINT  IN  QUESTION,  IN
!           RADIANS (0 TO 2.*PI)
!     SA(2) DECLINATION (GEOCENTRIC LATITUDE) OF THE  POINT,
!           IN RADIANS (-PI TO PI).
!     SA(3) GEOCENTRIC ALTITUDE OF THE POINT IN  METERS,  IN
!           THE RANGE 110000. TO 2000000.M.
!     SU(1) RIGHT ASCENSION OF THE SUN AT THE DATE, IN RADI-
!           ANS (O TO 2.*PI).
!     SU(2) SUN DECLINATION IN RADIANS (-PI TO PI)
!     RJUD  MODIFIED JULIAN DATE (IF OUT OF RANGE, THE  SUB-
!           ROUTINE WILL PRINT A MESSAGE AND STOP).
!     DAFR  TIME (UT) OF THE DAY, IN SECONDS.
!     GSTI  GREENWICH SIDERAL TIME, IN RADIANS (0 TO 2.*PI),
!           AT THE TIME DAFR OF THE DATE RJUD.(NOT USED. FOR
!           COMPATIBILITY PURPOSE WITH  OTHER MODELS ONLY)
!
! OUTPUTS:
!
!     TE(1) EXOSPHERIC  TEMPERATURE  ABOVE  THE   POINT   IN
!           QUESTION, AS DEFINED IN REFERENCE (1), IN KELVIN
!     TE(2) LOCAL TEMPERATURE AROUND THE POINT, IN KELVIN.
!     AD(1) LOGARITHM IN BASE 10 OF THE HE NUMBER-DENSITY.
!     AD(2) LOGARITHM IN BASE 10 OF THE O2 NUMBER-DENSITY.
!     AD(3) LOGARITHM IN BASE 10 OF THE N2 NUMBER-DENSITY.
!     AD(4) LOGARITHM IN BASE 10 OF THE AR NUMBER-DENSITY.
!     AD(5) LOGARITHM IN BASE 10 OF THE O  NUMBER-DENSITY.
!     AD(6) LOGARITHM IN BASE 10 OF THE H  NUMBER-DENSITY.
!     WMOL  MEAN-MOLECULAR-WEIGHT OF THE ATMOSPHERE  AT  THE
!           POINT IN KG/KGMOL.
!     RHOD  MEAN-MASS-DENSITY OF THE ATMOSPHERE, IN KG/M/M/M
!
! OBS:
!
!     HE    HELIUM
!     O2    MOLECULAR OXYGEN
!     N2    MOLECULAR NITROGEN
!     AR    ARGON
!     O     ATOMIC OXYGEN
!     H     ATOMIC HYDROGEN
!
! REFERENCES:
!
!     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
!         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE
!         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
!         VERLAG, 1972. PART 3, P. 227-338.
!
!     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
!         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
!         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
!         DEC. 1971.
!
! AUTHORS:
!
!     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!
!     APR. 1989              V. 1.0
!
   IMPLICIT NONE

   real(dp) Ad , dafl , Dafr , Gsti , outr , Rhod , rjfl , Rjud , Sa , sd , sf , Su , tauo , Te , Wmol
   INTEGER int , nd
   DIMENSION Sa(3) , Su(2) , Te(2) , Ad(6)
   DIMENSION sf(3) , sd(15)
!
!------
!

   rjfl = 0.0_dp  ! this was uninitiated in the original code

   rjfl = rjfl - 1.0_dp
   IF ( Dafr<61344.0_dp ) THEN
      rjfl = rjfl - 1.0_dp
      dafl = Dafr + 25056.0_dp
   ELSE
      dafl = Dafr - 61344.0_dp
   ENDIF

   CALL soflud(rjfl,dafl,sd,outr)

   IF ( outr/=0. ) THEN
      WRITE (6,*) ' ERROR IN ROUTINE JDYMOS: OUTR = ' , int(outr)
      STOP
   ENDIF

   tauo = 6.696_dp
   nd = (Dafr/3600.0_dp-sd(6)+12.0_dp-tauo)/3.0_dp
   sf(1) = sd(9)
   sf(2) = sd(11)
   sf(3) = sd(nd)

   CALL rsdamo(Sa,Su,sf,Rjud,Dafr,Gsti,Te,Ad,Wmol,Rhod)

END SUBROUTINE rdymos

SUBROUTINE rsmods(Altu,Rjud,Dafr,Te,Al,Wmol,Rhod)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE RSMODS USES THE ROBERTS VERSION (2)  OF
!     THE JACCHIA 70 STATIC AND DYNAMIC MODEL (1), WITH  THE
!     SOLAR FLUX DATA TO OBTAIN THE  TEMPERATURE,  MOLECULAR
!     WEIGHT AND DENSITY OF LOCAL ATMOSPHERE.
!
! INPUTS:
!
!     ALTU  ALTITUDE  OF  THE  POINT  IN  METERS  (110000 TO
!           2000000).
!     RJUD  MODIFIED  JULIAN  DATE  (IF  OUT  OF  RANGE  THE
!           ROUTINE WILL PRINT A MESSAGE AND STOP).
!     DAFR  TIME (UT) OF THE DAY, IN SECONDS (0  TO  86400).
!
! OUTPUTS:
!
!     TE(1) EXOSPHERIC  TEMPERATURE  ABOVE  THE   POINT   IN
!           QUESTION, AS DEFINED IN REFERENCE (1), IN KELVIN
!     TE(2) LOCAL TEMPERATURE AROUND THE POINT, IN KELVIN.
!     AD(1) LOGARITHM IN BASE 10 OF THE HE NUMBER-DENSITY.
!     AD(2) LOGARITHM IN BASE 10 OF THE O2 NUMBER-DENSITY.
!     AD(3) LOGARITHM IN BASE 10 OF THE N2 NUMBER-DENSITY.
!     AD(4) LOGARITHM IN BASE 10 OF THE AR NUMBER-DENSITY.
!     AD(5) LOGARITHM IN BASE 10 OF THE O  NUMBER-DENSITY.
!     AD(6) LOGARITHM IN BASE 10 OF THE H  NUMBER-DENSITY.
!     WMOL  MEAN-MOLECULAR-WEIGHT OF THE ATMOSPHERE  AT  THE
!           POINT IN KG/KGMOL.
!     RHOD  MEAN-MASS-DENSITY OF THE ATMOSPHERE, IN KG/M/M/M
!
! SUBCALLS:
!
!     SOFLUD
!     RSMADE
!
! REFERENCES:
!
!     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
!         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE
!         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
!         VERLAG, 1972. PART 3, P. 227-338.
!
!     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
!         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
!         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
!         DEC. 1971.
!
! AUTHORS:
!
!     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!
!     APR. 1989              V. 1.0
!
   IMPLICIT NONE

   real(dp) Al , Altu , Dafr , outr , Rhod , &
            Rjud , sd , sf , Te , vari , Wmol
   INTEGER int
   DIMENSION Te(2) , Al(6)
   DIMENSION sf(3) , sd(15)
!
!------
!
   CALL soflud(Rjud,Dafr,sd,outr)

   IF ( outr/=0.0_dp ) THEN
      WRITE (6,*) ' ERROR IN ROUTINE ISMODS: OUTR = ' , int(outr)
      STOP
   ENDIF

   sf(1) = sd(9)
   sf(2) = sd(11)
   vari = 0.154_dp*sd(7)
   sf(3) = 1.89_dp*log(vari+sqrt(vari*vari+1.0_dp))

   CALL rsmade(Altu,sf,Te,Al,Wmol,Rhod)

END SUBROUTINE rsmods


SUBROUTINE rsdamo(Sa,Su,Sf,Rjud,Dafr,Gsti,Te,Ad,Wmol,Rhod)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE  RSDAMO  GIVES  THE  DENSITY, MOLECULAR
!     WEIGHT AND TEMPERATURE OF THE UPPER ATMOSPHERE,  USING
!     THE ROBERTS VERSION (2) OF THE JACCHIA 70  STATIC  AND
!     DYNAMIC ATMOSPHERIC MODEL.
!
! INPUTS:
!
!       SA(1)    RIGHT ASCENTION OF THE POINT, IN RADIANS.
!       SA(2)    DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
!                POINT, IN RADIANS (-PI TO PI).
!       SA(3)    GEOCENTRIC ALTITUDE IN METERS, BETWEEN  THE
!                RANGE 110,000-2,000,000.
!       SU(1)    RIGHT ASCENTION OF THE SUN AT THE DATE,  IN
!                RADIANS (0 TO 2*PI).
!       SU(2)    SUN DECLINATION IN RADIANS (-PI TO PI).
!       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
!                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
!                W/M/M/HZ.
!       SF(2)    AVERAGED DAILY OBSERVED FLUX  AS DEFINED BY
!                JACCHIA, IN 1E-22 W/M/M/HZ.
!       SF(3)    3-HOURLY PLANETARY GEOMAGNETIC INDEX KP, AT
!                THE TIME 0.279 DAYS EARLIER.
!       RJUD     MODIFIED JULIAN  DATE,  REFERED  TO  1950.0
!                (JULIAN DATE-2433282.5).
!       DAFR     TIME (UT) OF THE DAY, IN SECONDS.
!       GSTI     GREENWICH SIDEREAL TIME, IN RADIANS, AT THE
!                TIME DAFR OF THE DATE RJUD (0 TO 2*PI).(NOT
!                USED. FOR COMPATIBILITY PURPOSE WITH  OTHER
!                MODELS ONLY)
!
! OUTPUTS:
!
!       TE(1)    EXOSPHERIC TEMPERATURE ABOVE THE  POINT  AS
!                DEFINED BY JACCHIA'S 70 MODEL, IN KELVIN.
!       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
!                KELVIN.
!       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
!       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
!       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
!       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
!       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
!       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
!       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
!                THE POINT, IN KG/KGMOL.
!       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
!                POINT, IN KG/M/M/M.
! OBS:
!       HE       HELIUM
!       O2       MOLECULAR OXYGEN
!       N2       MOLECULAR NITROGEN
!       AR       ARGON
!       O        ATOMIC OXYGEN
!       H        ATOMIC HYDROGEN
!
! SUBCALLS:
!
!       DYJRMO
!
! REFERENCES:
!
!     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
!         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE
!         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
!         VERLAG, 1972. PART 3, P. 227-338.
!
!     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
!         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
!         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
!         DEC. 1971.
!
! AUTHORS:
!
!     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!
!     APR. 1989              V. 1.0
!
   IMPLICIT NONE

   real(dp) Ad , al , amjd , Dafr , Gsti , Rhod , Rjud , Sa , Sf , Su , Te , Wmol

   DIMENSION Sa(3) , Su(2) , Sf(3) , Te(2) , Ad(6)
   DIMENSION al(6)

!
!------
!

   amjd = Rjud + 33282.0_dp + Dafr/86400.0_dp

   CALL dyjrmo(amjd,Su,Sa,Sf,Te,al,Wmol,Rhod)

   Ad(1) = al(3)
   Ad(2) = al(4)
   Ad(3) = al(1)
   Ad(4) = al(2)
   Ad(5) = al(5)
   Ad(6) = al(6)

END SUBROUTINE rsdamo


SUBROUTINE rsmade(Altu,Sf,Te,Ad,Wmol,Rhod)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
!     HEIGHTS FROM 110 TO 2000 KM, USING THE ROBERTS VERSION
!     (2) OF THE JACCHIA 70  STATIC  MODEL  TO  COMPUTE  THE
!     ATMOSPHERIC DENSITY (1).
!
! INPUTS:
!
!       ALTU     ALTITUDE OF THE POINT IN METERS.
!       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
!                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
!                W/M/M/HZ.
!       SF(2)    AVERAGED DAILY OBSERVED FLUX  AS DEFINED BY
!                JACCHIA, IN 1E-22 W/M/M/HZ.
!
! OUTPUTS:
!
!       TE(1)    EXOSPHERIC TEMPERATURE ABOVE THE  POINT  AS
!                DEFINED BY JACCHIA'S 70 MODEL, IN KELVIN.
!       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
!                KELVIN.
!       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
!       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
!       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
!       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
!       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
!       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
!       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
!                THE POINT, IN KG/KGMOL.
!       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
!                POINT, IN KG/M/M/M.
!
! OBS:
!       HE       HELIUM
!       O2       MOLECULAR OXYGEN
!       N2       MOLECULAR NITROGEN
!       AR       ARGON
!       O        ATOMIC OXYGEN
!       H        ATOMIC HYDROGEN
!
! SUBCALLS:
!
!       RMOWEI
!
! REFERENCES:
!
!     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
!         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE
!         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
!         VERLAG, 1972. PART 3, P. 227-338.
!
!     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
!         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
!         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
!         DEC. 1971.
!
! AUTHORS:
!
!     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!
!     APR. 1987              V. 1.0
!     AUG. 2011              V. 1.1 (INCLUDED MISSING TEMPERATURE)
!
   IMPLICIT NONE

   real(dp) :: Ad(6) , al(6) , Altu , anac , anut , fbar , flux , &
            heig , Rhod , Sf(3) , Te(2) , thaf , tz , weig , Wmol
   INTEGER :: ic
!
!------
!

   flux = Sf(1)
   fbar = Sf(2)
   thaf = 379.0_dp + 3.24_dp*fbar + 1.3_dp*(flux-fbar)
   heig = Altu/1.0e3_dp

   CALL stjrmo(thaf,heig,tz,al)

   Ad(1) = al(3)
   Ad(2) = al(4)
   Ad(3) = al(1)
   Ad(4) = al(2)
   Ad(5) = al(5)
   Ad(6) = al(6)
   anut = 0.0_dp
   weig = 0.0_dp

   DO ic = 1 , 6
      anac = 10.0_dp**Ad(ic)
      weig = weig + wm(ic)*anac
      anut = anut + anac
   ENDDO

   Wmol = weig/anut
   Rhod = weig/avog
   Te(1) = thaf
   Te(2) = tz

END SUBROUTINE rsmade

SUBROUTINE rmowei(Tinf,Heig,Ad,Wmol,Rhod)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
!
!     HEIGHTS FROM 110 TO 2000 KM, USING THE ROBERTS VERSION
!                                            -
!     (2) OF THE JACCHIA'S 70 STATIC MODEL  TO  COMPUTE  THE
!                                    --
!     ATMOSPHERIC DENSITY AND MOLECULAR WEIGHT (1).
!                                       ---
!
! INPUTS:
!
!       TINF     EXOSPHERIC  TEMPERATURE   AS   DEFINED   BY
!                JACCHIA'S 1970 MODEL, IN KELVIN.
!       HEIG     ALTITUDE OF THE POINT IN KM.
!
! OUTPUTS:
!
!       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
!       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
!       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
!       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
!       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
!       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
!       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
!                THE POINT, IN KG/KGMOL.
!       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
!                POINT, IN KG/M/M/M.
!
! OBS:
!       HE       HELIUM
!       O2       MOLECULAR OXYGEN
!       N2       MOLECULAR NITROGEN
!       AR       ARGON
!       O        ATOMIC OXYGEN
!       H        ATOMIC HYDROGEN
!
! SUBCALLS:
!
!       STJRMO
!
! REFERENCES:
!
!     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
!         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE
!         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
!         VERLAG, 1972. PART 3, P. 227-338.
!
!     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
!         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
!         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
!         DEC. 1971.
!
! AUTHORS:
!
!     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!
!     APR. 1989              V. 1.0
!
   IMPLICIT NONE

   real(dp) :: Ad(6) , al(6) , anac , anut , Heig , Rhod , Tinf , tz , weig , Wmol
   INTEGER :: ic

!
!------
!

   CALL stjrmo(Tinf,Heig,tz,al)

   Ad(1) = al(3)
   Ad(2) = al(4)
   Ad(3) = al(1)
   Ad(4) = al(2)
   Ad(5) = al(5)
   Ad(6) = al(6)
   anut = 0.0_dp
   weig = 0.0_dp

   DO ic = 1 , 6
      anac = 10.0_dp**Ad(ic)
      weig = weig + wm(ic)*anac
      anut = anut + anac
   ENDDO

   Wmol = weig/anut
   Rhod = weig/avog
END SUBROUTINE rmowei


SUBROUTINE dyjrmo(Djm,Sun,Sat,Geo,Temp,Dn,Amw,Dens)
!
!-----
!
!  PURPOSE : COMPUTATION OF THE ATMOSPHERIC PROPERTIES
!            ACCORDING TO THE ANALYTICAL ROBERTS(1972)
!            METHOD APPLIED TO THE JACCHIA(1971) MODEL
!
! INPUTS : DJM ... MODIFIED JULIAN DATE DJM=JD-2400000.5
!          SUN(1). RIGHT ASCENSION OF SUN (RAD)
!          SUN(2). DECLINATION     OF SUN (RAD)
!          SAT(1). RIGHT ASCENSION OF THE POINT (RAD)
!          SAT(2). DECLINATION     OF THE POINT (RAD)
!          SAT(3). ALTITUDE        OF THE POINT (M)
!          GEO(1). 10.7 CM SOLAR FLUX,IN UNITS OF
!                  1.E-22 WATTS M**2/HERTZ , FOR A
!                  TABULAR TIME 1.71 DAYS EARLIER
!          GEO(2). 10.7 CM SOLAR FLUX AVERAGED OVER
!                  FOUR SOLAR ROTATIONS,CENTERED ON
!                  THE PRESENT TIME
!          GEO(3). GEOMAGNETIC PLANETARY THREE HOUR
!                  RANGE  INDEX "KP" FOR A TABULAR
!                  TIME 0.279 DAYS EARLIER
!
! OUTPUTS : TEMP(1)...EXOSPHERIC TEMPERATURE (KELVIN)
!           TEMP(2)...LOCAL TEMPERATURE      (KELVIN)
!           DN(1) ... LOG10 OF N2 DENSITY NUMBER (M**-3)
!           DN(2) ... LOG10 OF A  DENSITY NUMBER
!           DN(3) ... LOG10 OF HE DENSITY NUMBER
!           DN(4) ... LOG10 OF O2 DENSITY NUMBER
!           DN(5) ... LOG10 OF O  DENSITY NUMBER
!           DN(6) ... LOG10 OF H  DENSITY NUMBER
!           AMW   ... MEAN MOLECULAR WEIGHT (KG/KGMOL)
!           DENS  ... ATMOSPHERIC DENSITY (KG/M**3)
!
! REF. KUGA,H.K. "REFORMULACAO COMPUTACIONAL DO
!         MODELO DE JACCHIA-ROBERTS PARA A DEN-
!         SIDADE ATMOSFERICA".INPE,SAO JOSE DOS
!         CAMPOS, OCT.1985 (INPE-3691-RPE/493).
!
!      JACCHIA,L.G."REVISED STATIC MODELS OF THE
!         THERMOSPHERE AND EXOSPHERE WITH EMPIRICAL
!         TEMPERATURE PROFILES."SAO,CAMBRIDGE,MA,
!         1971.SAO SPECIAL REPORT NO. 332
!
!      ROBERTS JR,C.E."AN ANALYTIC MODEL FOR UPPER
!         ATMOSPHERE DENSITIES BASED UPON JACCHIA'S
!         1970 MODELS."CELESTIAL MECHANICS 4:368-
!         377,1971
!
! AUTHOR : HELIO KOITI KUGA - JUNE 1985 -INPE-DMC/DDO
!
! DIMENSION ARRAY'S, VARIABLES AND CONSTANTS
!
   IMPLICIT NONE

   real(dp) abs , Amw , capphi , cos , d1 , d2 , d3 , d4 , d5 , d6 , Dens , df , Djm , dlhe , dlr , dlr20 , dlrgm , dlrsa , &
            dlrsl
   real(dp) Dn , dtg , dtg18 , dtg20 , eta , expkp , f , fdfz , fs , fsm , gdft , Geo , h , pk , s , Sat
   real(dp) sat1 , sat2 , sat3 , sign , sin , sumn , sumnm , Sun , sun1 , sun2 , tanh , tau , Temp , theta , tinf , tsubc , tsubl ,  &
            tz
   INTEGER mod
   DIMENSION Sun(2) , Sat(3) , Geo(3) , Temp(2) , Dn(6)

   real(dp),parameter :: pi = acos(-1.0_dp)
   real(dp),parameter :: piv2 = 2.0_dp * pi
   real(dp),parameter :: piv4 = 4.0_dp * pi
   real(dp),parameter :: pid4 = pi / 4.0_dp
   real(dp),parameter :: cons25 = sin(pid4)**3
!
   sun1 = Sun(1)
   sun2 = Sun(2)
   sat1 = Sat(1)
   sat2 = Sat(2)
   sat3 = Sat(3)/1000.0_dp
   fs = Geo(1)
   fsm = Geo(2)
   pk = Geo(3)
!
!       MINIMUM NIGHT-TIME TEMPERATURE OF THE GLOBAL
!       EXOSPHERIC TEMPERATURE DISTRIBUTION WHEN THE
!       GEOMAGNETIC ACTIVITY INDEX KP = 0
!       EQUATION 14J
!
   tsubc = 379.0_dp + 3.24_dp*fsm + 1.3_dp*(fs-fsm)
!
!       EQUATION 15J
!
   eta = 0.5_dp*abs(sat2-sun2)
   theta = 0.5_dp*abs(sat2+sun2)
!
!       EQUATION 16J
!
   h = sat1 - sun1
   tau = h - 0.64577182_dp + 0.10471976_dp*sin(h+0.75049158_dp)
!
!       EXOSPHERIC TEMPERATURE TSUBL WITHOUT CORRECTION
!       FOR GEOMAGNETIC ACTIVITY
!       EQUATION 17J
!
   s = sin(theta)**2.2_dp
   df = s + (cos(eta)**2.2_dp-s)*abs(cos(0.5*tau))**3
   tsubl = tsubc*(1.0_dp+0.3_dp*df)
!
!       EQUATION 18J
!
   expkp = exp(pk)
   dtg18 = 28.0_dp*pk + 0.03_dp*expkp
!
!       EQUATION 20J
!
   dtg20 = 14.0_dp*pk + 0.02_dp*expkp
   dlr20 = 0.012_dp*pk + 1.2e-05_dp*expkp
!
!       THE FOLLOWING STATEMENTS EFFECT A CONTINUOUS
!       TRANSITION FROM EQ. 20J AT HEIGHTS WELL BELOW
!       350 KM TO EQ. 18J AT HEIGHTS WELL ABOVE
!       350 KM .
!
   f = 0.5_dp*(tanh(0.04_dp*(sat3-350.0_dp))+1.0_dp)
   dlrgm = dlr20*(1.0_dp-f)
   dtg = dtg20*(1.0_dp-f) + dtg18*f
!
!       EXOSPHERIC TEMPERATURE
!
   tinf = tsubl + dtg
!
!   STATIC MODEL OF JACCHIA-ROBERTS FOR THE
!   ATMOSPHERIC DENSITY
!
   CALL stjrmo(tinf,sat3,tz,Dn)
!
!   EQ. 23J   PHASE OF THE SEMI-ANNUAL VARIATION
!
   capphi = mod((Djm-36204.0_dp)/365.2422_dp,1.0_dp)
!
!   EQ. 22J
!
   tau = capphi + 0.09544_dp*((0.5_dp+0.5_dp*sin(piv2*capphi+6.035_dp))**1.650_dp-0.5_dp)
   gdft = 0.02835_dp + 0.3817_dp*(1.0_dp+0.4671_dp*sin(piv2*tau+4.137_dp))*sin(piv4*tau+4.259_dp)
   fdfz = (5.876e-07_dp*sat3**2.331_dp+0.06328_dp)*exp(-2.868e-03_dp*sat3)
!
!   EQ. 21J  SEMI-ANNUAL VARIATION
!
   dlrsa = fdfz*gdft
!
!   EQ. 24J  SEASONAL-LATITUDINAL VARIATION OF THE
!            LOWER THERMOSPHERE
!
   dlrsl = 0.014_dp*(sat3-90.0_dp)*exp(-0.0013_dp*(sat3-90.0_dp)**2)*sign(1.0_dp,sat2)*sin(piv2*capphi+1.72_dp)*sin(sat2)**2
!
!   SUM THE CORRECTIONS AND APPLY TO THE
!   NUMBER DENSITIES
!
   dlr = dlrgm + dlrsa + dlrsl
   Dn(1) = Dn(1) + dlr
   Dn(2) = Dn(2) + dlr
   Dn(3) = Dn(3) + dlr
   Dn(4) = Dn(4) + dlr
   Dn(5) = Dn(5) + dlr
   Dn(6) = Dn(6) + dlr
!
!   EQ. 25J  SEASONAL-LATITUDINAL VARIATION
!            OF HELIUM
!
   dlhe = 0.65_dp*abs(sun2/0.4091609_dp)*(sin(pid4-0.5_dp*sat2*sign(1.0_dp,sun2))**3-cons25)
   Dn(3) = Dn(3) + dlhe
!
!  COMPUTE DENSITY AND MEAN MOLECULAR WEIGHT
!
   d1 = 10.0_dp**Dn(1)
   d2 = 10.0_dp**Dn(2)
   d3 = 10.0_dp**Dn(3)
   d4 = 10.0_dp**Dn(4)
   d5 = 10.0_dp**Dn(5)
   d6 = 10.0_dp**Dn(6)
   sumn = d1 + d2 + d3 + d4 + d5 + d6
   sumnm = 28.0134_dp*d1 + 39.9480_dp*d2 + 4.0026_dp*d3 + 31.9988_dp*d4 + 15.9994_dp*d5 + 1.00797_dp*d6
   Amw = sumnm/sumn
   Dens = sumnm/6.02257e+26_dp
   Temp(1) = tinf
   Temp(2) = tz
!
END SUBROUTINE dyjrmo

SUBROUTINE stjrmo(Tinf,Sat3,Tz,Dn)
!
!-----
!
!  PURPOSE : STATIC MODEL FOR CALCULATION OF
!            ATMOSPHERIC PROPERTIES AT A GIVEN
!            ALTITUDE.
!
!  INPUTS : TINF...EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3...ALTITUDE               (KM)
!
!  OUTPUTS : TZ  ... LOCAL TEMPERATURE    (KELVIN)
!            DN  ... LOG10 OF DENSITY NUMBERS (M**-3)
!            DN(1) ... N2
!            DN(2) ... A
!            DN(3) ... HE
!            DN(4) ... O2
!            DN(5) ... O
!            DN(6) ... H
!
!  AUTHOR : HELIO KOITI KUGA - JUNE 1985
!           INPE - DMC/DDO
!
   IMPLICIT NONE

   real(dp) Dn , Sat3 , Tinf , Tz
   DIMENSION Dn(6)

   IF ( Sat3>125.0_dp ) THEN
!
      CALL stjr03(Tinf,Sat3,Tz,Dn)
      RETURN
   ELSEIF ( Sat3>100.0_dp ) THEN
!
      CALL stjr02(Tinf,Sat3,Tz,Dn)
      RETURN
   ELSEIF ( Sat3<90.0_dp ) THEN
!
      PRINT 99001
99001 FORMAT (1X,'ATENCA0 : MENSAGEM DA ROTINA DE ',/,1X,'*******   CALCULO DA DENSIDADE AT-',/,1X,'          MOSFERICA.',//,3X,   &
             &'ALTITUDE DO SATELITE MENOR QUE 90 KM')
      RETURN
   ENDIF
   CALL stjr01(Tinf,Sat3,Tz,Dn)

END SUBROUTINE stjrmo


SUBROUTINE stjr01(Tinf,Sat3,Tl2,Al10n)
!
!-----
!
!  PURPOSE : STATIC MODEL FOR CALCULATION OF
!            ATMOSPHERIC PROPERTIES FOR THE
!            BAND FROM 90 TO 100 KM.
!
!  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3 ... ALTITUDE               (KM)
!
!  OUTPUTS : TL2 ... LOCAL TEMPERATURE (KELVIN)
!            AL10N . ALOG10 OF DENSITY NUMBERS (M**-3)
!
!            AL10N(1) ... N2
!            AL10N(2) ... A
!            AL10N(3) ... HE
!            AL10N(4) ... O2
!            AL10N(5) ... O
!            AL10N(6) ... H
!
!  AUTHOR : HELIO KOITI KUGA - JUNE 1985
!           INPE - DMC/DDO
!
!  REF. : JACCHIA,L.G.-"ATMOSPHERIC MODELS IN THE
!              REGION FROM 110 TO 2000 KM".IN :
!              CIRA 1972,PART.3,PP. 225-338
!
   IMPLICIT NONE

   real(dp) ain , al , Al10n , am1 , am2 , an , anm , dens , dz , dzx , exp , fact1 , fact2 , gx , gz , Sat3 , sum1 ,  &
            sum2
   real(dp) Tinf , tl1 , Tl2 , tx , z , zd , zend , zr
   INTEGER i , int , j , n
   DIMENSION Al10n(6)

   real(dp),parameter :: ra = 6356.766_dp !! POLAR EARTH RADIUS (KM)
   real(dp),dimension(5),parameter :: wt = 2.0_dp/45.0_dp*[7.0_dp,32.0_dp,12.0_dp,32.0_dp,7.0_dp] !! WEIGHTS FOR THE NEWTON-COTES FIVE POINT QUADRATURE FORMULAE
!
   tx = 371.668_dp + 0.0518806_dp*Tinf - 294.3503_dp*exp(-0.00216222_dp*Tinf)
   gx = 0.054285714_dp*(tx-183.0_dp)
   al = log(Sat3/90.0_dp)
   n = int(al/0.050_dp) + 1
   zr = exp(al/real(n,dp))
   am1 = 28.82678_dp
   tl1 = 183.0_dp
   zend = 90.0_dp
   sum2 = 0.0_dp
   ain = am1*9.534750028_dp/tl1
!
   DO i = 1 , n
      z = zend
      zend = zr*z
      dz = 0.25_dp*(zend-z)
      sum1 = 0.31111111_dp*ain
      DO j = 2 , 5
         z = z + dz
!
!       MOLECULAR WEIGHT FOR Z BETWEEN 90 KM AND
!       100 KM . ACCORDING TO JACCHIA 1971,EQ.1J
!
         zd = z - 90.0_dp
         am2 = 28.82678_dp - 7.40066e-02_dp*zd +                                                                                         &
               zd*(-1.19407e-02_dp*zd+zd*(4.51103e-04_dp*zd+zd*(-8.21895e-06_dp*zd+zd*(1.07561e-05_dp*zd-6.97444e-07_dp*zd*zd))))
!
!       TEMPERATURE FOR Z BETWEEN 90 AND 100 KM
!       EQ. 5R
!
         dzx = z - 125.0_dp
         Tl2 = tx + ((-9.8204695e-06_dp*dzx-7.3039742e-04_dp)*dzx*dzx+1.0_dp)*dzx*gx
         gz = 9.80665_dp*(ra/(ra+z))**2
         ain = am2*gz/Tl2
         sum1 = sum1 + wt(j)*ain
      ENDDO
      sum2 = sum2 + dz*sum1
   ENDDO
!
   fact1 = 0.12027444181_dp
   dens = 3.46e-06_dp*am2*tl1*exp(-fact1*sum2)/am1/Tl2
   anm = 6.02257e+26_dp*dens
   an = anm/am2
   fact2 = anm/28.960_dp
!
   Al10n(1) = log10(0.78110_dp*fact2)
   Al10n(2) = log10(9.3432e-03_dp*fact2)
   Al10n(3) = log10(6.1471e-06_dp*fact2)
   Al10n(4) = log10(1.20955_dp*fact2-an)
   Al10n(5) = log10(2.0_dp*(an-fact2))
   Al10n(6) = Al10n(5) - 15.0_dp
!
END SUBROUTINE stjr01


SUBROUTINE stjr02(Tinf,Sat3,Tz,Dn)
!
!-----
!
!  PURPOSE : STATIC MODEL FOR CALCULATION OF
!            ATMOSPHERIC PROPERTIES FOR THE
!            BAND FROM 100 TO 125 KM.
!
!  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3 ... ALTITUDE               (KM)
!
!  OUTPUTS : TZ ... LOCAL TEMPERATURE (KELVIN)
!            DN ... ALOG10 OF DENSITY NUMBERS (M**-3)
!
!            DN(1) ... N2
!            DN(2) ... A
!            DN(3) ... HE
!            DN(4) ... O2
!            DN(5) ... O
!            DN(6) ... H
!
!  AUTHOR : HELIO KOITI KUGA - JUNE 1985
!           INPE - DMC/DDO
!
!  REF. : KUGA,H.K. "REFORMULACAO COMPUTACIONAL DO
!              MODELO DE JACCHIA-ROBERTS PARA A
!              DENSIDADE ATMOSFERICA".INPE,SAO JOSE
!              DOS CAMPOS,1985.A SER PUBLICADO
!
!         ROBERTS JR,C.E."AN ANALYTICAL MODEL FOR
!              UPPER ATMOSPHERE DENSITIES BASED UPON
!              JACCHIA'S 1970 MODELS."CELESTIAL
!              MECHANICS 4:368-377,1971.
!
   IMPLICIT NONE

   real(dp) abs , am100 , atan , aux , aux1 , aux2 , c0a , cx , d1 , d2 , d3 , d4 , d5 , de100 , deavog , dife , Dn , dpz1 ,  &
        & dpz2
   real(dp) dzx , f3 , f4 , gsubx , h2 , h3 , h4 , prod , pz1 , pz2 , q1 , q2 , q3 , q4 , q5 , q6 , r , r1 , r1n
   real(dp) r2 , r2n , ra , ras , Sat3 , sksf , sksf34 , soma , sqrt , t100 , t100tz , temp , Tinf , tsubx , txmt0 , Tz , ur1 ,      &
        & ur1h2 , ur2 , ur2h3
   real(dp) vra , wr1 , wr2 , x , x2y2 , y , z
   INTEGER i
   DIMENSION Dn(6)
!
!-----
!       R ... UNIVERSAL GAS CONSTANT(JOULES/K MOLE)
!       RA... POLAR EARTH RADIUS    (KM)
!       RAS.. RA**2                 (KM**2)
!
   DATA r/8.31432_dp/
   DATA ra , ras/6356.766_dp , 4.04084739788e+07_dp/
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
   tsubx = 371.6678_dp + 0.0518806_dp*Tinf - 294.3503_dp*exp(-0.00216222_dp*Tinf)
!
!       EQ. 11J
!
   txmt0 = tsubx - 183.0_dp
   gsubx = 0.054285714_dp*txmt0
!
!       VALUE OF SMALL K <= SK  AND SMALL F <= SF
!
   sksf = 9.80665_dp/(r*txmt0)*1500625.0_dp*ras/0.8_dp
!
!       VALUE OF  C0* <= C0A FOR COMPOSING THE
!       FOURTH DEGREE POLYNOMIAL
!
   c0a = -87783750.0_dp + 274614375.0_dp/txmt0
!
!       NEWTON-RAPHSON PROCEDURE FOR OBTAINING
!       THE TWO REAL ROOTS OF THE QUARTIC
!       POLYNOMIAL C4*P(Z) , EQ. 10R
!
!               INITIAL GUESSES
!
   temp = (tsubx-300.0_dp)/200.0_dp
   r1 = 167.77_dp - 3.35_dp*temp
   r2 = 57.34_dp + 7.95_dp*temp
!
   DO i = 1 , 7
      pz1 = c0a + 3542400.0_dp*r1 + r1*(r1*(-52687.5_dp+340.5_dp*r1-0.8_dp*r1*r1))
      pz2 = c0a + 3542400.0_dp*r2 + r2*(r2*(-52687.5_dp+340.5_dp*r2-0.8_dp*r2*r2))
      dpz1 = 3542400.0_dp - 105375.0_dp*r1 + r1*(1021.5_dp*r1-3.2_dp*r1*r1)
      dpz2 = 3542400.0_dp - 105375.0_dp*r2 + r2*(1021.5_dp*r2-3.2_dp*r2*r2)
      r1n = r1 - pz1/dpz1
      r2n = r2 - pz2/dpz2
      IF ( abs(r1n-r1)<1.0e-07_dp .AND. abs(r2n-r2)<1.0e-07_dp ) EXIT
      r1 = r1n
      r2 = r2n
   ENDDO
   r1 = r1n
   r2 = r2n
!
!       COMPLEX ROOTS OR X & X**2+Y**2
!
   soma = r1 + r2
   prod = r1*r2
   dife = r1 - r2
   x = -0.5*(soma-425.625_dp)
   x2y2 = -c0a/(0.8_dp*prod)
!
!       CALCULATE U(R1),U(R2),W(R1),W(R2),CX(CAPITAL X),
!                 AND V(-RA)
!
!  EXPRESSION OF W CORRECTED ACCORDING TO GSFC (NASA,1976)
!
   h2 = r1 + ra
   h3 = r2 + ra
   h4 = ras + 2.0_dp*x*ra + x2y2
!
   ur1h2 = h2*(r1*r1-2.0_dp*x*r1+x2y2)*dife
   ur2h3 = h3*(r2*r2-2.0_dp*x*r2+x2y2)*dife
   wr1 = ra + x2y2/r1
   wr2 = ra + x2y2/r2
   vra = h4*h2*h3
   cx = -h4 - h4
!
   de100 = (((((0.7026942e-32_dp*Tinf*Tinf-0.7734110e-28_dp*Tinf)*Tinf+0.3727894e-24_dp*Tinf)*Tinf-0.1021474e-20_dp*Tinf) &
           *Tinf+0.1711735e-17_dp*Tinf)*Tinf-0.1833490e-14_dp*Tinf+0.1985549e-10_dp)
   t100 = tsubx - 0.94585589_dp*txmt0
   z = Sat3
!
!      NUMBER OF PARTICLES PER M**3 AT 100 KM
!
!         D  ... TOTAL NUMBER
!         D1 ... N2 NITROGEN
!         D2 ... AR ARGON
!         D3 ... HE HELIUM
!         D4 ... O2 DIATOMYC OXYGEN
!         D5 ... O  MONOATOMYC OXYGEN
!
   am100 = 27.6396281382_dp
   deavog = de100*6.02257e+29_dp
   d1 = 0.78110_dp*deavog
   d2 = 0.0093432_dp*deavog
   d3 = 6.1471e-06_dp*deavog
   d4 = (1.20955_dp-28.96_dp/am100)*deavog
   d5 = 2.*(28.96_dp-am100)/am100*deavog
!
!      Q(I) PARAMETERS
!
   ur1 = h2*ur1h2
   ur2 = h3*ur2h3
   q2 = 1/ur1
   q3 = -1/ur2
   q5 = 1/vra
   q4 = (1.0_dp/(prod*ra)+(ra-x2y2/ra)/vra+wr1/ur1h2-wr2/ur2h3)/cx
   q6 = -q5 - 2.0_dp*(x+ra)*q4 + 1.0_dp/ur2h3 - 1.0_dp/ur1h2
   q1 = -q4 - q4 - q3 - q2
!
!       TEMPERATURE FOR Z BETWEEN 100 AND 125 KM
!       EQ. 5R
!
   dzx = z - 125.0_dp
   Tz = tsubx + ((-9.8204695e-06_dp*dzx-7.3039742e-04_dp)*dzx*dzx+1.0_dp)*dzx*gsubx
!
   aux = z - 100.0_dp
   y = sqrt(x2y2-x*x)
   aux1 = z + ra
   aux2 = ra + 100.0_dp
!
   f3 = log(aux1/aux2)*q1 + &
        log((z-r1)/(100.0_dp-r1))*q2 + &
        log((z-r2)/(100.0_dp-r2))*q3 + &
        log((z*z-2.0_dp*x*z+x2y2)/(10000.0_dp-200.0_dp*x+x2y2))*q4
!
   f4 = q5*aux/(aux1*aux2) + q6/y*atan(y*aux/(x2y2+100.0_dp*z-(100.0_dp+z)*x))
!
!      DENSITY NUMBERS D(I) : N2,AR,HE,O2,O,H
!      EQ. 20R
!
   t100tz = t100/Tz
   sksf34 = sksf*(f3+f4)
!
   Dn(1) = log10(d1*t100tz*exp(28.0134_dp*sksf34))
   Dn(2) = log10(d2*t100tz*exp(39.9480_dp*sksf34))
   Dn(3) = log10(d3*t100tz**0.62*exp(4.0026_dp*sksf34))
   Dn(4) = log10(d4*t100tz*exp(31.9988_dp*sksf34))
   Dn(5) = log10(d5*t100tz*exp(15.9994_dp*sksf34))
!
END SUBROUTINE stjr02


SUBROUTINE stjr03(Tinf,Sat3,Tz,Dn)
!
!-----
!
!  PURPOSE : STATIC MODEL FOR CALCULATION OF
!            ATMOSPHERIC PROPERTIES FOR THE
!            BAND ABOVE 125 KM.
!
!  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3 ... ALTITUDE               (KM)
!
!  OUTPUTS : TZ ... LOCAL TEMPERATURE (KELVIN)
!            DN ... ALOG10 OF DENSITY NUMBERS (M**-3)
!
!            DN(1) ... N2
!            DN(2) ... A
!            DN(3) ... HE
!            DN(4) ... O2
!            DN(5) ... O
!            DN(6) ... H
!
!  AUTHOR : HELIO KOITI KUGA - JUNE 1985
!           INPE - DMC/DDO
!
!  REFS. KUGA,H.K."REFORMULACAO COMPUTACIONAL DO
!              MODELO ATMOSFERICO DE JACCHIA-
!              ROBERTS PARA A DENSIDADE ATMOSFERICA".
!              INPE,SAO JOSE DOS CAMPOS,1985.A SER
!              PUBLICADO.
!
!        ROBERTS JR,C.E."AN ANALYTIC MODEL FOR UPPER
!              ATMOSPHERE DENSITIES BASED UPON
!              JACCHIA'S 1970 MODELS."CELESTIAL
!              MECHANICS 4:368-377,1971
!
   IMPLICIT NONE

   real(dp) a1 , a1a2a , a2 , al , aux , aux1 , aux2 , d1 , d2 , d3 , d4 , d5 , d6 , Dn , exp , h500 , r , ra , ras , Sat3
   real(dp) tetx , Tinf , tsubx , txmt0 , Tz , tz500 , z
   INTEGER Icount
   COMMON /ncall / Icount
   DIMENSION Dn(6)
!
!-----
!
!       R....UNIVERSAL GAS CONSTANT (JOULES/K MOLE)
!       RA...POLAR EARTH RADIUS     (KM)
!       RAS..RA**2                  (KM**2)
!
   DATA r/8.31432_dp/
   DATA ra , ras/6356.766_dp , 4.04084739788e+07_dp/
   Icount = Icount + 1
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
   tsubx = 371.6678_dp + 0.0518806_dp*Tinf - 294.3503_dp*exp(-0.00216222_dp*Tinf)
!
   d1 = ((((-0.2296182e-19_dp*Tinf*Tinf+0.1969715e-15_dp*Tinf)*Tinf-0.7139785e-12_dp*Tinf)*Tinf+0.1420228e-08_dp*Tinf)*Tinf-0.1677341e-05_dp*Tinf)&
      & *Tinf + 0.1186783e-02_dp*Tinf + 0.1093155e+02_dp
   d2 = ((((-0.4837461e-19_dp*Tinf*Tinf+0.4127600e-15_dp*Tinf)*Tinf-0.1481702e-11_dp*Tinf)*Tinf+0.2909714e-08_dp*Tinf)*Tinf-0.3391366e-05_dp*Tinf)&
      & *Tinf + 0.2382822e-02_dp*Tinf + 0.8049405e+01_dp
   d3 = (((-0.1270838e-16_dp*Tinf*Tinf+0.9451989e-13_dp*Tinf)*Tinf-0.2894886e-09_dp*Tinf)*Tinf+0.4694319e-06_dp*Tinf) &
      & *Tinf - 0.4383486D-03*Tinf + 0.7646886D+01
   d4 = ((((-0.3131808e-19_dp*Tinf*Tinf+0.2698450e-15_dp*Tinf)*Tinf-0.9782183e-12_dp*Tinf)*Tinf+0.1938454e-08_dp*Tinf)*Tinf-0.2274761e-05_dp*Tinf)&
      & *Tinf + 0.1600311e-02_dp*Tinf + 0.9924237e+01_dp
   d5 = (((0.5116298e-17_dp*Tinf*Tinf-0.3490739e-13_dp*Tinf)*Tinf+0.9239354e-10_dp*Tinf)*Tinf-0.1165003e-06_dp*Tinf) &
      & *Tinf + 0.6118742e-04_dp*Tinf + 0.1097083e+02_dp
   z = Sat3
   al = ((0.2462708e-09_dp*Tinf*Tinf-0.1252487e-05_dp*Tinf)*Tinf+0.1579202e-02_dp*Tinf)*Tinf + 0.2341230e+01_dp*Tinf + 0.1031445e+05_dp
!
!      TEMPERATURE PROFILE EQ. 23R
!
   txmt0 = tsubx - 183.0_dp
   tetx = Tinf - tsubx
   Tz = tetx*exp(-txmt0/tetx*(z-125.0_dp)/35.0_dp*al/(z+ra))
!
!      PARAMETERS G(I) : N2,AR,HE,O2,O EQ. 25'R
!
   aux = 9.80665_dp*ras/(r*al*Tinf)*tetx/txmt0*35.0_dp/6481.766_dp
!
   aux1 = tsubx/(Tinf-Tz)
   aux2 = Tz/tetx
   a1 = log10(aux1)
   a2 = log10(aux2)
   a1a2a = (a1+a2)*aux
!
!        DENSITY NUMBERS D(I) EQ.25R
!
   d1 = d1 + a1 + 28.0134_dp*a1a2a
   d2 = d2 + a1 + 39.9480_dp*a1a2a
   d3 = d3 + 0.62_dp*a1 + 4.0026_dp*a1a2a
   d4 = d4 + a1 + 31.9988_dp*a1a2a
   d5 = d5 + a1 + 15.9994_dp*a1a2a
!
!       CALCULATE TZ(500) EQ.23R
!       T(Z) ALREADY CALCULATED
!
   tz500 = Tinf - tetx*exp(-txmt0/tetx*(375.0_dp/35.0_dp*al/(500.0_dp+ra)))
!
!       INCLUSION OF HYDROGEN
!
!       DENSITY NUMBER FROM EQS. 26R,27R
!
   aux1 = log10(Tinf)
   h500 = 73.13_dp - (39.4_dp-5.5_dp*aux1)*aux1
   a1 = log10(tz500/(Tinf-Tz))
   a2 = log10(Tz/(Tinf-tz500))
   d6 = h500 + a1 + 1.00797_dp*aux*(a1+a2)
!
!        LOAD ALOG10 OF DENSITY NUMBERS
!        IN M**-3
!
   Dn(1) = d1 + 6.0_dp
   Dn(2) = d2 + 6.0_dp
   Dn(3) = d3 + 6.0_dp
   Dn(4) = d4 + 6.0_dp
   Dn(5) = d5 + 6.0_dp
   Dn(6) = d6 + 6.0_dp
   Tz = Tinf - Tz
!
END SUBROUTINE stjr03

FUNCTION temlo(Altu,C)
!
!------
!
! PURPOSE:
!
!     THE FUNCTION TEMLO  EVALUATES  THE  TEMPERATURE AT THE
!     LOCAL REGION.
!
! INPUTS:
!       ALTU     ALTITUDE OF THE POINT, IN KM (90 TO 2000).
!       C        ARRAY THAT CONTAINS THE PARAMETERS USED  IN
!                THE EVALUATION OF EQUATIONS 3 AND 4, AS GI-
!                VEN BY JACCHIA.
!
! OUTPUTS:
!       TEMLO    LOCAL TEMPERATURE AS  DEFINED  IN EQUATIONS
!                3 AND 4, IN KELVIN.
!
! REFERENCES:
!       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
!                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
!                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
!                375).
!
! AUTHOR:    VALDEMIR CARRARA  APR/87    V 1.0
!
   IMPLICIT NONE

   real(dp) Altu , auxi , C , higo , higx , temlo , to , zo , zx

   DIMENSION C(7)
!
!------
!
   DATA zx/125.0_dp/
   DATA zo/90.0_dp/
   DATA to/188.0_dp/

   higx = Altu - zx
   higo = Altu - zo
   temlo = to

   IF ( higo==0.0_dp ) RETURN

   IF ( higx>0.0_dp ) THEN
      temlo = C(7) + C(4)*atan(C(5)*higx+C(6)*higx*higx*higx)
   ELSE
      auxi = higx/higo
      temlo = C(7) + C(1)*atan(C(2)*higx+C(3)*higx*auxi*auxi)
   ENDIF

END FUNCTION temlo


SUBROUTINE sealat(Tyfr,Sudc,Rlat,Altu,Al)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE SEALAT  OBTAINS  THE VARIATIONS  ON THE
!     NUMBER DENSITY OF THE ATMOSPHERE, DUE TO THE SEAZONAL-
!     LATITUDINAL EFFECT.
!     ---
!
! INPUTS:
!
!       TYFR     FRACTION OF THE TROPIC YEAR, IN  THE  RANGE
!                0. TO 1., STARTING ON JAN. 1ST.
!       SUDC     SUN DECLINATION IN RADIANS (-PI TO PI).
!       RLAT     DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
!                POINT, IN RADIANS (-PI TO PI).
!       ALTU     GEOCENTRIC ALTITUDE  IN KM,  IN  THE  RANGE
!                90. TO 2000..
!
! OUTPUTS:
!
!       AL       ARRAY CONTAINING  THE  SEAZONAL-LATITUDINAL
!                VARIATIONS FOR THE HE (HELIUM), O2 (MOLECU-
!                LAR OXYGEN), N2  (MOLECULAR  NITROGEN),  AR
!                (ARGON), O (ATOMIC OXYGEN)  AND  H  (ATOMIC
!                HYDROGEN) NUMBER DENSITY, RESPECTIVELY.
!
!
! REFERENCES:
!
!       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
!                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
!                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
!                375).
!
! AUTHOR:    VALDEMIR CARRARA       APR/87             V 1.0
!            BENTO SILVA DE MATOS
!
   IMPLICIT NONE

   real(dp) Al(6) , Altu , cr(6) , delz , dslm , dslt , esse , pcap , pitw , Rlat , sila , Sudc , Tyfr
   INTEGER i
!
!------
!
   DATA pitw/6.28318530718_dp/

   DATA cr/ - 0.79_dp , 0_dp , 0_dp , 0_dp , -.16_dp , 0_dp/

   sila = sin(Rlat)
   dslt = Sudc*sila/0.409157536545_dp
   delz = Altu - 91.0_dp
   esse = 0.014_dp*delz*exp(-0.0013_dp*delz*delz)
   pcap = sin(pitw*Tyfr+1.72_dp)
   dslm = sign(sila*sila*esse*pcap,Rlat)

   DO i = 1 , 6
      Al(i) = dslt*cr(i) + dslm
   ENDDO

END SUBROUTINE sealat


SUBROUTINE semian(Tyfr,Altu,Alco)
!
!------
!
! PURPOSE:
!
!     THE SUBROUTINE SEMIAN GIVES THE CORRECTION FACTOR ALCO
!     FOR   THE  ATMOSPHERE   NUMBER  DENSITY,  DUE  TO  THE
!     SEMIANNUAL EFFECT.
!     ------
!
! INPUTS:
!
!       TYFR     FRACTION OF THE TROPIC  YEAR,  IN THE RANGE
!                0. TO 1., STARTING ON JAN. 1ST.
!       ALTU     GEOCENTRIC ALTITUDE  IN KM,  IN  THE  RANGE
!                90. TO 2000..
!
! OUTPUTS:
!
!       ALCO     THE SEMIANNUAL VARIATION OF THE  ATMOSPHERE
!                NUMBER DENSITY.
!
! REFERENCES:
!
!       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
!                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
!                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
!                375).
!
! AUTHOR:    VALDEMIR CARRARA        APR/87            V 1.0
!            BENTO SILVA DE MATOS
!
   IMPLICIT NONE

   real(dp) Alco , Altu , auxi , foft , goft , pitw , tauc , Tyfr
!
!------
!
   DATA pitw/6.28318530718_dp/

   auxi = 0.04_dp*Altu*Altu/1.0e+04_dp + 0.05_dp
   foft = auxi*exp(-0.25e-02_dp*Altu)
   auxi = (0.5_dp+0.5_dp*sin(pitw*Tyfr+6.04_dp))**1.65_dp
   tauc = 0.0954_dp*(auxi-0.5_dp) + Tyfr
   auxi = sin(2.0_dp*pitw*tauc+4.26_dp)*(1.0_dp+0.467_dp*sin(pitw*tauc+4.14_dp))
   goft = auxi*0.382_dp + 0.0284_dp
   Alco = foft*goft

END SUBROUTINE semian


!  Replacement for the SOFLUD solar flux library, backed by the new
!  Fortran `space_weather_module` (reads CSSI/Celestrak format files).
!
!  The original SOFLUD is part of an INPE atmospheric model library by
!  Valdemir Carrara; its source is not available here.
!
!  ## Usage
!
!      call soflud_init(filename, status)     ! once, to load the CSSI file
!      call soflud(rjud_1950, dafr, sd, outr) ! then call as before
!
!  ## sd array layout (this implementation)
!
!   - `sd(1..8)` = Kp for 8 3-hour periods (periods 1..8 = 0-3 h, 3-6 h, …, 21-24 h)
!   - `sd(6)`    = 0.0  (time-reference constant; see rdymos note below)
!   - `sd(7)`    = daily average Ap index (for rsmods Ap->Kp conversion)
!   - `sd(9)`    = F10.7 daily observed solar flux
!   - `sd(11)`   = F10.7 81-day centred average
!
!  ## Compatibility notes
!
!  **rsmods path** — uses `sd(7)` (Ap), `sd(9)` (F10.7), `sd(11)` (F10.7a).
!  All three are populated correctly; rsmods works perfectly.
!
!  **rdymos path** — has two issues that preclude direct use:
!
!   1. *Uninitialized variable bug*: `rdymos` contains `rjfl = rjfl - 1.`
!      but `rjfl` is never initialised from `Rjud`. The line should read
!      `rjfl = Rjud - 1.`. Without this fix, `rdymos` produces garbage dates.
!
!   2. *Kp period conflict*: `rdymos` selects Kp via
!      `nd = int((Dafr/3600 - sd(6) + 12 - 6.696) / 3)` then `sf(3) = sd(nd)`.
!      With `sd(6) = 0`, `nd` maps UT 0-3 h -> 1, 3-6 h -> 2, …, 21-24 h -> 8.
!      However:
!      - `nd = 6` at UT 15-18 h -> `sd(6) = 0.0` (wrong; Kp treated as 0)
!      - `nd = 7` at UT 18-21 h -> `sd(7) = Ap` (wrong unit)
!      For all other UT hours (18 of 24) the correct Kp period is returned.
!
!  **Recommended alternative for rdymos**: use `rdymos_cssi` (below), which
!  calls `dyjrmo` directly with F10.7 and Kp obtained from `prepare_flux_data`.
!  This avoids both issues and matches the timing logic of the new Fortran model.

!---------------------------------------------------------------------------
!> Initialize the soflud wrapper by loading a CSSI space weather file.
!>
!> Must be called once before any call to `soflud` or `rdymos_cssi`.
!>
!> @param[in]  filename  Path to the CSSI space weather file
!> @param[out] status    0 = success, non-zero = error

   subroutine soflud_init(filename, status)
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: status

      call sw_global%initialize(filename, status)
      if (status == 0) sw_initialized = .true.

   end subroutine soflud_init

!---------------------------------------------------------------------------
!>
!  Drop-in replacement for the missing SOFLUD library subroutine.
!
!  Retrieves solar flux and geomagnetic data from the CSSI space weather
!  file loaded by `soflud_init`.

   subroutine soflud(rjud_1950, dafr, sd, outr)
      real(dp), intent(in)  :: rjud_1950 !! Modified Julian Date referenced to 1950.0
                                         !! (JD − 2433282.5).  Standard MJD = rjud_1950 + 33282.
      real(dp), intent(in)  :: dafr !! Time of day in seconds (not used for data lookup;
                                    !! retained for interface compatibility).
      real(dp), intent(out) :: sd(15) !! 15-element output array (see module header).
      real(dp), intent(out) :: outr !! Status: 0.0 = success, non-zero = error code.

      real(dp) :: mjd
      type(flux_data_type) :: flux_data
      logical :: status

      sd   = 0.0_dp
      outr = 0.0_dp

      if (.not. sw_initialized) then
         write(*,'(A)') 'ERROR (soflud): not initialised — call soflud_init first.'
         outr = 1.0_dp
         return
      end if

      ! Convert MJD-1950 to standard MJD
      mjd = rjud_1950 + MJD_1950_OFFSET

      call sw_global%get_flux_data(mjd, flux_data, status)
      if (.not. status) then
         outr = 1.0_dp
         write(*,'(A)') 'ERROR (soflud): space weather data lookup failed.'
         return
      end if

      ! Kp for 8 3-hour periods (sd(6) and sd(7) are overridden below)
      sd(1:8) = flux_data%kp(1:8)

      ! sd(6) = 0.0 — time-reference constant used by rdymos nd formula.
      ! This overwrites the period-6 (15-18 h) Kp slot; see module header.
      sd(6)  = 0.0_dp

      ! sd(7) = daily average Ap — used by rsmods for Ap->Kp conversion.
      ! This overwrites the period-7 (18-21 h) Kp slot; see module header.
      sd(7)  = flux_data%ap_avg

      ! F10.7 solar flux
      sd(9)  = flux_data%f107_obs

      ! F10.7 81-day centred average
      sd(11) = flux_data%f107a_obs_ctr

   end subroutine soflud

   !---------------------------------------------------------------------------
   !>
   !  Drop-in replacement for `rdymos` in `roberts.f90`, using the new
   !  `space_weather_module` instead of the unavailable SOFLUD library.
   !
   !  Fixes two bugs in the original `rdymos`:
   !   1. Uninitialised `rjfl` (should be `rjfl = Rjud - 1.`)
   !   2. Unreliable Kp period selection via `sd(nd)` (see module header)
   !
   !  F10.7 and Kp are obtained through `prepare_flux_data`, applying the
   !  same timing as the new model (6.7 h Kp lag, previous-day F10.7,
   !  detected-day F10.7a).  The actual computation is delegated to `rsdamo`
   !  (in roberts.f90), which calls `dyjrmo` and reorders the output densities.
   !
   !  This subroutine has the **same interface as `rdymos`** and can be used
   !  as a direct replacement in any code that calls `rdymos`.

   subroutine rdymos_cssi(Sa, Su, Rjud, Dafr, Gsti, Te, Ad, Wmol, Rhod, status)

      real(dp), intent(in)  :: Sa(3) !! Sa(1)=RA (rad), Sa(2)=geocentric lat (rad), Sa(3)=alt (m)
      real(dp), intent(in)  :: Su(2) !! Su(1)=sun RA (rad), Su(2)=sun dec (rad)
      real(dp), intent(in)  :: Rjud !! Modified Julian Date referred to 1950.0 (JD − 2433282.5)
      real(dp), intent(in)  :: Dafr !! Time of day UT in seconds (0-86400)
      real(dp), intent(in)  :: Gsti !! Greenwich Sidereal Time in radians (not used; compatibility only)
      real(dp), intent(out) :: Te(2) !! Te(1)=T∞ (K), Te(2)=local T (K)
      real(dp), intent(in)  :: Ad(6) !! log10 number densities (m^-3):
                                     !! Ad(1)=He, Ad(2)=O₂, Ad(3)=N₂, Ad(4)=Ar, Ad(5)=O, Ad(6)=H
      real(dp), intent(in)  :: Wmol !! Mean molecular weight (kg/kgmol)
      real(dp), intent(in)  :: Rhod !! Mass density (kg/m^3)
      integer, intent(out) :: status !! 0 = success, non-zero = space weather lookup error

      type(flux_data_type) :: flux_data
      real(dp) :: utc_mjd, kp, f107, f107a
      real(dp)  :: sf(3)
      logical  :: sw_status

      if (.not. sw_initialized) then
         write(*,'(A)') 'ERROR (rdymos_cssi): not initialised — call soflud_init first.'
         status = 1
         return
      end if

      ! Convert from MJD-1950 + fractional day to standard MJD
      utc_mjd = Rjud + MJD_1950_OFFSET + Dafr / 86400.0_dp

      ! Retrieve today's flux record
      call sw_global%get_flux_data(utc_mjd, flux_data, sw_status)
      if (sw_status) then
         status = 0
      else
         write(*,'(A)') 'ERROR (rdymos_cssi): space weather data lookup failed.'
         status = 1
         return
      end if

      ! Apply timing offsets: 6.7 h Kp lag, previous-day F10.7, detected-day F10.7a
      call prepare_flux_data(sw_global, flux_data, utc_mjd, kp, f107, f107a)

      ! Pack SF array for rsdamo (same layout as soflud output to rdymos)
      sf(1) = f107     ! F10.7 daily (1.71-day lag)
      sf(2) = f107a    ! F10.7 81-day average
      sf(3) = kp       ! Kp (6.7-h lag)

      ! Delegate computation to rsdamo, which calls dyjrmo and reorders output
      call rsdamo(Sa, Su, sf, Rjud, Dafr, Gsti, Te, Ad, Wmol, Rhod)

      status = 0

   end subroutine rdymos_cssi

end module inpe_roberts_module