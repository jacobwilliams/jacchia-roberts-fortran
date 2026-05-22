!>
!  Test code for the Jacchia-Roberts model.
!  Uses the older INPE code for the Jacchia-Roberts model
!  as a reference for testing the new implementation.

module inpe_roberts_module

   use space_weather_module,       only: sw_data_type, flux_data_type
   use jacchia_roberts_module,     only: prepare_flux_data
   use jacchia_roberts_kinds,      only: dp

   implicit none

   !> Standard MJD = MJD-1950 + 33282  (JD-2400000.5 = (JD-2433282.5) + 33282)
   real(dp), parameter :: MJD_1950_OFFSET = 33282.0_dp

   type(sw_data_type), save :: sw_global
   logical,            save :: sw_initialized = .false.

   public :: soflud_init
   public :: soflud
   public :: rdymos_cssi

   contains

!     LIBRARY DENSITY
!
!----
!
!     The library DENSITY includes several routines to
!     compute the high atmospheric properties, using the
!     Robert's version of the Jacchia 1970 model.
!
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
! SUBCALLS:
!
!     SOFLUD
!     RSDAMO
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

   REAL*8 Ad , dafl , Dafr , Gsti , outr , Rhod , rjfl , Rjud , Sa , sd , sf , Su , tauo , Te , Wmol
   INTEGER int , nd
   DIMENSION Sa(3) , Su(2) , Te(2) , Ad(6)
   DIMENSION sf(3) , sd(15)
!
!------
!

   rjfl = rjfl - 1.
   IF ( Dafr<61344. ) THEN
      rjfl = rjfl - 1.
      dafl = Dafr + 25056.
   ELSE
      dafl = Dafr - 61344.
   ENDIF

   CALL soflud(rjfl,dafl,sd,outr)

   IF ( outr/=0. ) THEN
      WRITE (6,*) ' ERROR IN ROUTINE JDYMOS: OUTR = ' , int(outr)
      STOP
   ENDIF

   tauo = 6.696
   nd = (Dafr/3600.-sd(6)+12.-tauo)/3.
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

   REAL*8 Al , Altu , Dafr , dlog , dsqrt , outr , Rhod , Rjud , sd , sf , Te , vari , Wmol
   INTEGER int
   DIMENSION Te(2) , Al(6)
   DIMENSION sf(3) , sd(15)
!
!------
!
   CALL soflud(Rjud,Dafr,sd,outr)

   IF ( outr/=0. ) THEN
      WRITE (6,*) ' ERROR IN ROUTINE ISMODS: OUTR = ' , int(outr)
      STOP
   ENDIF

   sf(1) = sd(9)
   sf(2) = sd(11)
   vari = .154*sd(7)
   sf(3) = 1.89*dlog(vari+dsqrt(vari*vari+1.D0))

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

   REAL*8 Ad , al , amjd , Dafr , Gsti , Rhod , Rjud , Sa , Sf , Su , Te , Wmol

   DIMENSION Sa(3) , Su(2) , Sf(3) , Te(2) , Ad(6)
   DIMENSION al(6)

!
!------
!

   amjd = Rjud + 33282. + Dafr/86400.

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

   REAL*8 Ad , al , Altu , anac , anut , avog , fbar , flux , heig , Rhod , Sf , Te , thaf , tz , weig , wm , Wmol
   INTEGER ic

   DIMENSION Sf(3) , Te(2) , Ad(6)
   DIMENSION wm(6) , al(6)
!
!------
!
   DATA avog/6.02217D+26/
   DATA wm/4.0026 , 31.9988 , 28.0134 , 39.948 , 15.9994 , 1.00797/

   flux = Sf(1)
   fbar = Sf(2)
   thaf = 379.0 + 3.24*fbar + 1.3*(flux-fbar)
   heig = Altu/1.D3

   CALL stjrmo(thaf,heig,tz,al)

   Ad(1) = al(3)
   Ad(2) = al(4)
   Ad(3) = al(1)
   Ad(4) = al(2)
   Ad(5) = al(5)
   Ad(6) = al(6)
   anut = 0.
   weig = 0.

   DO ic = 1 , 6
      anac = 10.**Ad(ic)
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

   REAL*8 Ad , al , anac , anut , avog , Heig , Rhod , Tinf , tz , weig , wm , Wmol
   INTEGER ic

   DIMENSION Ad(6)
   DIMENSION al(6) , wm(6)
!
!------
!
   DATA avog/6.02217D+26/
   DATA wm/4.0026 , 31.9988 , 28.0134 , 39.948 , 15.9994 , 1.00797/

   CALL stjrmo(Tinf,Heig,tz,al)

   Ad(1) = al(3)
   Ad(2) = al(4)
   Ad(3) = al(1)
   Ad(4) = al(2)
   Ad(5) = al(5)
   Ad(6) = al(6)
   anut = 0.
   weig = 0.

   DO ic = 1 , 6
      anac = 10.**Ad(ic)
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

   REAL*8 abs , Amw , capphi , cons25 , cos , d1 , d2 , d3 , d4 , d5 , d6 , Dens , df , Djm , dlhe , dlr , dlr20 , dlrgm , dlrsa , &
        & dlrsl
   REAL*8 Dn , dtg , dtg18 , dtg20 , eta , exp , expkp , f , fdfz , fs , fsm , gdft , Geo , h , pid4 , piv2 , piv4 , pk , s , Sat
   REAL*8 sat1 , sat2 , sat3 , sign , sin , sumn , sumnm , Sun , sun1 , sun2 , tanh , tau , Temp , theta , tinf , tsubc , tsubl ,  &
        & tz
   INTEGER mod
   DIMENSION Sun(2) , Sat(3) , Geo(3) , Temp(2) , Dn(6)
!-----
   DATA piv2 , piv4 , pid4 , cons25/6.2831853D0 , 12.566371D0 , 0.78539816D0 , 0.35355339D0/
!
!      PIV2 = 2 * PI
!      PIV4 = 4 * PI
!      PID4 = PI / 4
!      CONS25 = SIN (PI/4) **3
!
   sun1 = Sun(1)
   sun2 = Sun(2)
   sat1 = Sat(1)
   sat2 = Sat(2)
   sat3 = Sat(3)/1000.
   fs = Geo(1)
   fsm = Geo(2)
   pk = Geo(3)
!
!       MINIMUM NIGHT-TIME TEMPERATURE OF THE GLOBAL
!       EXOSPHERIC TEMPERATURE DISTRIBUTION WHEN THE
!       GEOMAGNETIC ACTIVITY INDEX KP = 0
!       EQUATION 14J
!
   tsubc = 379. + 3.24*fsm + 1.3*(fs-fsm)
!
!       EQUATION 15J
!
   eta = 0.5*abs(sat2-sun2)
   theta = 0.5*abs(sat2+sun2)
!
!       EQUATION 16J
!
   h = sat1 - sun1
   tau = h - 0.64577182 + 0.10471976*sin(h+0.75049158)
!
!       EXOSPHERIC TEMPERATURE TSUBL WITHOUT CORRECTION
!       FOR GEOMAGNETIC ACTIVITY
!       EQUATION 17J
!
   s = sin(theta)**2.2
   df = s + (cos(eta)**2.2D0-s)*abs(cos(0.5*tau))**3
   tsubl = tsubc*(1.+0.3*df)
!
!       EQUATION 18J
!
   expkp = exp(pk)
   dtg18 = 28.*pk + 0.03*expkp
!
!       EQUATION 20J
!
   dtg20 = 14.*pk + 0.02*expkp
   dlr20 = 0.012*pk + 1.2D-05*expkp
!
!       THE FOLLOWING STATEMENTS EFFECT A CONTINUOUS
!       TRANSITION FROM EQ. 20J AT HEIGHTS WELL BELOW
!       350 KM TO EQ. 18J AT HEIGHTS WELL ABOVE
!       350 KM .
!
   f = 0.5*(tanh(0.04*(sat3-350.))+1.)
   dlrgm = dlr20*(1.-f)
   dtg = dtg20*(1.-f) + dtg18*f
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
   capphi = mod((Djm-36204.)/365.2422,1.D0)
!
!   EQ. 22J
!
   tau = capphi + 0.09544*((0.5+0.5*sin(piv2*capphi+6.035))**1.650-0.5)
   gdft = 0.02835 + 0.3817*(1.+0.4671*sin(piv2*tau+4.137))*sin(piv4*tau+4.259)
   fdfz = (5.876D-07*sat3**2.331D0+0.06328)*exp(-2.868D-03*sat3)
!
!   EQ. 21J  SEMI-ANNUAL VARIATION
!
   dlrsa = fdfz*gdft
!
!   EQ. 24J  SEASONAL-LATITUDINAL VARIATION OF THE
!            LOWER THERMOSPHERE
!
   dlrsl = 0.014*(sat3-90.)*exp(-0.0013*(sat3-90.)**2)*sign(1.D0,sat2)*sin(piv2*capphi+1.72)*sin(sat2)**2
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
   dlhe = 0.65*abs(sun2/0.4091609)*(sin(pid4-0.5*sat2*sign(1.D0,sun2))**3-cons25)
   Dn(3) = Dn(3) + dlhe
!
!  COMPUTE DENSITY AND MEAN MOLECULAR WEIGHT
!
   d1 = 10.**Dn(1)
   d2 = 10.**Dn(2)
   d3 = 10.**Dn(3)
   d4 = 10.**Dn(4)
   d5 = 10.**Dn(5)
   d6 = 10.**Dn(6)
   sumn = d1 + d2 + d3 + d4 + d5 + d6
   sumnm = 28.0134*d1 + 39.9480*d2 + 4.0026*d3 + 31.9988*d4 + 15.9994*d5 + 1.00797*d6
   Amw = sumnm/sumn
   Dens = sumnm/6.02257D+26
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

   REAL*8 Dn , Sat3 , Tinf , Tz
   DIMENSION Dn(6)
!-----
   IF ( Sat3>125. ) THEN
!
      CALL stjr03(Tinf,Sat3,Tz,Dn)
      RETURN
   ELSEIF ( Sat3>100. ) THEN
!
      CALL stjr02(Tinf,Sat3,Tz,Dn)
      RETURN
   ELSEIF ( Sat3<90. ) THEN
!
      PRINT 99001
99001 FORMAT (1X,'ATENCA0 : MENSAGEM DA ROTINA DE ',/,1X,'*******   CALCULO DA DENSIDADE AT-',/,1X,'          MOSFERICA.',//,3X,   &
             &'ALTITUDE DO SATELITE MENOR QUE 90 KM')
      RETURN
   ENDIF
   CALL stjr01(Tinf,Sat3,Tz,Dn)
   RETURN
!
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

   REAL*8 ain , al , Al10n , am1 , am2 , an , anm , dens , dfloat , dz , dzx , exp , fact1 , fact2 , gx , gz , ra , Sat3 , sum1 ,  &
        & sum2
   REAL*8 Tinf , tl1 , Tl2 , tx , wt , z , zd , zend , zr
   INTEGER i , int , j , log , log10 , n
   DIMENSION wt(5) , Al10n(6)
!
!-----
!
!   RA = POLAR EARTH RADIUS (KM)
!   WT = WEIGHTS FOR THE NEWTON-COTES
!        FIVE POINT QUADRATURE FORMULAE
!
   DATA ra , wt/6356.766 , 0.31111111 , 1.4222222 , 0.53333333 , 1.4222222 , 0.31111111/
!
   tx = 371.668 + 0.0518806*Tinf - 294.3503*exp(-0.00216222*Tinf)
   gx = 0.054285714*(tx-183.)
   al = log(Sat3/90.)
   n = int(al/0.050) + 1
   zr = exp(al/dfloat(n))
   am1 = 28.82678
   tl1 = 183.
   zend = 90.
   sum2 = 0.
   ain = am1*9.534750028/tl1
!
   DO i = 1 , n
      z = zend
      zend = zr*z
      dz = 0.25*(zend-z)
      sum1 = 0.31111111*ain
      DO j = 2 , 5
         z = z + dz
!
!       MOLECULAR WEIGHT FOR Z BETWEEN 90 KM AND
!       100 KM . ACCORDING TO JACCHIA 1971,EQ.1J
!
         zd = z - 90.
         am2 = 28.82678 - 7.40066D-02*zd +                                                                                         &
             & zd*(-1.19407D-02*zd+zd*(4.51103D-04*zd+zd*(-8.21895D-06*zd+zd*(1.07561D-05*zd-6.97444D-07*zd*zd))))
!
!       TEMPERATURE FOR Z BETWEEN 90 AND 100 KM
!       EQ. 5R
!
         dzx = z - 125.
         Tl2 = tx + ((-9.8204695D-06*dzx-7.3039742D-04)*dzx*dzx+1.)*dzx*gx
         gz = 9.80665*(ra/(ra+z))**2
         ain = am2*gz/Tl2
         sum1 = sum1 + wt(j)*ain
      ENDDO
      sum2 = sum2 + dz*sum1
   ENDDO
!
   fact1 = 0.12027444181
   dens = 3.46D-06*am2*tl1*exp(-fact1*sum2)/am1/Tl2
   anm = 6.02257D+26*dens
   an = anm/am2
   fact2 = anm/28.960
!
   Al10n(1) = log10(0.78110*fact2)
   Al10n(2) = log10(9.3432D-03*fact2)
   Al10n(3) = log10(6.1471D-06*fact2)
   Al10n(4) = log10(1.20955*fact2-an)
   Al10n(5) = log10(2.*(an-fact2))
   Al10n(6) = Al10n(5) - 15.
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

   REAL*8 abs , am100 , atan , aux , aux1 , aux2 , c0a , cx , d1 , d2 , d3 , d4 , d5 , de100 , deavog , dife , dlog , Dn , dpz1 ,  &
        & dpz2
   REAL*8 dzx , exp , f3 , f4 , gsubx , h2 , h3 , h4 , prod , pz1 , pz2 , q1 , q2 , q3 , q4 , q5 , q6 , r , r1 , r1n
   REAL*8 r2 , r2n , ra , ras , Sat3 , sksf , sksf34 , soma , sqrt , t100 , t100tz , temp , Tinf , tsubx , txmt0 , Tz , ur1 ,      &
        & ur1h2 , ur2 , ur2h3
   REAL*8 vra , wr1 , wr2 , x , x2y2 , y , z
   INTEGER i , log , log10
   DIMENSION Dn(6)
!
!-----
!       R ... UNIVERSAL GAS CONSTANT(JOULES/K MOLE)
!       RA... POLAR EARTH RADIUS    (KM)
!       RAS.. RA**2                 (KM**2)
!
   DATA r/8.31432D0/
   DATA ra , ras/6356.766D0 , 4.04084739788D+07/
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
   tsubx = 371.6678 + 0.0518806*Tinf - 294.3503*exp(-0.00216222*Tinf)
!
!       EQ. 11J
!
   txmt0 = tsubx - 183.
   gsubx = 0.054285714*txmt0
!
!       VALUE OF SMALL K <= SK  AND SMALL F <= SF
!
   sksf = 9.80665/(r*txmt0)*1500625.*ras/0.8
!
!       VALUE OF  C0* <= C0A FOR COMPOSING THE
!       FOURTH DEGREE POLYNOMIAL
!
   c0a = -87783750. + 274614375./txmt0
!
!       NEWTON-RAPHSON PROCEDURE FOR OBTAINING
!       THE TWO REAL ROOTS OF THE QUARTIC
!       POLYNOMIAL C4*P(Z) , EQ. 10R
!
!               INITIAL GUESSES
!
   temp = (tsubx-300.)/200.
   r1 = 167.77 - 3.35*temp
   r2 = 57.34 + 7.95*temp
!
   SPAG_Loop_1_1: DO i = 1 , 7
      pz1 = c0a + 3542400.*r1 + r1*(r1*(-52687.5+340.5*r1-0.8*r1*r1))
      pz2 = c0a + 3542400.*r2 + r2*(r2*(-52687.5+340.5*r2-0.8*r2*r2))
      dpz1 = 3542400. - 105375.*r1 + r1*(1021.5*r1-3.2*r1*r1)
      dpz2 = 3542400. - 105375.*r2 + r2*(1021.5*r2-3.2*r2*r2)
      r1n = r1 - pz1/dpz1
      r2n = r2 - pz2/dpz2
      IF ( abs(r1n-r1)<1.D-07 .AND. abs(r2n-r2)<1.D-07 ) EXIT SPAG_Loop_1_1
      r1 = r1n
      r2 = r2n
   ENDDO SPAG_Loop_1_1
   r1 = r1n
   r2 = r2n
!
!       COMPLEX ROOTS OR X & X**2+Y**2
!
   soma = r1 + r2
   prod = r1*r2
   dife = r1 - r2
   x = -0.5*(soma-425.625)
   x2y2 = -c0a/(0.8*prod)
!
!       CALCULATE U(R1),U(R2),W(R1),W(R2),CX(CAPITAL X),
!                 AND V(-RA)
!
!  EXPRESSION OF W CORRECTED ACCORDING TO GSFC (NASA,1976)
!
   h2 = r1 + ra
   h3 = r2 + ra
   h4 = ras + 2.*x*ra + x2y2
!
   ur1h2 = h2*(r1*r1-2.*x*r1+x2y2)*dife
   ur2h3 = h3*(r2*r2-2.*x*r2+x2y2)*dife
   wr1 = ra + x2y2/r1
   wr2 = ra + x2y2/r2
   vra = h4*h2*h3
   cx = -h4 - h4
!
   de100 = (((((0.7026942D-32*Tinf*Tinf-0.7734110D-28*Tinf)*Tinf+0.3727894D-24*Tinf)*Tinf-0.1021474D-20*Tinf)                      &
         & *Tinf+0.1711735D-17*Tinf)*Tinf-0.1833490D-14*Tinf+0.1985549D-10)
   t100 = tsubx - 0.94585589*txmt0
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
   am100 = 27.6396281382
   deavog = de100*6.02257D+29
   d1 = 0.78110*deavog
   d2 = 0.0093432*deavog
   d3 = 6.1471D-06*deavog
   d4 = (1.20955-28.96/am100)*deavog
   d5 = 2.*(28.96-am100)/am100*deavog
!
!      Q(I) PARAMETERS
!
   ur1 = h2*ur1h2
   ur2 = h3*ur2h3
   q2 = 1/ur1
   q3 = -1/ur2
   q5 = 1/vra
   q4 = (1./(prod*ra)+(ra-x2y2/ra)/vra+wr1/ur1h2-wr2/ur2h3)/cx
   q6 = -q5 - 2.*(x+ra)*q4 + 1./ur2h3 - 1./ur1h2
   q1 = -q4 - q4 - q3 - q2
!
!       TEMPERATURE FOR Z BETWEEN 100 AND 125 KM
!       EQ. 5R
!
   dzx = z - 125.
   Tz = tsubx + ((-9.8204695D-06*dzx-7.3039742D-04)*dzx*dzx+1.)*dzx*gsubx
!
   aux = z - 100.
   y = sqrt(x2y2-x*x)
   aux1 = z + ra
   aux2 = ra + 100.
!
   f3 = dlog(aux1/aux2)*q1 + log((z-r1)/(100.-r1))*q2 + log((z-r2)/(100.-r2))*q3 + log((z*z-2.*x*z+x2y2)/(10000.-200.*x+x2y2))*q4
!
   f4 = q5*aux/(aux1*aux2) + q6/y*atan(y*aux/(x2y2+100.*z-(100.+z)*x))
!
!      DENSITY NUMBERS D(I) : N2,AR,HE,O2,O,H
!      EQ. 20R
!
   t100tz = t100/Tz
   sksf34 = sksf*(f3+f4)
!
   Dn(1) = log10(d1*t100tz*exp(28.0134*sksf34))
   Dn(2) = log10(d2*t100tz*exp(39.9480*sksf34))
   Dn(3) = log10(d3*t100tz**0.62*exp(4.0026*sksf34))
   Dn(4) = log10(d4*t100tz*exp(31.9988*sksf34))
   Dn(5) = log10(d5*t100tz*exp(15.9994*sksf34))
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

   REAL*8 a1 , a1a2a , a2 , al , aux , aux1 , aux2 , d1 , d2 , d3 , d4 , d5 , d6 , Dn , exp , h500 , r , ra , ras , Sat3
   REAL*8 tetx , Tinf , tsubx , txmt0 , Tz , tz500 , z
   INTEGER Icount , log10
   COMMON /ncall / Icount
   DIMENSION Dn(6)
!
!-----
!
!       R....UNIVERSAL GAS CONSTANT (JOULES/K MOLE)
!       RA...POLAR EARTH RADIUS     (KM)
!       RAS..RA**2                  (KM**2)
!
   DATA r/8.31432D0/
   DATA ra , ras/6356.766D0 , 4.04084739788D+07/
   Icount = Icount + 1
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
   tsubx = 371.6678 + 0.0518806*Tinf - 294.3503*exp(-0.00216222*Tinf)
!
   d1 = ((((-0.2296182D-19*Tinf*Tinf+0.1969715D-15*Tinf)*Tinf-0.7139785D-12*Tinf)*Tinf+0.1420228D-08*Tinf)*Tinf-0.1677341D-05*Tinf)&
      & *Tinf + 0.1186783D-02*Tinf + 0.1093155D+02
   d2 = ((((-0.4837461D-19*Tinf*Tinf+0.4127600D-15*Tinf)*Tinf-0.1481702D-11*Tinf)*Tinf+0.2909714D-08*Tinf)*Tinf-0.3391366D-05*Tinf)&
      & *Tinf + 0.2382822D-02*Tinf + 0.8049405D+01
   d3 = (((-0.1270838D-16*Tinf*Tinf+0.9451989D-13*Tinf)*Tinf-0.2894886D-09*Tinf)*Tinf+0.4694319D-06*Tinf)                          &
      & *Tinf - 0.4383486D-03*Tinf + 0.7646886D+01
   d4 = ((((-0.3131808D-19*Tinf*Tinf+0.2698450D-15*Tinf)*Tinf-0.9782183D-12*Tinf)*Tinf+0.1938454D-08*Tinf)*Tinf-0.2274761D-05*Tinf)&
      & *Tinf + 0.1600311D-02*Tinf + 0.9924237D+01
   d5 = (((0.5116298D-17*Tinf*Tinf-0.3490739D-13*Tinf)*Tinf+0.9239354D-10*Tinf)*Tinf-0.1165003D-06*Tinf)                           &
      & *Tinf + 0.6118742D-04*Tinf + 0.1097083D+02
   z = Sat3
   al = ((0.2462708D-09*Tinf*Tinf-0.1252487D-05*Tinf)*Tinf+0.1579202D-02*Tinf)*Tinf + 0.2341230D+01*Tinf + 0.1031445D+05
!
!      TEMPERATURE PROFILE EQ. 23R
!
   txmt0 = tsubx - 183.
   tetx = Tinf - tsubx
   Tz = tetx*exp(-txmt0/tetx*(z-125.)/35.*al/(z+ra))
!
!      PARAMETERS G(I) : N2,AR,HE,O2,O EQ. 25'R
!
   aux = 9.80665*ras/(r*al*Tinf)*tetx/txmt0*35./6481.766
!
   aux1 = tsubx/(Tinf-Tz)
   aux2 = Tz/tetx
   a1 = log10(aux1)
   a2 = log10(aux2)
   a1a2a = (a1+a2)*aux
!
!        DENSITY NUMBERS D(I) EQ.25R
!
   d1 = d1 + a1 + 28.0134*a1a2a
   d2 = d2 + a1 + 39.9480*a1a2a
   d3 = d3 + 0.62*a1 + 4.0026*a1a2a
   d4 = d4 + a1 + 31.9988*a1a2a
   d5 = d5 + a1 + 15.9994*a1a2a
!
!       CALCULATE TZ(500) EQ.23R
!       T(Z) ALREADY CALCULATED
!
   tz500 = Tinf - tetx*exp(-txmt0/tetx*(375./35.*al/(500.+ra)))
!
!       INCLUSION OF HYDROGEN
!
!       DENSITY NUMBER FROM EQS. 26R,27R
!
   aux1 = log10(Tinf)
   h500 = 73.13 - (39.4-5.5*aux1)*aux1
   a1 = log10(tz500/(Tinf-Tz))
   a2 = log10(Tz/(Tinf-tz500))
   d6 = h500 + a1 + 1.00797*aux*(a1+a2)
!
!        LOAD ALOG10 OF DENSITY NUMBERS
!        IN M**-3
!
   Dn(1) = d1 + 6.
   Dn(2) = d2 + 6.
   Dn(3) = d3 + 6.
   Dn(4) = d4 + 6.
   Dn(5) = d5 + 6.
   Dn(6) = d6 + 6.
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

   REAL*8 Altu , auxi , C , datan , higo , higx , temlo , to , zo , zx

   DIMENSION C(7)
!
!------
!
   DATA zx/125.D0/
   DATA zo/90.0D0/
   DATA to/188.D0/

   higx = Altu - zx
   higo = Altu - zo
   temlo = to

   IF ( higo==0.D0 ) RETURN

   IF ( higx>0.D0 ) THEN
      temlo = C(7) + C(4)*datan(C(5)*higx+C(6)*higx*higx*higx)
   ELSE
      auxi = higx/higo
      temlo = C(7) + C(1)*datan(C(2)*higx+C(3)*higx*auxi*auxi)
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

   REAL*8 Al , Altu , cr , delz , dexp , dsign , dsin , dslm , dslt , esse , pcap , pitw , Rlat , sila , Sudc , Tyfr
   INTEGER i

   DIMENSION Al(6)
   DIMENSION cr(6)
!
!------
!
   DATA pitw/6.28318530718D0/

   DATA cr/ - 0.79D0 , 0.D0 , 0.D0 , 0.D0 , -.16D0 , 0.D0/

   sila = dsin(Rlat)
   dslt = Sudc*sila/0.409157536545D0
   delz = Altu - 91.D0
   esse = 0.014D0*delz*dexp(-0.0013D0*delz*delz)
   pcap = dsin(pitw*Tyfr+1.72D0)
   dslm = dsign(sila*sila*esse*pcap,Rlat)

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

   REAL*8 Alco , Altu , auxi , dexp , dsin , foft , goft , pitw , tauc , Tyfr
!
!------
!
   DATA pitw/6.28318530718D0/

   auxi = 0.04D0*Altu*Altu/1.D+04 + 0.05D0
   foft = auxi*dexp(-0.25D-02*Altu)
   auxi = (0.5D0+0.5D0*dsin(pitw*Tyfr+6.04D0))**1.65D0
   tauc = 0.0954D0*(auxi-0.5D0) + Tyfr
   auxi = dsin(2.D0*pitw*tauc+4.26D0)*(1.D0+0.467D0*dsin(pitw*tauc+4.14D0))
   goft = auxi*0.382D0 + 0.0284D0
   Alco = foft*goft

END SUBROUTINE semian

! FUNCTION datanh(X)
! !
! !------
! !
! ! PURPOSE:
! !
! !     THE DATANH FUNCTION CALCULATES THE HYPERBOLIC  TANGENT
! !     ARC OF AN ARGUMENT X, IN DOUBLE PRECISION.
! !
! ! INPUTS:
! !       X        HYPERBOLIC TANGENT VALUE (-1<X<1).
! !
! ! OUTPUTS:
! !       DATANH   HYPERBOLIC TANGENT ARC.
! !
! ! AUTHOR:    VALDEMIR CARRARA  APR/87    V 1.0
! !
!    IMPLICIT NONE

!    REAL*8 argu , dabs , datanh , X
!    INTEGER log
! !
! !------
! !
!    IF ( dabs(X)>1.D0 ) THEN
!       WRITE (*,*) 'INVALID HYPERBOLIC TANGENT ARC ARGUMENT'
!       STOP
!    ENDIF

!    argu = (1.D0+X)/(1.D0-X)

!    datanh = log(argu)/2.D0

! END FUNCTION datanh



!> @file soflud_module.f90
!>
!> Replacement for the SOFLUD solar flux library, backed by the new
!> Fortran `space_weather_module` (reads CSSI/Celestrak format files).
!>
!> The original SOFLUD is part of an INPE atmospheric model library by
!> Valdemir Carrara; its source is not available here.
!>
!> ## Usage
!>
!>     call soflud_init(filename, status)     ! once, to load the CSSI file
!>     call soflud(rjud_1950, dafr, sd, outr) ! then call as before
!>
!> ## sd array layout (this implementation)
!>
!>  - `sd(1..8)` = Kp for 8 3-hour periods (periods 1..8 = 0–3 h, 3–6 h, …, 21–24 h)
!>  - `sd(6)`    = 0.0  (time-reference constant; see rdymos note below)
!>  - `sd(7)`    = daily average Ap index (for rsmods Ap→Kp conversion)
!>  - `sd(9)`    = F10.7 daily observed solar flux
!>  - `sd(11)`   = F10.7 81-day centred average
!>
!> ## Compatibility notes
!>
!> **rsmods path** — uses `sd(7)` (Ap), `sd(9)` (F10.7), `sd(11)` (F10.7a).
!> All three are populated correctly; rsmods works perfectly.
!>
!> **rdymos path** — has two issues that preclude direct use:
!>
!>  1. *Uninitialized variable bug*: `rdymos` contains `rjfl = rjfl - 1.`
!>     but `rjfl` is never initialised from `Rjud`. The line should read
!>     `rjfl = Rjud - 1.`. Without this fix, `rdymos` produces garbage dates.
!>
!>  2. *Kp period conflict*: `rdymos` selects Kp via
!>     `nd = int((Dafr/3600 - sd(6) + 12 - 6.696) / 3)` then `sf(3) = sd(nd)`.
!>     With `sd(6) = 0`, `nd` maps UT 0–3 h → 1, 3–6 h → 2, …, 21–24 h → 8.
!>     However:
!>     - `nd = 6` at UT 15–18 h → `sd(6) = 0.0` (wrong; Kp treated as 0)
!>     - `nd = 7` at UT 18–21 h → `sd(7) = Ap` (wrong unit)
!>     For all other UT hours (18 of 24) the correct Kp period is returned.
!>
!> **Recommended alternative for rdymos**: use `rdymos_cssi` (below), which
!> calls `dyjrmo` directly with F10.7 and Kp obtained from `prepare_flux_data`.
!> This avoids both issues and matches the timing logic of the new Fortran model.


   !---------------------------------------------------------------------------
   !> Initialise the soflud wrapper by loading a CSSI space weather file.
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
!> Drop-in replacement for the missing SOFLUD library subroutine.
!>
!> Retrieves solar flux and geomagnetic data from the CSSI space weather
!> file loaded by `soflud_init`.

   subroutine soflud(rjud_1950, dafr, sd, outr)
      real(8), intent(in)  :: rjud_1950 !! Modified Julian Date referenced to 1950.0
                                        !! (JD − 2433282.5).  Standard MJD = rjud_1950 + 33282.
      real(8), intent(in)  :: dafr   !! Time of day in seconds (not used for data lookup;
                                        !! retained for interface compatibility).
      real(8), intent(out) :: sd(15) !! 15-element output array (see module header).
      real(8), intent(out) :: outr !! Status: 0.0 = success, non-zero = error code.

      real(dp) :: mjd
      type(flux_data_type) :: flux_data
      logical :: status

      sd   = 0.0d0
      outr = 0.0d0

      if (.not. sw_initialized) then
         write(*,'(A)') 'ERROR (soflud): not initialised — call soflud_init first.'
         outr = 1.0d0
         return
      end if

      ! Convert MJD-1950 to standard MJD
      mjd = real(rjud_1950, dp) + MJD_1950_OFFSET

      call sw_global%get_flux_data(mjd, flux_data, status)
      if (.not. status) then
         outr = 1.0d0
         write(*,'(A)') 'ERROR (soflud): space weather data lookup failed.'
         return
      end if

      ! Kp for 8 3-hour periods (sd(6) and sd(7) are overridden below)
      sd(1:8) = real(flux_data%kp(1:8), 8)

      ! sd(6) = 0.0 — time-reference constant used by rdymos nd formula.
      ! This overwrites the period-6 (15–18 h) Kp slot; see module header.
      sd(6)  = 0.0d0

      ! sd(7) = daily average Ap — used by rsmods for Ap→Kp conversion.
      ! This overwrites the period-7 (18–21 h) Kp slot; see module header.
      sd(7)  = real(flux_data%ap_avg, 8)

      ! F10.7 solar flux
      sd(9)  = real(flux_data%f107_obs, 8)

      ! F10.7 81-day centred average
      sd(11) = real(flux_data%f107a_obs_ctr, 8)

   end subroutine soflud

   !---------------------------------------------------------------------------
   !> Drop-in replacement for `rdymos` in `roberts.f90`, using the new
   !> `space_weather_module` instead of the unavailable SOFLUD library.
   !>
   !> Fixes two bugs in the original `rdymos`:
   !>  1. Uninitialised `rjfl` (should be `rjfl = Rjud - 1.`)
   !>  2. Unreliable Kp period selection via `sd(nd)` (see module header)
   !>
   !> F10.7 and Kp are obtained through `prepare_flux_data`, applying the
   !> same timing as the new model (6.7 h Kp lag, previous-day F10.7,
   !> detected-day F10.7a).  The actual computation is delegated to `rsdamo`
   !> (in roberts.f90), which calls `dyjrmo` and reorders the output densities.
   !>
   !> This subroutine has the **same interface as `rdymos`** and can be used
   !> as a direct replacement in any code that calls `rdymos`.
   !>
   !> @param[in]  Sa     Sa(1)=RA (rad), Sa(2)=geocentric lat (rad), Sa(3)=alt (m)
   !> @param[in]  Su     Su(1)=sun RA (rad), Su(2)=sun dec (rad)
   !> @param[in]  Rjud   Modified Julian Date referred to 1950.0 (JD − 2433282.5)
   !> @param[in]  Dafr   Time of day UT in seconds (0–86400)
   !> @param[in]  Gsti   Greenwich Sidereal Time in radians (not used; compatibility only)
   !> @param[out] Te     Te(1)=T∞ (K), Te(2)=local T (K)
   !> @param[out] Ad     log₁₀ number densities (m⁻³):
   !>                    Ad(1)=He, Ad(2)=O₂, Ad(3)=N₂, Ad(4)=Ar, Ad(5)=O, Ad(6)=H
   !> @param[out] Wmol   Mean molecular weight (kg/kgmol)
   !> @param[out] Rhod   Mass density (kg/m³)
   !> @param[out] status 0 = success, non-zero = space weather lookup error

   subroutine rdymos_cssi(Sa, Su, Rjud, Dafr, Gsti, Te, Ad, Wmol, Rhod, status)

      real(8), intent(in)  :: Sa(3), Su(2)
      real(8), intent(in)  :: Rjud, Dafr, Gsti
      real(8), intent(out) :: Te(2), Ad(6), Wmol, Rhod
      integer, intent(out) :: status

      type(flux_data_type) :: flux_data
      real(dp) :: utc_mjd, kp, f107, f107a
      real(8)  :: sf(3)
      logical  :: sw_status

      if (.not. sw_initialized) then
         write(*,'(A)') 'ERROR (rdymos_cssi): not initialised — call soflud_init first.'
         status = 1
         return
      end if

      ! Convert from MJD-1950 + fractional day to standard MJD
      utc_mjd = real(Rjud, dp) + MJD_1950_OFFSET + real(Dafr, dp) / 86400.0_dp

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
      sf(1) = real(f107,  8)   ! F10.7 daily (1.71-day lag)
      sf(2) = real(f107a, 8)   ! F10.7 81-day average
      sf(3) = real(kp,    8)   ! Kp (6.7-h lag)

      ! Delegate computation to rsdamo, which calls dyjrmo and reorders output
      call rsdamo(Sa, Su, sf, Rjud, Dafr, Gsti, Te, Ad, Wmol, Rhod)

      status = 0

   end subroutine rdymos_cssi

end module inpe_roberts_module