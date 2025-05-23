#include "SHELFICE_OPTIONS.h"

CBOP
      SUBROUTINE SHELFICE_INIT_VARIA( myThid )
C     *============================================================*
C     | SUBROUTINE SHELFICE_INIT_VARIA
C     | o Routine to initialize SHELFICE variables.
C     *============================================================*
C     | Initialize SHELFICE parameters and variables.
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SHELFICE.h"
#ifdef ALLOW_COST
# include "SHELFICE_COST.h"
#endif /* ALLOW_COST */
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
# include "SURFACE.h"
#endif

C     === Routine arguments ===
C     myThid -  Number of this instance of SHELFICE_INIT_VARIA
      INTEGER myThid
CEndOfInterface

#ifdef ALLOW_SHELFICE
C     === Local variables ===
C     i,j,bi,bj - Loop counters
      INTEGER i, j, bi, bj
#ifdef ALLOW_SHELFICE_GROUNDED_ICE
      INTEGER ks
#endif

CEOP

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
           shelficeForcingT      (i,j,bi,bj) = 0. _d 0
           shelficeForcingS      (i,j,bi,bj) = 0. _d 0
           shelficeHeatFlux      (i,j,bi,bj) = 0. _d 0
           shelficeFreshWaterFlux(i,j,bi,bj) = 0. _d 0
           shiTransCoeffT        (i,j,bi,bj) = 0. _d 0
           shiTransCoeffS        (i,j,bi,bj) = 0. _d 0
           shelficeMass          (i,j,bi,bj) = 0. _d 0
           shiCDragFld           (i,j,bi,bj) = 0. _d 0
           shiDragQuadFld        (i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#ifdef ALLOW_AUTODIFF
C- Note: To simplify dependency & avoid recomputations, when compiling
C        pkg/autodiff, we always call SHELFICE_INIT_VARIA to initialise control
C        variables (as done above) without condition on useShelfIce.
C        Therefore, in this case, the "If useShelfIce" is added back here:
      IF ( useShelfIce ) THEN
#endif

      IF ( .NOT. SHELFICEuseGammaFrict ) THEN
       IF ( SHELFICETransCoeffTFile .NE. ' ' ) THEN
        CALL READ_FLD_XY_RL( SHELFICETransCoeffTFile, ' ',
     &                       shiTransCoeffT, 0, myThid )
        _EXCH_XY_RL( shiTransCoeffT, myThid )
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            shiTransCoeffS(i,j,bi,bj) = SHELFICEsaltToHeatRatio *
     &                                  shiTransCoeffT(i,j,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ELSE
C     set coefficients to constant values
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            shiTransCoeffT(i,j,bi,bj) = SHELFICEheatTransCoeff
            shiTransCoeffS(i,j,bi,bj) = SHELFICEsaltTransCoeff
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
      ENDIF
C     Update the spatially varying version of the drag coefficient
C     shiCDragFld used in shelfice_thermodynamics
C     shiCDragQuadFld used in shelfice_u/v_drag_coeff
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
          shiCDragFld(i,j,bi,bj) = shiCDrag
          shiDragQuadFld(i,j,bi,bj) = SHELFICEDragQuadratic
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy, sNy+OLy
         DO i = 1-OLx, sNx+OLx
           shelficeMass(i,j,bi,bj) = shelficeMassInit(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( SHELFICEMassStepping .AND.
     &     ( nIter0.NE.0 .OR. startTime.NE.baseTime
     &                   .OR. pickupSuff.NE.' ') ) THEN
        CALL SHELFICE_READ_PICKUP( 1, nIter0, myThid )
      ENDIF

#ifdef ALLOW_SHELFICE_REMESHING
      IF ( SHI_update_kTopC ) THEN
C     set/update "kTopC" to agree with changing ice-shelf mass
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
             IF ( kSurfC(i,j,bi,bj).LE.Nr .AND.
     &            shelficeMass(i,j,bi,bj).GT.zeroRL ) THEN
               kTopC(i,j,bi,bj) = kSurfC(i,j,bi,bj)
             ELSE
               kTopC(i,j,bi,bj) = 0
             ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDIF
#endif /* ALLOW_SHELFICE_REMESHING */

#ifdef ALLOW_COST
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           cMeanSHIforT (i,j,bi,bj) = 0. _d 0
           cMeanSHIforS (i,j,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
       ENDDO
      ENDDO
c--
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
         objf_shelfice(bi,bj) = 0. _d 0
         objf_shifwflx(bi,bj) = 0. _d 0
         num_shifwflx(bi,bj)  = 0. _d 0
       ENDDO
      ENDDO
#endif /* ALLOW_COST */

#ifdef ALLOW_SHELFICE_GROUNDED_ICE

C-- Compute the mimimum water column thickness under 
C-- grounded ice
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          IF(maskInC(i,j,bi,bj).eq.1) THEN
           R_MWCT(i,j,bi,bj) =
     &      4.0*hfacmin*drF(kLowC(i,j,bi,bj))
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C-- Compute the mimimum value of r_surf (used for computing hFac_surfC)
C --This is meant to override an earlier def'n (CHECK THIS)
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1,sNy
         DO i=1,sNx
          ks = kSurfC(i,j,bi,bj)
           IF (ks.le.Nr) THEN
            Rmin_surf(i,j,bi,bj) =
!     &       R_low(i,j,bi,bj) + hFacInf*drF(ks)
     &       MAX(rF(ks+1),R_low(i,j,bi,bj)) + hFacInf*drF(ks)
           ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      _EXCH_XY_RL (Rmin_surf, myThid)

#endif
#ifdef ALLOW_AUTODIFF
C-    end if useShelfIce block
      ENDIF
#endif

#endif /* ALLOW_SHELFICE */

      RETURN
      END
