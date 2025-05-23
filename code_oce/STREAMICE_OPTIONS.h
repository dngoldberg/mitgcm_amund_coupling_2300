C $Header: /u/gcmpack/MITgcm_contrib/verification_other/shelfice_remeshing/code/STREAMICE_OPTIONS.h,v 1.3 2016/04/02 23:13:10 jmc Exp $
C $Name:  $

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C CPP options file for MYPACKAGE
C
C Use this file for selecting options within package "streamice"

#ifndef STREAMICE_OPTIONS_H
#define STREAMICE_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_STREAMICE

#include "CPP_OPTIONS.h"

C Place CPP define/undef flag here

#define STREAMICE_CONSTRUCT_MATRIX
#define STREAMICE_HYBRID_STRESS
#define USE_ALT_RLOW
#define STREAMICE_GEOM_FILE_SETUP
#define STREAMICE_SMOOTH_FLOATATION2
!#define ALLOW_PETSC
!#define STREAMICE_FLOWLINE_BUTTRESS
#define STREAMICE_COULOMB_SLIDING

#endif /* ALLOW_MYPACKAGE */
#endif /* MYPACKAGE_OPTIONS_H */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
