! Copyright 2012-2017, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> This module contains definitions of provenance data about e.g. how the
!> program was built, passed to Puffin from CMake


module PuffProvenance

use paratype

!  character(200), parameter  :: timeStamp =  gitBranch, puffVersion, fortCompiler, &
!                              fortVersion, fortFlags, buildHost, hostType

  character(200), parameter  :: timeStamp = "2026-07-17 13:11", &
                                gitBranch = 'master : f4263ae9361db940efbf721b865529c9d476fbfc', &
  puffVersion = '1.9.0', &
                    fortCompiler = 'gfortran', &
                    fortVersion = '8.5.0', &
                    fortFlags = " -pipe  -fopenmp", &
                    buildHost = 'sdfiana005', &
                    hostType = 'Linux-4.18.0-372.32.1.el8_6.x86_64'

end module PuffProvenance
