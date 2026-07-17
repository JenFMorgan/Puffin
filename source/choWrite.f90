! ###############################################
! Copyright 2012-2018, University of Strathclyde
! Authors: Lawrence T. Campbell
! License: BSD-3-Clause
! ###############################################

!> @author
!> Lawrence Campbell,
!> University of Strathclyde,
!> Glasgow, UK
!> @brief
!> This module contains top-level subroutines to write data in SDDS or HDF5
!> format.

module dummyf

USE lattice
USE RK4int
use hdf5_puff
use ParaField
use cwrites

implicit none

contains


subroutine writeIM(sZ, sZl, &
                   iStep, iCstep, iL, iWriteNthSteps, &
                   iIntWriteNthSteps, nSteps, qOK)


! Subroutine to write data, making the necessary
! preparations before writing due to the layout
! of data when integrating in the undulator module.
!
! Lawrence Campbell
! University of Strathclyde
! Jan 2015

  implicit none

  real(kind=wp), intent(inout) :: sZ, sZl
  integer(kind=ip), intent(in) :: iStep, iWriteNthSteps, iIntWriteNthSteps, nSteps
  integer(kind=ip), intent(in) :: iCstep, iL
  integer(kind=ip) :: nslices
  logical, intent(inout) :: qOK

  integer error

  logical :: qOKL, qWriteInt, qWriteFull

  qOK = .false.


  call int_or_full(istep, iCstep, iIntWriteNthSteps, iWriteNthSteps, &
                   qWriteInt, qWriteFull, qOK)


  call wr_cho(sZ, sZl, &
              iStep, iCstep, iL, iWriteNthSteps, &
              iIntWriteNthSteps, nSteps, qWriteInt, qWriteFull, qOK)

!              Set error flag and exit

  qOK = .TRUE.

  goto 2000

1000  call Error_log('Error in writeIM',tErrorLog_G)

2000 continue



end subroutine writeIM






subroutine wr_cho(sZ, sZl, &
                  iStep, iCstep, iL, iWriteNthSteps, &
                  iIntWriteNthSteps, nSteps, qWriteInt, qWriteFull, qOK)


! Subroutine to write data, choosing either sdds or hdf5 (or both!)
!
! Lawrence Campbell
! University of Strathclyde
! Jan 2017

  implicit none

  real(kind=wp), intent(inout) :: sZ, sZl
  integer(kind=ip), intent(in) :: iStep, iWriteNthSteps, iIntWriteNthSteps, nSteps
  integer(kind=ip), intent(in) :: iCstep, iL
  logical, intent(in) :: qWriteInt, qWriteFull
  logical, intent(inout) :: qOK

  integer(kind=ip) :: nslices
  integer error

  logical :: qOKL

  if (qhdf5_G) then

    if (fieldMesh == iTemporal) then
      nslices=ceiling( (sLengthOfElmZ2_G*NZ2_G)/(4*pi*srho_g))
    else
      nslices=ceiling( (sLengthOfElmZ2_G * real((NZ2_G-1_ip),kind=wp) )/(4*pi*srho_g)) ! + 30_ip
    end if

    call wr_h5(sZ, szl, tArrayA, tArrayE, tArrayZ, iL, &
               iIntWriteNthSteps, iWriteNthSteps, qSeparateStepFiles_G, &
               qWriteFull, &
               qWriteInt, nslices, qOK)

  end if

end subroutine wr_cho



!          ##################################################################



  subroutine int_or_full(istep, iCsteps, iIntWr, iWr, &
                         qWriteInt, qWriteFull, qOK)

    implicit none

!   Figure out whether to write integrated data or
!   full particle dump


    integer(kind=ip), intent(in) :: istep, iCsteps
    integer(kind=ip), intent(in) :: iIntWr, iWr
    logical, intent(inout) :: qWriteInt, qWriteFull, qOK

    integer(kind=ip) :: iw

    logical ::  qOKL

    qOK = .false.

    qWriteInt = .false.
    qWriteFull = .false.
    
    

    if (qInitWrLat_G) then

      if ((mod(iStep,iIntWr)==0) .or. (iStep == nSteps) .or. (iStep == 0) ) then

        qWriteInt = .true.

      end if


      if ((mod(iStep,iWr)==0) .or. (iStep == nSteps) .or. (iStep == 0) ) then

        qWriteFull = .true.

      end if

    else

      if ((mod(iCsteps,iIntWr)==0) .or. (iCsteps == nSteps) .or. (iCsteps == 0) ) then

        qWriteInt = .true.

      end if


      if ((mod(iCsteps,iWr)==0) .or. (iCsteps == 0) ) then

        qWriteFull = .true.

      end if

    end if

  if (qWrArray_G) then

    do iw = 1, size(wrarray)

      if (wrarray(iw) == iCsteps) then
        qWriteFull = .true.
        qWriteInt = .true.
      end if

    end do

  end if

  end subroutine int_or_full


!########################################################################


function qWriteq(iStep, iCsteps, iWriteNthSteps, iIntWriteNthSteps, nSteps)



  implicit none

  logical :: qWriteq
  integer(kind=ip) :: iStep, iCsteps, iWriteNthSteps, iIntWriteNthSteps, nSteps
  integer(kind=ip) :: iw


  if (qInitWrLat_G) then

    if ((mod(iStep,iIntWriteNthSteps)==0) .or. (iStep == nSteps) &
                 .or. (mod(iStep,iWriteNthSteps)==0)) then

      qWriteq = .true.

    else

      qWriteq = .false.

    end if

  else

    if ((mod(iCsteps,iIntWriteNthSteps)==0) .or. (iCsteps == nSteps) &
                 .or. (mod(iCsteps,iWriteNthSteps)==0)) then

      qWriteq = .true.

    else

      qWriteq = .false.

    end if

  end if


  if (qWrArray_G) then

    do iw = 1, size(wrarray)

      if (wrarray(iw) == iCsteps) qWriteq = .true.

    end do

  end if


end function qWriteq

!> @author
!> Jenny Morgan,
!> SLAC,
!> Menlo Park, USA
!> @brief
!> function to check if energy change should be applied.

function qEshiftq(iStep, iCsteps, iwakefieldNthSteps, ispacechargeNthSteps, nSteps)
  implicit none
  
  logical :: qEshiftq
  integer(kind=ip) :: iStep, iCsteps, iwakefieldNthSteps, ispacechargeNthSteps, nSteps
  integer(kind=ip) :: iw

  if ((qspacecharge_G .or. qwake_G))  then
     if (qInitWrLat_G) then

        if ((mod(iStep,ispacechargeNthSteps)==0) .or. (iStep == nSteps) &
                     .or. (mod(iStep,iwakefieldNthSteps)==0)) then

          qEshiftq = .true.

       else

          qEshiftq = .false.

        end if

      else

        if ((mod(iCsteps,ispacechargeNthSteps)==0) .or. (iCsteps == nSteps) &
                     .or. (mod(iCsteps,iwakefieldNthSteps)==0)) then

         qEshiftq = .true.

       else

         qEshiftq = .false.

       end if

     end if
    
    else

      qEshiftq = .false.


    end if 



!  if (qWrArray_G) then
!
!    do iw = 1, size(wrarray)
!
!      if (wrarray(iw) == iCsteps) qWriteq = .true.
!
!    end do
! not sure what this does shall ignore for now. 
  
!  end if
end function qEshiftq



subroutine EnergyshiftIM(sZ, sZl, &
  iStep, iCstep, iL, iwakefieldNthSteps, & 
  ispacechargeNthSteps, &
  nSteps, qOK)


  ! Top subroutine to apply space charge or wakefeilds. Checking which to apply. 

   implicit none
 
   real(kind=wp), intent(inout) :: sZ, sZl
   integer(kind=ip), intent(in) :: iStep, iwakefieldNthSteps, ispacechargeNthSteps, nSteps
   integer(kind=ip), intent(in) :: iCstep, iL
   integer(kind=ip) :: nslices
   logical, intent(inout) :: qOK

   integer error

   logical :: qOKL, qApplyWake, qApplyspaceCharge

   qOK = .false.


   call wake_or_spacecharge(istep, iCstep, iwakefieldNthSteps, ispacechargeNthSteps, &
   qApplyWake, qApplyspaceCharge, qOK)


   call Energy_cho(sZ, sZl, &
   iStep, iCstep, iL, iwakefieldNthSteps, &
   ispacechargeNthSteps, nSteps, qApplyWake, qApplyspaceCharge, qOK)

   qOK = .TRUE.

   goto 2000

1000  call Error_log('Error in  EnergyshiftIM',tErrorLog_G)

2000 continue



end subroutine EnergyshiftIM

subroutine wake_or_spacecharge(istep, iCsteps, iApwake, iApspace, &
  qApplyWake, qApplyspaceCharge, qOK)

   implicit none

   !   Figure out whether to apply wake or space charge


   integer(kind=ip), intent(in) :: istep, iCsteps
   integer(kind=ip), intent(in) :: iApwake, iApspace
   logical, intent(inout) :: qApplyWake, qApplyspaceCharge, qOK

   logical ::  qOKL
 
   qOK = .false.

   qApplyWake = .false.
   qApplyspaceCharge = .false.
 
   if (qInitWrLat_G) then

   if (qwake_G .and. ((mod(iStep,iApwake)==0) .or. (iStep == nSteps) .or. (iStep == 0)) ) then

     qApplyWake = .true.

   end if


   if ( qspacecharge_G .and. ((mod(iStep,iApspace)==0) .or. (iStep == nSteps) .or. (iStep == 0)) ) then

      qApplyspaceCharge = .true.

   end if

   else

   if (qwake_G .and. ((mod(iCsteps,iApwake)==0) .or. (iCsteps == nSteps) .or. (iCsteps == 0)) ) then

      qApplyWake = .true.

   end if


   if (qspacecharge_G .and. ((mod(iCsteps,iApspace)==0) .or. (iCsteps == 0)) ) then

      qApplyspaceCharge = .true.

   end if

   end if

!if (qWrArray_G) then

!do iw = 1, size(wrarray)

!if (wrarray(iw) == iCsteps) then
!qWriteFull = .true.
!qWriteInt = .true.
!end if

!end do

!end if

end subroutine wake_or_spacecharge

subroutine Energy_cho(sZ, sZl, &
  iStep, iCstep, iL, iwakefieldNthSteps, &
  ispacechargeNthSteps, nSteps, qApplyWake, qApplyspaceCharge, qOK)


! Subroutine to apply wakefeilds space charge or both


implicit none

real(kind=wp), intent(inout) :: sZ, sZl
integer(kind=ip), intent(in) :: iStep, iwakefieldNthSteps, ispacechargeNthSteps, nSteps
integer(kind=ip), intent(in) :: iCstep, iL
logical, intent(in) :: qApplyWake, qApplyspaceCharge
logical, intent(inout) :: qOK

integer(kind=ip) :: nslices
integer error

logical :: qOKL


if (fieldMesh == iTemporal) then
  nslices=ceiling( (sLengthOfElmZ2_G*NZ2_G)/(4*pi*srho_g))
else
  nslices=ceiling( (sLengthOfElmZ2_G * real((NZ2_G-1_ip),kind=wp) )/(4*pi*srho_g)) ! + 30_ip
end if

call ApplyEnergyshift(sZ, szl, iL, &
                      iwakefieldNthSteps, ispacechargeNthSteps, &
                      qApplyWake, &
                      qApplyspaceCharge, nslices, qOK)



end subroutine Energy_cho


subroutine Energy_cho_drift(sZ, &
  iL, iwakefieldNthSteps, &
  ispacechargeNthSteps, qApplyWake, qApplyspaceCharge, qOK, del_dr_z)


! Subroutine to apply wakefeilds space charge or both


implicit none

real(kind=wp), intent(inout) :: sZ
integer(kind=ip), intent(in) :: iwakefieldNthSteps, ispacechargeNthSteps
integer(kind=ip), intent(in) :: iL
logical, intent(in) :: qApplyWake, qApplyspaceCharge
logical, intent(inout) :: qOK
real(kind=wp) :: szl
integer(kind=ip) :: nslices
integer error
real(kind=wp), intent(in) :: del_dr_z
logical :: qOKL

szl = del_dr_z

if (fieldMesh == iTemporal) then
  nslices=ceiling( (sLengthOfElmZ2_G*NZ2_G)/(4*pi*srho_g))
else
  nslices=ceiling( (sLengthOfElmZ2_G * real((NZ2_G-1_ip),kind=wp) )/(4*pi*srho_g)) ! + 30_ip
end if

zUndType_G = 'Drift'

call ApplyEnergyshift(sZ, szl, iL, &
                      iwakefieldNthSteps, ispacechargeNthSteps, &
                      qApplyWake, &
                      qApplyspaceCharge, nslices, qOK)



end subroutine Energy_cho_drift



end module dummyf
