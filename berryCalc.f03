      program berryCalc
!
!     This program computes the geometric phase and number of enclosed conical intersections
!     of different electronic structure methods.
!
!     L. M. Thompson, 2022
!
      use MQC_Gaussian
      use omp_lib
      use iso_fortran_env, only: int32, int64, real64
!
!****x* Main/WH_PUHF
!*    NAME
!*      Geometric phase calculator. 
!*
!*    SYNOPSIS
!*      Computes the geometric phase and number of enclosed conical intersections of different 
!*      electronic structure methods.
!
      implicit none
      character(len=:),allocatable::command,fileName,help_path
      type(mqc_gaussian_unformatted_matrix_file)::fileInfo
      character(len=256),dimension(:),allocatable::fileList
      type(mqc_molecule_data),dimension(:),allocatable::molecule_list
      integer(kind=int64)::iOut=output_unit,iIn=input_unit,iPrint=1,iUnit,flag,stat_num,nullSize,&
        step,i,j,numFile=0,detLen=0
      type(mqc_scf_integral),dimension(:),allocatable::overlap_list
      type(mqc_scf_integral),dimension(:,:),allocatable::mo_list
      type(mqc_scf_integral)::overlap_matrix,shalf1,shalf2
      type(mqc_scalar)::NIJ,pNIJ,prodNIJ,Beta,NCI,nalpha,nbeta,nbasis,NIJint
      character(len=256)::method='slater'
      real(kind=real64),parameter::zero=0.0d0,half=0.5d0,one=1.0d0,zero_thresh=1.0d-8
      type(mqc_vector),dimension(:),allocatable::weights
      type(mqc_matrix)::NIJmat
!
!*    USAGE
!*      berryCalc [-f <matrix_file>] [--print-level <print_level>] [--method <method>] 
!*        [--help] 
!*
!*    OPTIONS
!* 
!
!     Print program information.
!
      write(IOut,'(*(A))') NEW_LINE('a'),' ',repeat('*',73),NEW_LINE('a'), &
          '  Geometric Phase Calculator',NEW_LINE('a'), &
          ' ',repeat('*',73),NEW_LINE('a'), &
          NEW_LINE('a'),repeat(' ',30),'Version 22.06.1',NEW_LINE('a'),NEW_LINE('a'),&
          ' L. M. Thompson, Louisville KY, 2022.',NEW_LINE('a')
!
!     Parse input options.
!
!*   1. Input/output
!*
      j = 1
      do i=1,command_argument_count()
        if(i.ne.j) cycle
        call mqc_get_command_argument(i,command)
        if(command.eq.'-f') then
!
!*      -f matrix_file                   Input file giving wavefunction at each point on closed loop. 
!*                                       The first line contains the number of matrix files in the 
!*                                       input, and then on each line is a separate matrix file.
!*
          call mqc_get_command_argument(i+1,fileName)
          j = i+2
        elseif(command.eq.'--print-level') then
!
!*      --print-level print_level        Verbosity of output. Default print level is 1. Options
!*                                       0-4.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I1)') iPrint
          j = i + 2
        elseif(command.eq.'--determinant-length') then
!
!*      --determinant-length number      Determinants are stored on the input matrix file in records 
!*                                       titled MO GRID POINT # (preceded by either ALPHA or BETA), 
!*                                       where # is a number up to the specified determinant length. 
!*                                       A vector with determinant weights should also be stored under 
!*                                       the record GRID POINT WEIGHTS, where the vector is of the same 
!*                                       length as the determinant length. Default is zero which just 
!*                                       reads from the standard MO coefficient record.
!*
          call mqc_get_command_argument(i+1,command)
          read(command,'(I10)') detLen
          j = i + 2
        elseIf(command.eq.'--help') then
!
!*      --help                           Output help documentation to terminal.
!*
          if(command_argument_count().gt.1) call mqc_error_I('Help output requested with multiple arguments',6, &
            'command_argument_count()',command_argument_count())
          call mqc_get_command_argument(0,help_path)
          help_path = 'less ' // trim(help_path(1:scan(help_path,'/',.true.))) // 'doc/TDCIS.txt'
          call execute_command_line(help_path,exitstat=flag)
          if(flag.ne.0) call mqc_error('Help output command failed')
          stop
        else
          call mqc_error_A('Unrecognised input flag',6,'command',command)
        endIf
      endDo
!
!     Parse input file and extract required data from matrix files.
!
      if(.not.allocated(fileName)) call mqc_error('No input file provided',iOut)
      open(newunit=iUnit,file=fileName,status='old',iostat=stat_num)
      if(stat_num/=0) call mqc_error_a('Error opening file',iOut,'fileName',fileName)
      read(unit=iUnit,fmt='(i20)',iostat=stat_num) numFile
      if(stat_num/=0) call mqc_error('Error reading file number',iOut)
      allocate(fileList(numFile))
      do i = 1, numFile
        read(unit=iUnit,fmt='(A)',iostat=stat_num) fileList(i)
        if((stat_num<0).and.(i<=numFile)) call mqc_error('File EOF reached early',iOut)
      endDo
      close(unit=iUnit)
      allocate(overlap_list(numFile))
      if(detLen.gt.0) then
        allocate(mo_list(numFile,detLen),weights(numFile))
      else
        allocate(mo_list(numfile,1),weights(numFile))
      endIf
      if(iPrint.ge.3) allocate(molecule_list(numfile))
      do i = 1, numFile
        if(iPrint.ge.3) then
          call fileInfo%getMolData(molecule_list(i),filename=fileList(i))
          call molecule_list(i)%print(iOut)
          call mqc_print(mqc_get_nuclear_repulsion(molecule_list(i)),iOut,'Vnn from geometry '//trim(num2char(i)))
        endIf
        call fileInfo%getESTObj('overlap',est_integral=overlap_list(i),filename=fileList(i))
        if(iPrint.ge.4) call overlap_list(i)%print(iOut,'Overlap matrix from geometry '//trim(num2char(i)))
        if(detLen.eq.0) then
          call fileInfo%getESTObj('mo coefficients',est_integral=mo_list(i,1),filename=fileList(i))
          if(iPrint.ge.4) call mo_list(i,1)%print(iOut,'MO coefficients from geometry '//trim(num2char(i)))
          weights(i) = [1.0]
        else
          do j = 1, detLen
            call fileInfo%getESTObj('MO GRID POINT '//trim(num2char(j)),est_integral=mo_list(i,j),filename=fileList(i))
            if(iPrint.ge.4) call mo_list(i,j)%print(iOut,'MO coefficients from geometry '//trim(num2char(i))//&
              ', determinant '//trim(num2char(j)))
          endDo
          call fileInfo%getArray('GRID POINT WEIGHTS',vectorOut=weights(i),filename=fileList(i))
        endIf
        if(iPrint.ge.4) call weights(i)%print(iOut,'Determinant weights from geometry '//trim(num2char(i)))
      endDo
      call fileInfo%load(fileList(1))
      nBasis = fileInfo%getVal('nBasis')
      nAlpha = fileInfo%getVal('nAlpha')
      nBeta = fileInfo%getVal('nBeta')

      prodNIJ = cmplx(One,Zero)
      do step= 1,numFile
      
        shalf1 = overlap_list(step)
        shalf2 = overlap_list(mod(step,numfile)+1)
        call shalf1%power(half)
        call shalf2%power(half)
        overlap_matrix = matmul(shalf1,shalf2)

        call NIJmat%init(max(detLen,1),max(detlen,1))
        do i = 1, max(detLen,1)
          do j = 1, max(detLen,1)
            call get_nij(NIJint,pnIJ,nullSize,mo_list(step,i),mo_list(mod(step,numfile)+1,j),overlap_matrix,nBasis,nAlpha,nBeta)
            call NIJmat%put(NIJInt,i,j)
          endDo
        endDo
        if(step.ne.numfile) then 
          call phase_rotate(iout,iprint,weights(step+1),weights(step),NIJmat,NIJ)
        else
          call compute_NIJ(max(detLen,1),weights(step),weights(mod(step,numfile)+1),NIJmat,NIJ)
        endIf

        if(iPrint.ge.2) call NIJ%print(iOut,'NIJ at step '//trim(num2char(step)),FormatStr='F20.15')
        prodNIJ = prodNIJ*NIJ

      end do

      if(iPrint.ge.2) call prodNIJ%print(iOut,'Pancharatnam product of overlaps')
      Beta = (-1)*cmplx(zero,one)*log(prodNIJ)
      call Beta%print(iOut,'Geometric phase',FormatStr='F20.15')
      NCI = Beta/pi
      call NCI%print(iOut,'Number of conical intersections enclosed in loop',FormatStr='F20.15')
      
      contains

!
!     PROCEDURE get_nij
!     
!     get_nij is a subroutine that returns the overlap and psuedooverlap of two
!     (potentially nonorthogonal) Slater determinants, as well as the dimension
!     of the overlap null space, the overlap and psuedo-overlap matrix elements. The
!     sign of nullSize gives the multiple of matrix elements accounting for antisymmetry
!     due to permutation of orbitals in the SVD.
!
      subroutine get_nij(nIJ,pnIJ,nullSize,mo_I,mo_J,overlap,nBasis,nAlpha,nBeta)

      implicit none

!     input/output variables
      type(mqc_scalar),intent(inOut)::nIJ,pnIJ
      integer(kind=int64)::nullSize
      type(mqc_scf_integral),intent(in)::mo_I,mo_J,overlap
      type(mqc_scalar),intent(in)::nBasis,nAlpha,nBeta

!     subroutine variables
      type(mqc_matrix)::mo_I_occ,mo_J_occ,mIJ,uMat,vMat
      type(mqc_vector)::sigmaMat
      logical::orthflag
      integer(kind=int64)::i

!
      mo_I_occ = mqc_integral_output_block(mo_I%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      mo_J_occ = mqc_integral_output_block(mo_J%orbitals('occupied',[int(nAlpha)],[int(nBeta)]),'full') 
      
      mIJ = matmul(matmul(dagger(mo_I_occ),overlap%getBlock('full')),mo_J_occ)
!      call mIJ%print(6,'LMTLMT mij')
      nIJ = mIJ%det()

      orthflag = .false.
      if((nIJ%abs()).lt.zero_thresh) then
        call mIJ%svd(EVals=sigmaMat,EUVecs=uMat,EVVecs=vMat)
        if(minval(abs(sigmaMat)).lt.zero_thresh) orthflag = .true.
      endIf

      if(orthflag) then
        pnIJ = 1.0
        nullSize = 0
        do i = 1,size(sigmaMat)
          if(sigmaMat%at(i).gt.zero_thresh) then
            pnIJ = pnIJ*sigmaMat%at(i)
          else
            nullSize = nullSize + 1
          endIf
        endDo
        nullSize = sign(1.0,real(uMat%det()))*sign(1.0,real(vMat%det()))*nullSize
      else
        pnIJ = nIJ
        nullSize = 0
      endIf
!
      end subroutine get_nij
!    
!
      subroutine phase_rotate(iOut,iPrint,weights_rotate,weights_fixed,NIJ_matrix,NIJ)
!
      implicit none
      type(mqc_vector),intent(inOut)::weights_rotate
      type(mqc_vector),intent(in)::weights_fixed
      type(mqc_matrix),intent(in)::NIJ_matrix
      integer,intent(in)::iOut,iPrint
      type(mqc_scalar),intent(out)::NIJ
      type(mqc_vector)::thetas,norms
      type(mqc_scalar)::threshzero,trial_angle,norm,step,err
      integer::k,max_iters=150,iters
      real(kind=real64),parameter::zero=0.0d0
!
      threshzero = 1.0e-12
!
      thetas = [0.0,pi/2,2*pi]
      call norms%init(3,(100.0,0.0))
      do k=1,3
        trial_angle = thetas%at(k)
        if(iPrint.ge.3) call trial_angle%print(iout,'theta')
        call compute_NIJ(size(NIJ_matrix,1),weights_fixed,exp(cmplx(zero,float(trial_angle)))*weights_rotate,NIJ_matrix,norm)
        call norms%put(norm,k)
      endDo
      step = ((3-sqrt(5.0))/2)
      iters = 1
      do while (iters.le.max_iters)
        if(iPrint.ge.3) then
          call thetas%print(iout,'Theta list')
          call norms%print(iout,'Norm list')
        endIf
        if(abs(thetas%at(2)-thetas%at(1)).ge. &
          abs(thetas%at(2)-thetas%at(3))) then
          trial_angle = thetas%at(2) + &
            step*(thetas%at(2)-thetas%at(1))
          if(iPrint.ge.3) call trial_angle%print(iout,'new theta')
          call compute_NIJ(size(NIJ_matrix,1),weights_fixed,exp(cmplx(zero,float(trial_angle)))*weights_rotate,NIJ_matrix,err)
          if(iPrint.ge.3) call err%print(iOut,'new norm')
          if(err.gt.norms%at(2)) then
            call thetas%put(thetas%at(2),3)
            call thetas%put(trial_angle,2)
            call norms%put(norms%at(2),3)
            call norms%put(err,2)
          else
            call thetas%put(trial_angle,1)
            call norms%put(err,1)
          endIf
        else
          trial_angle = thetas%at(2) + &
            step*(thetas%at(3)-thetas%at(2))
          if(iPrint.ge.3) call trial_angle%print(iout,'new theta')
          call compute_NIJ(size(NIJ_matrix,1),weights_fixed,exp(cmplx(zero,float(trial_angle)))*weights_rotate,NIJ_matrix,err)
          if(iPrint.ge.3) call err%print(iOut,'new norm')
          if(err.gt.norms%at(2)) then
            call thetas%put(thetas%at(2),1)
            call thetas%put(trial_angle,2)
            call norms%put(norms%at(2),1)
            call norms%put(err,2)
          else
            call thetas%put(trial_angle,3)
            call norms%put(err,3)
          endIf
        endIf
        iters = iters + 1
        if((abs(norms%at(1)-norms%at(2)).lt.threshzero).and. &
          (abs(norms%at(2)-norms%at(3)).lt.threshzero)) then
          if(iPrint.ge.3) write(iOut,'(A)') ' Phase rotation converged'
          if(iPrint.ge.3) call trial_angle%print(iOut,'Greatest overlap achieved with phase angle')
          if(iPrint.ge.3) call weights_rotate%print(iout,'Initial weights_rotate')
          weights_rotate = exp(cmplx(zero,float(trial_angle)))*weights_rotate
          call compute_NIJ(size(NIJ_matrix,1),weights_fixed,weights_rotate,NIJ_matrix,NIJ)
          if(iPrint.ge.3) call weights_rotate%print(iout,'Rotated weights_rotate')
          if(iPrint.ge.3) call weights_fixed%print(iout,'Compare to weights_fixed')
          exit
        endIf
        if(iters.eq.max_iters) then
          if(iPrint.ge.1) write(iOut,'(A)') ' Phase rotation did not converge'
          if(iPrint.ge.1) write(iOut,'(A)') ' Continuing with input phases'
          exit
        endIf
      endDo

      end subroutine phase_rotate 
!
!
      subroutine compute_NIJ(detLen,bra_weights,ket_weights,NIJ_matrix,NIJ)
!
      implicit none
      type(mqc_vector),intent(in)::bra_weights,ket_weights
      type(mqc_matrix),intent(in)::NIJ_matrix
      integer,intent(in)::detLen
      type(mqc_scalar),intent(out)::NIJ
      integer::i,j
      real(kind=real64),parameter::zero=0.0d0
!
      NIJ = zero
      do i = 1, max(detLen,1)
        do j = 1, max(detLen,1)
          NIJ = NIJ + conjg(bra_weights%at(i))*ket_weights%at(j)*NIJ_matrix%at(i,j)
        endDo
      endDo
!
      end subroutine compute_NIJ
!
!
!*    NOTES
!*      Compilation of this program requires the MQC library (https://github.com/MQCPack/mqcPack)
!*      and the gauopen utility (http://gaussian.com/g16/gauopen.zip) and compilation with the
!*      f08 standard.
!*
!*      Compilation tested using: gfortran 9.2.0
!*
!*      Note that subroutine Wr_LCBuf needs modifying in gauopen/qcmatrix.F as follows:
!*        line 58:       LenBX = (LenBuf/(2*abs(NR)))*abs(NR)
!*        line 60:       Call Wr_CBuf(IU,NTot*abs(NR),LenBX,X)
!*
!*      Documentation generated with robodoc. To update documentation edit robodoc.rc to
!*      determine documentation output type and then run robodoc at the command line in the
!*      main directory.
!*
!*    AUTHORS
!*      Lee M. Thompson, University of Louisville, lee.thompson.1@lousiville.edu
!*
!*    COPYRIGHT
!*      (c) 2022 by Lee M. Thompson distributed under terms of the MIT license.
!*
!****
!
      end program berryCalc
