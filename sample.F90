#ifdef MPI
# define _MPI_
# define CLUSTER
# undef MPI
#endif

module mpiconf
  !> Arguments: idtask, nproc, wid
  use mpi
  integer  :: idtask
  integer  :: nproc
  logical  :: wid
  integer  :: MPIerror

  private :: MPIerror
  public  :: idtask, nproc, wid
  public  :: mpiconf_init
  save
contains
  subroutine mpiconf_init()
#ifdef _MPI_
    call mpi_init(ierr)
    call MPI_Comm_Rank( MPI_Comm_World, idtask, MPIerror )
    call MPI_Comm_Size( MPI_Comm_World, nproc, MPIerror )
    if (MPIerror /= 0) error stop " Error in initializing MPI "
    wid = (idtask==0)
#endif
  end subroutine mpiconf_init
end module mpiconf




module mpi_write
  use mpiconf,      only: wid
  use iso_fortran_env, only: OUTPUT_UNIT
  implicit none
  integer, parameter   :: sp = kind(1.0)
  integer, parameter   :: dp = kind(1.0d0)
  integer              :: ounit = OUTPUT_UNIT  ! for the time being
  character(len=100)   :: fmt_single  = '(A, T40, ":: ", T50, G0)'
  character(len=100)   :: fmt_double  = '(A, T40, ":: ", T50, G0)'
  character(len=100)   :: fmt_int     = '(A, T40, ":: ", T50, I0)'
  character(len=100)   :: fmt_string  = '(A, T40, ":: ", T50, A)'
  character(len=100)   :: fmt_logical = '(A, T40, ":: ", T50, L1)'
  character(len=100)   :: fmt_list_int    = '(A, T40, ":: ", T50, 100(I0,:,", "))'
  character(len=100)   :: fmt_list_single = '(A, T40, ":: ", T50, 100(G0,:,", "))'
  character(len=100)   :: fmt_list_double = '(A, T40, ":: ", T50, 100(G0,:,", "))'
  character(len=100)   :: fmt_message = '(A)'
  private
  public :: mpiwrite, wid
  save ! this is very important

  interface mpiwrite
    module procedure mpi_write_integer
    module procedure mpi_write_real
    module procedure mpi_write_realdp
    module procedure mpi_write_string
    module procedure mpi_write_logical
    module procedure mpi_write_list_integer
    module procedure mpi_write_list_single
    module procedure mpi_write_list_double
    module procedure mpi_write_message
  end interface mpiwrite

contains
  integer function mpi_write_integer(message, value )
    character(len=*), intent(in)  :: message
    integer, intent(in)           :: value
    integer                       :: iostat
    if (wid) then
      write(ounit,fmt=fmt_int,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing integer to the output file "
    mpi_write_integer = iostat
    return
  end function mpi_write_integer

  integer function mpi_write_real(message, value )
    character(len=*), intent(in)  :: message
    real(sp), intent(in)          :: value
    integer                       :: iostat
    if (wid) then
      write(ounit,fmt_single,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing float to the output file "
    mpi_write_real = iostat
    return
  end function mpi_write_real

  integer function mpi_write_realdp(message, value )
    character(len=*), intent(in)  :: message
    real(dp), intent(in)          :: value
    integer                       :: iostat
    if (wid) then
      write(ounit,fmt_double,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing double precision to the output file "
    mpi_write_realdp = iostat
    return
  end function mpi_write_realdp

  integer function mpi_write_string(message, value )
    character(len=*), intent(in)            :: message
    character(len=*), intent(in)            :: value
    integer                                 :: iostat
    if (wid) then
      write(ounit,fmt_string,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing string to the output file "
    mpi_write_string = iostat
    return
  end function mpi_write_string

  integer function mpi_write_logical(message, value )
    character(len=*), intent(in)  :: message
    logical, intent(in)           :: value
    integer                       :: iostat
    if (wid) then
      write(ounit,fmt_logical,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing boolean to the output file "
    mpi_write_logical = iostat
    return
  end function mpi_write_logical

  integer function mpi_write_list_integer(message, value)
    character(len=*), intent(in)    :: message
    integer, intent(in),dimension(:):: value
    integer                         :: iostat
    if (wid) then
      write(ounit,fmt=fmt_list_int,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing list of integers to the output file "
    mpi_write_list_integer = iostat
    return
  end function mpi_write_list_integer

  integer function mpi_write_list_single(message, value)
    character(len=*), intent(in)    :: message
    real, intent(in),dimension(:)   :: value
    integer                         :: iostat
    if (wid) then
      write(ounit,fmt=fmt_list_single,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing list of floats to the output file "
    mpi_write_list_single = iostat
    return
  end function mpi_write_list_single

  integer function mpi_write_list_double(message, value)
    character(len=*), intent(in)        :: message
    real(dp), intent(in),dimension(:)   :: value
    integer                             :: iostat
    if (wid) then
      write(ounit,fmt=fmt_list_double,iostat=iostat) message, value
    endif
    if (iostat .ne. 0) stop "Error in writing list of doubles to the output file "
    mpi_write_list_double = iostat
    return
  end function mpi_write_list_double

  integer function mpi_write_message(message)
    character(len=*), intent(in)  :: message
    integer                       :: iostat
    if (wid) then
      write(ounit,fmt=fmt_message,iostat=iostat) message
    endif
    if (iostat .ne. 0) stop "Error in writing a message to the output file "
    mpi_write_message = iostat
    return
  end function mpi_write_message


end module mpi_write




PROGRAM SAMPLE
  USE fdf
  USE prec
  use mpiconf, only: idtask, nproc, wid
  use mpiconf, only: mpiconf_init
  use mpi
  use mpi_write
  implicit none
!--------------------------------------------------------------- Local Variables
  integer, parameter         :: maxa = 100
  logical                    :: doit, debug, mpiflag, exist, opened, closed
  character(len=100)         :: fname, file_data
  character(len=2)           :: axis, status
  integer(sp)                :: i, j, ia, na, external_entry, iunit, iostat
  integer(sp)                :: isa(maxa), MPIerror
  real(sp)                   :: wmix
  real(dp)                   :: cutoff, phonon_energy, factor
  real(dp)                   :: xa(3, maxa), weights(10)
  real(dp)                   :: listdp(maxa),  var1, var2, var3
  real(dp)                   :: listr(maxa)
  real                       :: calibration
  type(block_fdf)            :: bfdf
  type(parsed_line), pointer :: pline

  logical, save              :: initialized=.false.
  integer                    :: err, natoms, ierr
  real(dp), allocatable      :: cent(:,:)
  character(2), allocatable  :: symbol(:)



!------------------------------------------------------------------------- BEGIN

  call mpiconf_init()


! Initialize
  call fdf_init('sample.fdf', 'sample.out')


! Handle/Use fdf structure
  ! if (fdf_defined('new-style')) then
  !   err = mpiwrite("[strings] just a string without value")
  ! endif

  err = mpiwrite(" node id         ", idtask)
  err = mpiwrite(" total processes ", nproc)

  na = fdf_integer('NumberOfAtoms', 0)
!  write(*,'(a,i0,a,i0)') " node id ", idtask, " of ", nproc
  err = mpiwrite("[integer] Number of atoms (via calls) ", na)
  if (wid) write(*,*) "[integer] Number of atoms ", na

  fname = fdf_string('NameOfFile', 'whatever')
  err = mpiwrite("[strings] Name of the file", fname)

  cutoff = fdf_physical('MeshCutoff', 8.d0, 'Ry')
  err = mpiwrite("[floats dp] energy cutoff", cutoff)

  calibration = fdf_single('calibration', 0.0)
  err = mpiwrite("[floats sp] calibration", calibration)

  debug = fdf_boolean('Debug', .TRUE.)
  err = mpiwrite("[boolean] debug flag", debug)

  if (fdf_block('AtomicCoordinatesAndAtomicSpecies', bfdf)) then
    ia = 1
    do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
      do i= 1, 3
        xa(i,ia) = fdf_breals(pline, i)
      enddo
      isa(ia) = fdf_bintegers(pline, 1)
      ia = ia + 1
    enddo
  endif

  if (wid) then
    write(6,*) 'Atomic Coordinates from the MPI process id :', idtask
    do ia= 1, na
      write(6,'(3F10.6,I5)') (xa(i,ia),i=1,3), isa(ia)
    enddo
  endif




!list of integers
  if ( fdf_islist('MyList') ) then
     na = -1
     call fdf_list('MyList',na,isa)
     call fdf_list('MyList',na,isa)
     err = mpiwrite("[list][integers] list of integers   ", isa(1:na))
  else
     write(*,*)'MyList was not recognized'
  end if


  if (fdf_block('molecule', bfdf)) then
    !   External file reading
!        write(6,*) 'Beginning of molecular coordinates block  '
        ia = 1

        do while((fdf_bline(bfdf, pline)))
!         get the integer from the first line
          if ((pline%id(1) .eq. "i") .and. (pline%ntokens .eq. 1)) then        ! check if it is the only integer present in a line
            natoms = fdf_bintegers(pline, 1)
            if (wid) write(*,*) "Molecule block :: number of atoms  ", natoms
          endif

          if (.not. allocated(cent)) allocate(cent(3,natoms))
          if (.not. allocated(symbol)) allocate(symbol(natoms))

          if (pline%ntokens == 4) then
            symbol(ia) = fdf_bnames(pline, 1)
            do i= 1, 3
              cent(i,ia) = fdf_bvalues(pline, i)
            enddo
            ia = ia + 1
          endif
        enddo


  endif

  if (wid) then
    write(6,*) 'Coordinates from Molecule block: '
    do ia= 1, natoms
      write(6,'(A4,3F10.6)') symbol(ia), (cent(i,ia),i=1,3)
    enddo

    write(6,'(A)')
    write(6,*) '------------------------------------------------------'
  endif


  if (wid) then
    inquire(file="data.txt", exist=exist, opened=opened)
    if (exist) then
        open (newunit=iunit,file="data.txt", iostat=iostat, action='read' )
        if (iostat .ne. 0) stop "Problem in opening the data file"
    else
        error stop " data file "// trim(file_data) // " does not exist."
    endif


    read(iunit,*) var1
    read(iunit,*) var2
    read(iunit,*) var3
  endif

  call MPI_Bcast(var1,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(var2,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)
  call MPI_Bcast(var3,1,MPI_Double_Precision,0,MPI_Comm_World,MPIerror)

  if (idtask == 2 .or. idtask == 3 ) then
    write(6,*) 'checking the broadcasted message var1 '
    write(6,'(3(G0, 4x),a,i0)') var1, var2, var3, " from idtask = ", idtask
    write(6,'(A)')
    write(6,*) '------------------------------------------------------'
  endif

! list of single precision floats NOT IMPLEMENTED YET
  if ( fdf_islreal('MyListReal') .and. fdf_islist('MyListReal') &
      .and. (.not. fdf_islinteger('MyListReal')) ) then
    na = -1
    call fdf_list('MyListReal',na,listr)
    if ( na < 2 ) stop 1
    call fdf_list('MyListReal',na,listr)
    err = mpiwrite("[list][floats] list of single floats   ", listr(1:na))
  else
    write(*,*)'MyListR was not recognized'
    stop 1
  end if



  if ( fdf_islreal('weights') .and. fdf_islist('weights') &
  .and. (.not. fdf_islinteger('weights')) ) then
i = -1
call fdf_list('weights',i,weights)
!write(ounit,'(tr1,a,i0,a)') 'weights has ',i,' entries'
if ( i < 2 ) stop 1
call fdf_list('weights',i,weights)
err = mpiwrite( 'Weights : ', weights(1:i))
  endif




!   phonon_energy = fdf_physical('phonon-energy', 0.01d0, 'eV')
!   write(6,*) 'Phonon Energy:', phonon_energy

!   i = fdf_integer('SomeInt', 34)
!   write(6,*) '#elems:', i

!   wmix = fdf_single('WmixValue', 0.55)
!   write(6,*) 'WmixValue:', wmix

!   factor = fdf_double('FactorValue', 1.d-10)
!   write(6,*) 'Factor:', factor


!   doit = fdf_boolean('DoIt', .FALSE.)
!   write(6,*) 'Doit:', doit

!   doit = fdf_defined('AtomicCoordinatesAndAtomicSpecies')
!   write(6,*) 'AtomCoordsBlockDefined:', doit


!   if (fdf_block('AtomicInfo', bfdf)) then
!     ia = 1
!     do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
!       do i= 1, 3
!         xa(i,ia) = fdf_breals(pline, i)
!       enddo
!       ia = ia + 1
!     enddo
!   endif

!   write(6,*) 'AtomicInfo:'
!   do ia= 1, na
!     write(6,'(3F10.6)') (xa(i,ia),i=1,3)
!   enddo

!   if (fdf_block('Other-Block', bfdf)) then

! !   Forward reading
!     ia = 1
!     do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
!       symbol(ia) = fdf_bnames(pline, 1)
!       do i= 1, na
!         xa(i,ia) = fdf_breals(pline, i)
!       enddo
!       ia = ia + 1
!     enddo

!     write(6,*) 'Other-Block (Forward):'
!     do ia= 1, na
!       write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
!     enddo

! !   Backward reading
!     ia = 1
!     do while((fdf_bbackspace(bfdf, pline)) .and. (ia .le. na))
!       symbol(ia) = fdf_bnames(pline, 1)
!       do i= 1, na
!         xa(i,ia) = fdf_breals(pline, i)
!       enddo
!       ia = ia + 1
!     enddo

!     write(6,*) 'Other-Block (Backward):'
!     do ia= 1, na
!       write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
!     enddo

! !   Forward reading
!     ia = 1
!     do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
!       symbol(ia) = fdf_bnames(pline, 1)
!       do i= 1, na
!         xa(i,ia) = fdf_breals(pline, i)
!       enddo
!       ia = ia + 1
!     enddo

!     write(6,*) 'Other-Block (Forward):'
!     do ia= 1, na
!       write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
!     enddo

! !   Forward reading with rewind
!     call fdf_brewind(bfdf)
!     ia = 1
!     do while((fdf_bline(bfdf, pline)) .and. (ia .le. na))
!       symbol(ia) = fdf_bnames(pline, 1)
!       do i= 1, na
!         xa(i,ia) = fdf_breals(pline, i)
!       enddo
!       ia = ia + 1
!     enddo

!     write(6,*) 'Other-Block (Forward-with-rewind):'
!     do ia= 1, na
!       write(6,'(A4,3F10.6)') symbol(ia), (xa(i,ia),i=1,3)
!     enddo
!   endif

!   if ( fdf_block('ListBlock',bfdf) ) then
!      i = 0
!      do while ( fdf_bline(bfdf,pline) )
!         i = i + 1
!         na = fdf_bnlists(pline)
!         write(*,'(2(a,i0),a)') 'Listblock line: ',i,' has ',na,' lists'
!         do ia = 1 , na
!            j = -1
!            call fdf_blists(pline,ia,j,isa)
!            write(*,'(tr5,2(a,i0),a)') 'list ',ia,' has ',j,' entries'
!            call fdf_blists(pline,ia,j,isa)
!            write(*,'(tr5,a,1000(tr1,i0))') 'list: ',isa(1:j)
!         end do
!      end do
!   end if



!   if ( fdf_islist('externalentry') ) then
!      write(*,*) 'externalentry is a list'
!   else
!      write(*,*) 'externalentry is not a list'
!   end if

!   external_entry = fdf_integer('externalentry', 60)
!   write(6,*) 'ExternalEntry:', external_entry

!   axis   = fdf_string('AxisXY', 'Cartesian')
!   status = fdf_string('StatusXY', 'Enabled')
!   write(6,*) 'Axis: ', TRIM(axis), ' | ', TRIM(status)

! Shutdown and deallocates fdf structure
  call fdf_shutdown()
  call mpi_finalize(ierr)
!----------------------------------------------------------------------------END
END PROGRAM SAMPLE
