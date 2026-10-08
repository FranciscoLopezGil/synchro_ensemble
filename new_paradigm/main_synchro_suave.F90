!=================================================================
! GHOST new_paradigm + SYNCHRO_SUAVE (continuous nudging)
!
! This file is the upstream src/main.F90 of the GHOST new_paradigm
! branch (commit afaf9ea, 2026-09-26) plus the minimal additions for
! the synchronized-ensemble experiments. Inside PROGRAM MAIN every
! addition is enclosed between the markers
!     ! >>> SYNCHRO   ...   ! <<< SYNCHRO
! so that "diff src/main.F90 main_synchro_suave.F90" shows only the
! synchronization logic. No other GHOST file is modified.
!
! Experiment: one reference run plus num_realizations members, all
! integrated with the same solver, stepper and forcing. Instead of
! overwriting the shells k <= kndg after each step, every member gets
! the extra (linear relaxation, "nudging") force
!     F_nudge(k) = nudge_lambda*[u_ref(k) - u_member(k)],  k <= kndg
! added to its copy of the forcing and integrated by the time stepper.
! F_nudge is evaluated with the states at the beginning of the step
! and kept fixed during the substeps (like the forcing), so on the
! nudged shells the difference relaxes as (1 - nudge_lambda*dt) per
! step: nudge_lambda*dt must be < 1 (checked at start-up). The shell of
! a mode is int(|k|/Dkk+.501), the criterion GHOST uses for its spectra
! (spectrumc), so kndg=4 means the shells 0,1,2,3,4.
!
! Input : parameter.inp     reference run + namelist &ensembleparams
!         parameterNNN.inp  one file per member, NNN = 001, 002, ...
!                           (&pdes, &HD with fidir/fodir/todir,
!                           &initial_conditions, &stepper)
! Output: per member, in its todir: the usual GHOST txt files plus
!         difference.txt     time, <|du|^2>, <|curl du|^2>, <du.f>
!         dspectrum.NNNN.txt k, E_du(k)   (same normalization as
!                            kspectrum: sum_k E(k).Dkk = <|du|^2>/2)
!         with du = u_member - u_reference (f is the physical forcing,
!         the members' balance.txt does not include the nudging power).
!         Binaries in its fodir. The reference also appends one line
!         to forcing_check.txt (in its todir) every time the run starts.
!
! &ensembleparams (in parameter.inp):
!   num_realizations   : number of members (1..999)
!   kndg               : last nudged shell (k <= kndg)
!   nudge_lambda       : relaxation rate lambda (1/time), > 0
!   ensemble_seed_base : seed of the members' initial small scales;
!                        member ir uses ensemble_seed_base+10007*ir
!   seed_forcing       : seed used ONLY to build the forcing, so that
!                        a restart rebuilds exactly the same forcing
!   ensemble_ic        : 'scramble' (default) member = reference with
!                        random phases (one per wave vector) in k>kndg
!                        'iclist' member = its own &initial_conditions
!                        (legacy protocol), then k<=kndg replaced once
!   ensemble_init      : 'auto' (default) new members if stat=0, read
!                        the member files if stat>0; 'new' always
!                        creates the members from the reference (e.g.
!                        to start an ensemble from a spun-up reference
!                        with stat>0); 'read' always reads them
!=================================================================

      MODULE synchro_suave_mod

      USE fprecision
      USE gstate_mod
      USE equationbase_mod
      USE gstepperbase_mod
      USE ic_factory
      IMPLICIT NONE

      INTEGER, PARAMETER :: SYNC_I8 = SELECTED_INT_KIND(18)
      PRIVATE :: hash32, mulmod32, mode_key, mode_phase

! One synchronized realization: its own solver (with its I/O dirs),
! states, IC chain and stepper, and its forcing (the forcing of the
! reference plus the nudging term). The workspace pool is that of the
! reference run.
      TYPE :: GSyncMember
         TYPE(GStateComp),    ALLOCATABLE :: field(:),field_nxt(:)
         TYPE(GStateComp),    ALLOCATABLE :: force_nudged(:)
         CLASS(EquationBase), ALLOCATABLE :: pde
         CLASS(GStepperBase), ALLOCATABLE :: stepper
         TYPE(icChain),       ALLOCATABLE :: iclist(:)
         CHARACTER(len=256)               :: infile
      END TYPE GSyncMember

      CONTAINS

!=================================================================
! Stops all MPI tasks with a message (printed by task 0)
!=================================================================
      SUBROUTINE sync_abort(msg)
      USE commtypes
      USE mpivars
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: msg
      IF (myrank.eq.0) THEN
         WRITE(*,'(A)') 'SYNCHRO ERROR: '//TRIM(msg)
         FLUSH(6)
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD,1)
      END SUBROUTINE sync_abort

!=================================================================
! Reads &ensembleparams from infile (task 0) and broadcasts it.
! Defaults are those of the previous Synchro mains.
!=================================================================
      SUBROUTINE read_ensemble_params(infile,variant,nreal,kcut,lambda,  &
                 sbase,sforc,icmode,initmode,newmembers)
      USE commtypes
      USE mpivars
      USE status, ONLY: stat
      IMPLICIT NONE
      CHARACTER(len=*) , INTENT(IN)  :: infile,variant
      INTEGER          , INTENT(OUT) :: nreal,kcut,sbase,sforc
      REAL(KIND=GP)    , INTENT(OUT) :: lambda
      CHARACTER(len=16), INTENT(OUT) :: icmode,initmode
      LOGICAL          , INTENT(OUT) :: newmembers
      INTEGER            :: num_realizations,kndg,ensemble_seed_base
      INTEGER            :: seed_forcing,ios
      REAL(KIND=GP)      :: nudge_lambda
      CHARACTER(len=16)  :: ensemble_ic,ensemble_init
      CHARACTER(len=256) :: msg
      NAMELIST / ensembleparams / num_realizations,kndg,nudge_lambda,   &
                                  ensemble_seed_base,seed_forcing,      &
                                  ensemble_ic,ensemble_init

      num_realizations   = 1
      kndg               = 1
      nudge_lambda       = 1.0_GP
      ensemble_seed_base = 50000
      seed_forcing       = 70001
      ensemble_ic        = 'scramble'
      ensemble_init      = 'auto'
      ios = 0
      msg = ''
      IF (myrank.eq.0) THEN
         OPEN(1,file=infile,status='old',form='formatted',iostat=ios,  &
              iomsg=msg)
         IF (ios.eq.0) THEN
            READ(1,NML=ensembleparams,iostat=ios,iomsg=msg)
            CLOSE(1)
         ENDIF
      ENDIF
      CALL MPI_BCAST(ios,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF (ios.ne.0) CALL sync_abort('cannot read &ensembleparams in '// &
                         TRIM(infile)//': '//TRIM(msg))
      CALL MPI_BCAST(num_realizations  ,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(kndg              ,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nudge_lambda      ,1 ,GC_REAL      ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ensemble_seed_base,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed_forcing      ,1 ,MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ensemble_ic       ,16,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ensemble_init     ,16,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      IF ((num_realizations.lt.1).or.(num_realizations.gt.999))         &
         CALL sync_abort('num_realizations must be in [1,999]')
      IF (kndg.lt.0) CALL sync_abort('kndg must be >= 0')
      IF ((TRIM(ensemble_ic).ne.'scramble').and.                         &
          (TRIM(ensemble_ic).ne.'iclist'))                               &
         CALL sync_abort("ensemble_ic must be 'scramble' or 'iclist'")
      SELECT CASE (TRIM(ensemble_init))
      CASE ('auto')
         newmembers = (stat.eq.0)
      CASE ('new')
         newmembers = .TRUE.
      CASE ('read')
         newmembers = .FALSE.
      CASE DEFAULT
         CALL sync_abort("ensemble_init must be 'auto', 'new' or 'read'")
      END SELECT
      IF ((.not.newmembers).and.(stat.eq.0))                             &
         CALL sync_abort("ensemble_init='read' needs stat > 0")

      nreal    = num_realizations
      kcut     = kndg
      lambda   = nudge_lambda
      sbase    = ensemble_seed_base
      sforc    = seed_forcing
      icmode   = ensemble_ic
      initmode = ensemble_init
      IF (myrank.eq.0) THEN
         WRITE(*,'(A)')      ' Synchro variant    : '//TRIM(variant)
         WRITE(*,'(A,I0)')   '   num_realizations : ',nreal
         WRITE(*,'(A,I0)')   '   kndg (shells<=)  : ',kcut
         WRITE(*,'(A,I0)')   '   ensemble_seed_base: ',sbase
         WRITE(*,'(A,I0)')   '   seed_forcing     : ',sforc
         WRITE(*,'(A)')      '   ensemble_ic      : '//TRIM(icmode)
         IF (newmembers) THEN
            WRITE(*,'(A)')   '   ensemble_init    : '//TRIM(initmode)// &
                             ' -> members created from the reference'
         ELSE
            WRITE(*,'(A)')   '   ensemble_init    : '//TRIM(initmode)// &
                             ' -> members read from their files'
         ENDIF
      ENDIF
      END SUBROUTINE read_ensemble_params

!=================================================================
! A member must be an HD run with the same physical parameters and
! number of fields as the reference (diagnostics use components 1-3)
!=================================================================
      SUBROUTINE check_member(m,ref)
      USE hd_mod, ONLY: HDSolver
      IMPLICIT NONE
      TYPE(GSyncMember)  , INTENT(IN) :: m
      CLASS(EquationBase), INTENT(IN) :: ref

      IF (m%pde%state_size().ne.ref%state_size())                       &
         CALL sync_abort(TRIM(m%infile)//': number of fields differs '//&
                         'from parameter.inp')
      SELECT TYPE (r => ref)
      TYPE IS (HDSolver)
         SELECT TYPE (p => m%pde)
         TYPE IS (HDSolver)
            IF ((p%traits_%nu.ne.r%traits_%nu).or.                       &
                (p%traits_%dorot.neqv.r%traits_%dorot).or.               &
                ANY(p%traits_%omega.ne.r%traits_%omega).or.              &
                (p%numpassive_.ne.r%numpassive_))                        &
               CALL sync_abort(TRIM(m%infile)//': &HD (nu, dorot, '//    &
                    'omega, npassive) differs from parameter.inp')
         CLASS DEFAULT
            CALL sync_abort(TRIM(m%infile)//': the member must use '//   &
                            'solver = "HD"')
         END SELECT
      CLASS DEFAULT
         CALL sync_abort('the Synchro mains support only solver = "HD"')
      END SELECT
      END SUBROUTINE check_member

!=================================================================
! To continue a member (stat>0) its IC chain must read its files
!=================================================================
      SUBROUTINE require_read_v(m)
      USE ic_velocity, ONLY: icRead_v
      IMPLICIT NONE
      TYPE(GSyncMember), INTENT(IN) :: m
      LOGICAL :: found
      INTEGER :: n
      found = .FALSE.
      DO n = 1,SIZE(m%iclist)
         SELECT TYPE (ic => m%iclist(n)%ic)
         TYPE IS (icRead_v)
            found = .TRUE.
         END SELECT
      END DO
      IF (.not.found) CALL sync_abort(TRIM(m%infile)//': iclist must '// &
          'contain read_v to continue the member (or use ensemble_init='//&
          "'new')")
      END SUBROUTINE require_read_v

!=================================================================
! Sets the host copies of a state to zero. Needed before random ICs
! and random forcings: they only write the modes of their band, and
! GState_alloc does not initialize memory.
!=================================================================
      SUBROUTINE zero_state_host(state)
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(INOUT) :: state(:)
      INTEGER :: ic
      DO ic = 1,SIZE(state)
         state(ic)%ccomp = CMPLX(0.0_GP,0.0_GP,KIND=GP)
      END DO
      END SUBROUTINE zero_state_host

!=================================================================
! dst = src on the HOST copies (GState_copy copies only the device
! copies in offload builds; this is used during initialization)
!=================================================================
      SUBROUTINE copy_state_host(dst,src)
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(INOUT) :: dst(:)
      TYPE(GStateComp), INTENT(IN)    :: src(:)
      INTEGER :: ic
      IF (SIZE(dst).ne.SIZE(src)) CALL sync_abort('copy_state_host: '// &
                                       'inconsistent state sizes')
      DO ic = 1,SIZE(src)
         dst(ic)%ccomp = src(ic)%ccomp
      END DO
      END SUBROUTINE copy_state_host

!=================================================================
! dst(k) = src(k) for every mode in the shells k <= kcut (all state
! components). Runs on the device copies while gdev_active is set
! (offload builds), on the host otherwise, following the kernel
! convention of GHOST (components addressed through pointers).
!=================================================================
      SUBROUTINE replace_scales(src,dst,kcut)
      USE mpivars
      USE grid
      USE kes
      USE boxsize
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(INOUT), TARGET :: src(:)
      TYPE(GStateComp), INTENT(INOUT), TARGET :: dst(:)
      INTEGER         , INTENT(IN)            :: kcut
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: ps,pd
      REAL(KIND=GP) :: rdk
      INTEGER       :: ic,i,j,k

      rdk = Dkk
      DO ic = 1,SIZE(src)
         ps => src(ic)%ccomp
         pd => dst(ic)%ccomp
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (int(sqrt(kk2(k,j,i))/rdk+.501).le.kcut) THEN
                     pd(k,j,i) = ps(k,j,i)
                  ENDIF
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE replace_scales

!=================================================================
! frc(k) = frc(k) + lambda*[ref(k) - mem(k)] for every mode in the
! shells k <= kcut (all state components). Same kernel convention as
! replace_scales (device copies while gdev_active is set).
!=================================================================
      SUBROUTINE add_nudge_forcing(ref,mem,kcut,lambda,frc)
      USE mpivars
      USE grid
      USE kes
      USE boxsize
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(INOUT), TARGET :: ref(:),mem(:)
      TYPE(GStateComp), INTENT(INOUT), TARGET :: frc(:)
      INTEGER         , INTENT(IN)            :: kcut
      REAL(KIND=GP)   , INTENT(IN)            :: lambda
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: pr,pm,pf
      REAL(KIND=GP) :: rdk,lam
      INTEGER       :: ic,i,j,k

      rdk = Dkk
      lam = lambda
      DO ic = 1,SIZE(ref)
         pr => ref(ic)%ccomp
         pm => mem(ic)%ccomp
         pf => frc(ic)%ccomp
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  IF (int(sqrt(kk2(k,j,i))/rdk+.501).le.kcut) THEN
                     pf(k,j,i) = pf(k,j,i)+lam*(pr(k,j,i)-pm(k,j,i))
                  ENDIF
               END DO
            END DO
         END DO
      END DO
      END SUBROUTINE add_nudge_forcing

!=================================================================
! Multiplies every mode in the shells k > kcut by exp(i.phi(k)), with
! ONE random phase per wave vector shared by all the components: the
! modal amplitudes (spectrum, helicity) and k.u = 0 are preserved.
! Hermitian symmetry is kept in the kx=0 plane (phi(-k) = -phi(k),
! self-conjugate modes untouched). phi(k) is a hash of (iseed, global
! mode index), so the result does not depend on the number of MPI
! tasks or threads. Runs on the host copies.
!=================================================================
      SUBROUTINE scramble_phases_above(fld,kcut,iseed)
      USE mpivars
      USE grid
      USE kes
      USE boxsize
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(INOUT) :: fld(:)
      INTEGER         , INTENT(IN)    :: kcut,iseed
      COMPLEX(KIND=GP) :: cph
      REAL(KIND=GP)    :: phase
      INTEGER          :: ic,i,j,k

!$omp parallel do collapse(2) private (k,ic,phase,cph)
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               IF (int(sqrt(kk2(k,j,i))/Dkk+.501).gt.kcut) THEN
                  phase = mode_phase(i,j,k,iseed)
                  cph   = CMPLX(COS(phase),SIN(phase),KIND=GP)
                  DO ic = 1,SIZE(fld)
                     fld(ic)%ccomp(k,j,i) = fld(ic)%ccomp(k,j,i)*cph
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      END SUBROUTINE scramble_phases_above

! Global index of the mode (k,j,i) of the (nz,ny,nx/2+1) array
      PURE FUNCTION mode_key(i,j,k) RESULT(key)
      USE grid, ONLY: ny,nz
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: i,j,k
      INTEGER(KIND=SYNC_I8) :: key
      key = (INT(i-1,SYNC_I8)*INT(ny,SYNC_I8)+INT(j-1,SYNC_I8))*         &
            INT(nz,SYNC_I8)+INT(k-1,SYNC_I8)
      END FUNCTION mode_key

! Random phase in [0,2.pi) of a mode, odd under k -> -k in kx=0
      PURE FUNCTION mode_phase(i,j,k,iseed) RESULT(phase)
      USE grid, ONLY: ny,nz
      USE var , ONLY: pi
      IMPLICIT NONE
      INTEGER, INTENT(IN)   :: i,j,k,iseed
      REAL(KIND=GP)         :: phase
      INTEGER(KIND=SYNC_I8) :: key,keyp,h0,h1
      REAL(KIND=GP)         :: sgn
      key = mode_key(i,j,k)
      sgn = 1.0_GP
      IF (i.eq.1) THEN               ! (k,j) and its conjugate partner
         keyp = mode_key(1,MOD(ny-j+1,ny)+1,MOD(nz-k+1,nz)+1)
         IF (keyp.eq.key) THEN       ! self-conjugate mode: keep it
            phase = 0.0_GP
            RETURN
         ELSE IF (keyp.lt.key) THEN  ! use the phase of the partner
            key = keyp
            sgn = -1.0_GP
         ENDIF
      ENDIF
      h0 = hash32(INT(iseed,SYNC_I8))
      h1 = hash32(IEOR(h0,key))
      h1 = hash32(IEOR(h1,hash32(h0+1_SYNC_I8)))
      phase = sgn*2.0_GP*pi*REAL((REAL(h1,KIND=8)+0.5D0)/4294967296.0D0, &
                                 KIND=GP)
      END FUNCTION mode_phase

! 32-bit integer hash (lowbias32, C. Wellons), in 64-bit arithmetic
      PURE FUNCTION hash32(x) RESULT(h)
      IMPLICIT NONE
      INTEGER(KIND=SYNC_I8), INTENT(IN) :: x
      INTEGER(KIND=SYNC_I8)             :: h
      INTEGER(KIND=SYNC_I8), PARAMETER  :: M32 = 4294967295_SYNC_I8
      h = IAND(x,M32)
      h = IEOR(h,ISHFT(h,-16))
      h = mulmod32(h,2146121005_SYNC_I8)          ! 0x7feb352d
      h = IEOR(h,ISHFT(h,-15))
      h = mulmod32(h,2221713035_SYNC_I8)          ! 0x846ca68b
      h = IEOR(h,ISHFT(h,-16))
      END FUNCTION hash32

! (a*c) mod 2^32 for 0 <= a,c < 2^32 without overflowing 64 bits
      PURE FUNCTION mulmod32(a,c) RESULT(r)
      IMPLICIT NONE
      INTEGER(KIND=SYNC_I8), INTENT(IN) :: a,c
      INTEGER(KIND=SYNC_I8)             :: r
      INTEGER(KIND=SYNC_I8), PARAMETER  :: M16 = 65535_SYNC_I8
      INTEGER(KIND=SYNC_I8), PARAMETER  :: M32 = 4294967295_SYNC_I8
      r = a*IAND(c,M16)+ISHFT(IAND(a*ISHFT(c,-16),M16),16)
      r = IAND(r,M32)
      END FUNCTION mulmod32

!=================================================================
! Same as GHOST's hdcheck (pseudospec_hd) for the difference field,
! writing 'difference.txt' instead of 'balance.txt':
!   time, <|du|^2>, <|curl du|^2>, <du.f>
! (and helicity.txt / divergence.txt of du if hel=1 / chk=1). Uses
! host-only temporaries like hdcheck (host copies of the fields).
!=================================================================
      SUBROUTINE diff_check(a,b,c,d,e,f,t,dt,hel,chk,path)
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      USE pseudospec_fluid
      USE class_GWorkspace3D, ONLY: gws
      IMPLICIT NONE
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3
      DOUBLE PRECISION    :: eng,ens,pot,khe
      DOUBLE PRECISION    :: div,tmp
      REAL(KIND=GP), INTENT(IN) :: dt
      REAL(KIND=GP)       :: tmq
      INTEGER, INTENT(IN) :: hel,chk,t
      INTEGER             :: i,j,k
      CHARACTER(len=*), INTENT(IN) :: path
      LOGICAL :: bret_

      CALL gws%get_complex_htmp(c1,bret_)
      CALL gws%get_complex_htmp(c2,bret_)
      CALL gws%get_complex_htmp(c3,bret_)
      div = 0.0D0
      tmp = 0.0D0
      tmq = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
      IF (chk.eq.1) THEN
         CALL derivk3(a,c1,1)
         CALL derivk3(b,c2,2)
         CALL derivk3(c,c3,3)
         IF (ista.eq.1) THEN
!$omp parallel private (k) reduction(+:tmp)
!$omp do
            DO j = 1,ny
               DO k = 1,nz
                  tmp = tmp+abs(c1(k,j,1)+c2(k,j,1)+c3(k,j,1))**2*tmq
               END DO
            END DO
!$omp end do
!$omp do collapse(2)
            DO i = 2,iend
               DO j = 1,ny
                  DO k = 1,nz
                     tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
                  END DO
               END DO
            END DO
!$omp end do
!$omp end parallel
         ELSE
!$omp parallel do collapse(2) private (k) reduction(+:tmp)
            DO i = ista,iend
               DO j = 1,ny
                  DO k = 1,nz
                     tmp = tmp+2*abs(c1(k,j,i)+c2(k,j,i)+c3(k,j,i))**2*tmq
                  END DO
               END DO
            END DO
         ENDIF
         CALL MPI_REDUCE(tmp,div,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,       &
                         MPI_COMM_WORLD,ierr)
      ENDIF
      CALL energy(a,b,c,eng,1)
      CALL energy(a,b,c,ens,0)
      IF (hel.eq.1) CALL helicity(a,b,c,khe)
      CALL cross(a,b,c,d,e,f,pot,1)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path)//'/difference.txt',position='append')
         WRITE(1,10) (t-1)*dt,eng,ens,pot
   10    FORMAT( E13.6,E26.18,E26.18,E26.18 )
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file=trim(path)//'/dhelicity.txt',position='append')
            WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,khe
            CLOSE(1)
         ENDIF
         IF (chk.eq.1) THEN
            OPEN(1,file=trim(path)//'/ddivergence.txt',position='append')
            WRITE(1,FMT='(E13.6,E26.18)') (t-1)*dt,div
            CLOSE(1)
         ENDIF
      ENDIF
      CALL gws%free_complex_htmp(c1)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c3)
      END SUBROUTINE diff_check

!=================================================================
! Energy spectrum of the difference field, 'dspectrum.NNNN.txt',
! with the normalization of GHOST's spectrum: E = sum[E(k).Dkk]
!=================================================================
      SUBROUTINE diff_spectrum(a,b,c,path,nmb)
      USE kes
      USE grid
      USE mpivars
      USE boxsize
      USE pseudospec_hd
      IMPLICIT NONE
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      CHARACTER(len=*), INTENT(IN) :: path,nmb
      DOUBLE PRECISION, DIMENSION(nmax/2+1) :: Ek,Hk
      INTEGER :: i

      CALL spectrumc(a,b,c,1,0,Ek,Hk)
      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(path)//'/dspectrum.'//nmb//'.txt')
         DO i = 1,nmax/2+1
            WRITE(1,FMT='(E13.6,E23.15)') Dkk*i,.5_GP*Ek(i)/Dkk
         END DO
         CLOSE(1)
      ENDIF
      END SUBROUTINE diff_spectrum

!=================================================================
! Global quantities of a member: the usual GHOST output of its solver
! (balance.txt, ... in its todir) and difference.txt for du. The
! difference field is built in normal pool temporaries, so that the
! host-only pool stays available for the diagnostics (as in hdcheck).
!=================================================================
      SUBROUTINE member_global(m,ref,force,t)
      USE status, ONLY: dt
      IMPLICIT NONE
      TYPE(GSyncMember), INTENT(IN) :: m
      TYPE(GStateComp) , INTENT(IN) :: ref(:),force(:)
      INTEGER          , INTENT(IN) :: t
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: dvx,dvy,dvz
      LOGICAL :: bret

      CALL m%pde%global(m%field_nxt,force,t)
      CALL m%pde%workspace_%get_complex_tmp(dvx,bret)
      CALL m%pde%workspace_%get_complex_tmp(dvy,bret)
      CALL m%pde%workspace_%get_complex_tmp(dvz,bret)
      dvx = m%field_nxt(1)%ccomp-ref(1)%ccomp
      dvy = m%field_nxt(2)%ccomp-ref(2)%ccomp
      dvz = m%field_nxt(3)%ccomp-ref(3)%ccomp
      CALL diff_check(dvx,dvy,dvz,force(1)%ccomp,force(2)%ccomp,         &
                      force(3)%ccomp,t,dt,0,0,m%pde%todir_)
      CALL m%pde%workspace_%free_complex_tmp(dvz)
      CALL m%pde%workspace_%free_complex_tmp(dvy)
      CALL m%pde%workspace_%free_complex_tmp(dvx)
      END SUBROUTINE member_global

!=================================================================
! Spectra of a member: the usual GHOST spectra of its solver and
! dspectrum.NNNN.txt for du (NNNN = sind, as for kspectrum)
!=================================================================
      SUBROUTINE member_spectra(m,ref)
      USE status, ONLY: sind
      USE filefmt
      IMPLICIT NONE
      TYPE(GSyncMember), INTENT(IN) :: m
      TYPE(GStateComp) , INTENT(IN) :: ref(:)
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: dvx,dvy,dvz
      LOGICAL :: bret

      CALL m%pde%spectra(m%field_nxt)
      CALL m%pde%workspace_%get_complex_tmp(dvx,bret)
      CALL m%pde%workspace_%get_complex_tmp(dvy,bret)
      CALL m%pde%workspace_%get_complex_tmp(dvz,bret)
      dvx = m%field_nxt(1)%ccomp-ref(1)%ccomp
      dvy = m%field_nxt(2)%ccomp-ref(2)%ccomp
      dvz = m%field_nxt(3)%ccomp-ref(3)%ccomp
      WRITE(ext,fmtext) sind
      CALL diff_spectrum(dvx,dvy,dvz,m%pde%todir_,ext)
      CALL m%pde%workspace_%free_complex_tmp(dvz)
      CALL m%pde%workspace_%free_complex_tmp(dvy)
      CALL m%pde%workspace_%free_complex_tmp(dvx)
      END SUBROUTINE member_spectra

!=================================================================
! Checks of the forcing built at start-up (host copy), printed by
! task 0 and appended to path/forcing_check.txt:
!  - <|f|^2> and a phase-sensitive checksum: equal values in two runs
!    mean the same forcing realization (compare restarts with them);
!  - whether tasks other than 0 drew random numbers while building
!    it (seed changed): then the realization depends on the MPI
!    decomposition (random_fv draws per task, from the same seed);
!  - largest forced shell, compared with kndg.
!=================================================================
      SUBROUTINE forcing_report(force,sforc,seed_after,fstatic,kcut,path)
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE boxsize
      USE status, ONLY: stat
      IMPLICIT NONE
      TYPE(GStateComp), INTENT(IN) :: force(:)
      INTEGER         , INTENT(IN) :: sforc,seed_after,kcut
      LOGICAL         , INTENT(IN) :: fstatic
      CHARACTER(len=*), INTENT(IN) :: path
      DOUBLE PRECISION :: loc(2),glo(2),w,th,tmp,fr,fi
      INTEGER          :: ic,i,j,k,kfl,kfg,ldraw,gdraw

      tmp = 1.0D0/(DBLE(nx)*DBLE(ny)*DBLE(nz))**2
      loc = 0.0D0
      kfl = -1
      DO ic = 1,SIZE(force)
         DO i = ista,iend
            w = 2.0D0
            IF (i.eq.1) w = 1.0D0
            DO j = 1,ny
               DO k = 1,nz
                  fr = DBLE(REAL (force(ic)%ccomp(k,j,i)))
                  fi = DBLE(AIMAG(force(ic)%ccomp(k,j,i)))
                  IF ((fr.ne.0.0D0).or.(fi.ne.0.0D0)) THEN
                     th = 0.1D0*DBLE(ic+3*i+7*j+11*k)
                     loc(1) = loc(1)+w*(fr**2+fi**2)*tmp
                     loc(2) = loc(2)+w*(fr*COS(th)+fi*SIN(th))*tmp
                     kfl = MAX(kfl,int(sqrt(kk2(k,j,i))/Dkk+.501))
                  ENDIF
               END DO
            END DO
         END DO
      END DO
      ldraw = 0
      IF ((myrank.ne.0).and.(seed_after.ne.sforc)) ldraw = 1
      CALL MPI_REDUCE(loc,glo,2,MPI_DOUBLE_PRECISION,MPI_SUM,0,          &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(kfl,kfg,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(ldraw,gdraw,1,MPI_INTEGER,MPI_SUM,0,               &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         WRITE(*,'(A,I0,A,I0,A,I0,A,L1)') ' Forcing check: stat=',stat,  &
              ' nprocs=',nprocs,' seed_forcing=',sforc,' static=',fstatic
         WRITE(*,'(A,ES23.15,A,ES23.15,A,I0)') '   <|f|^2>=',glo(1),      &
              '  checksum=',glo(2),'  last forced shell=',kfg
         IF (gdraw.gt.0) THEN
            WRITE(*,'(A,I0,A)') '   WARNING: ',gdraw,' task(s) other '// &
              'than 0 drew random numbers: this forcing realization '//  &
              'depends on the number of MPI tasks. Continue this run '// &
              'only with the same nprocs.'
         ELSE
            WRITE(*,'(A)') '   All random draws (if any) on task 0: '//  &
              'the forcing does not depend on the MPI decomposition.'
         ENDIF
         IF (.not.fstatic) WRITE(*,'(A)') '   WARNING: time-dependent '// &
              'forcing update: restarts rebuild the initial forcing but'//&
              ' not the later random updates (use constant_fv).'
         IF (kfg.gt.kcut) WRITE(*,'(A)') '   NOTE: the forcing acts on '//&
              'shells k > kndg (not synchronized).'
         OPEN(1,file=trim(path)//'/forcing_check.txt',position='append')
         WRITE(1,'(I8,I8,I12,2ES25.15,I6,I6,L3)') stat,nprocs,sforc,     &
              glo(1),glo(2),kfg,gdraw,fstatic
         CLOSE(1)
      ENDIF
      END SUBROUTINE forcing_report

      END MODULE synchro_suave_mod

!=================================================================
      PROGRAM MAIN
!=================================================================
! GHOST code: Geophysical High Order Suite for Turbulence
!
! Numerically integrates several fluid dynamics equations
! in 2 and 3 dimensions with periodic boundary conditions
! and external forcing. A pseudo-spectral method is used to
! compute spatial derivatives, while different time steppings
! can be used to evolve the system in the time domain.
!
! Notation: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!
! References:
! Mininni PD, Rosenberg DL, Reddy R, Pouquet A.; P.Comp.37, 123 (2011)
! Rosenberg DL, Mininni PD, Reddy R, Pouquet A.: Atmosph.11, 178 (2020)
!=================================================================

! Modules
      USE commtypes
      USE filefmt
      USE iovar
      USE fft
      USE threads
      USE gdevice
      USE boxsize
      USE status
      USE pstatus
      USE gtimer
      USE gbench
      USE ic_factory
      USE force_factory
      USE equation_factory
      USE particle_factory
      USE icp_factory
      USE stepper_factory
! >>> SYNCHRO
      USE synchro_suave_mod
! <<< SYNCHRO
      IMPLICIT NONE

! Arrays for the field and particle states, workspace, I/O, and solver classes
      TYPE   (GStateComp), ALLOCATABLE :: field(:),field_nxt(:),force (:)
      TYPE  (GPStateComp), ALLOCATABLE :: part (:),part_nxt (:)
      TYPE   (GWorkspace), TARGET      :: workspace
      TYPE       (ioplan)              :: planio
      CLASS(EquationBase), ALLOCATABLE :: fluid
      CLASS(ParticleBase), ALLOCATABLE :: particle
      CLASS(GStepperBase), ALLOCATABLE :: stepper
      CLASS     (icChain), ALLOCATABLE :: iclist(:)
      CLASS  (forceChain), ALLOCATABLE :: forcemethod(:)
      CLASS    (icpChain), ALLOCATABLE :: icplist(:)

! Auxiliary variables
      REAL(KIND=GP) :: time
      INTEGER       :: t, num_components
      LOGICAL       :: fstatic

! >>> SYNCHRO
! Ensemble members and parameters (see the header of this file)
      TYPE(GSyncMember), ALLOCATABLE, TARGET :: ensemble(:)
      REAL(KIND=GP)     :: nudge_lambda
      INTEGER           :: num_realizations,kndg,ensemble_seed_base
      INTEGER           :: seed_forcing,seed_reference,seed_after,ir
      LOGICAL           :: newmembers
      CHARACTER(len=16) :: ensemble_ic,ensemble_init
      CHARACTER(len=8)  :: outlabel
! <<< SYNCHRO

! Initialization
! Initializes the MPI and I/O libraries
      CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

! Initializes the grid. This must be done early to have nx, ny, nz.
      CALL grid_init('parameter.inp')

! Binds this MPI task to a GPU in offload builds (no-op otherwise).
! GHOST assumes the number of MPI tasks per node is a multiple of the
! number of GPUs per node; the user must ensure this condition.
      CALL device_init(myrank)

! Initialization of fluid and particles integration parameters
      CALL status_init ('parameter.inp')
      CALL pstatus_init('parameter.inp')

! >>> SYNCHRO
! The status seed is kept for the ICs of the reference; the forcing is
! built with its own seed (seed_forcing), see the initial states below.
      seed_reference = seed
      CALL read_ensemble_params('parameter.inp','SUAVE',                 &
           num_realizations,kndg,nudge_lambda,ensemble_seed_base,        &
           seed_forcing,ensemble_ic,ensemble_init,newmembers)
      IF (nudge_lambda.le.0.0_GP) CALL sync_abort('nudge_lambda must be > 0')
      IF (nudge_lambda*dt.ge.1.0_GP) CALL sync_abort('nudge_lambda*dt '// &
          'must be < 1 (the relaxation would overshoot every step)')
      IF (myrank.eq.0) WRITE(*,'(A,ES14.6,A,ES14.6,A,ES14.6)')           &
          '   nudge_lambda     : ',nudge_lambda,'  tau=1/lambda=',        &
          1.0_GP/nudge_lambda,'  lambda*dt=',nudge_lambda*dt
! <<< SYNCHRO

! I/O initialization
      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL io_init(myrank,(/nx,ny,nz/),ksta,kend,planio)

! Now we can initialize the PDE methods and read the particles status
      fluid        = init_pdes_from_file('parameter.inp')
      if (dopart) particle = init_particles_from_file('parameter.inp')
      CALL workspace%initialize_pool(NUMTMPREAL,NUMTMPCOMP,NUMTMPPART)
      CALL workspace%init_host_entries(NUMTMPHREAL,NUMTMPHCOMP)
      CALL fluid%Solver_ctor('parameter.inp',workspace,planio)
      num_components = fluid%state_size()
      CALL GState_alloc(field    , num_components)
      CALL GState_alloc(field_nxt, num_components)
      CALL GState_alloc(force    , num_components)
      iclist       = init_ic_from_file(     'parameter.inp')
      forcemethod  = init_forcing_from_file('parameter.inp',workspace)

! >>> SYNCHRO
! Ensemble members: solver, IC chain and states from parameterNNN.inp.
! They share the workspace pool, the I/O plan and the forcing of the
! reference. MOVE_ALLOC avoids copies of the polymorphic objects.
      ALLOCATE(ensemble(num_realizations))
      DO ir = 1,num_realizations
         WRITE(outlabel,'(I3.3)') ir
         ensemble(ir)%infile = 'parameter'//TRIM(outlabel)//'.inp'
         BLOCK
            CLASS(EquationBase), ALLOCATABLE :: tmp_pde
            TYPE (icChain)     , ALLOCATABLE :: tmp_ic(:)
            tmp_pde = init_pdes_from_file(TRIM(ensemble(ir)%infile))
            CALL MOVE_ALLOC(tmp_pde,ensemble(ir)%pde)
            CALL ensemble(ir)%pde%Solver_ctor(TRIM(ensemble(ir)%infile),  &
                                              workspace,planio)
            tmp_ic = init_ic_from_file(TRIM(ensemble(ir)%infile))
            CALL MOVE_ALLOC(tmp_ic,ensemble(ir)%iclist)
         END BLOCK
         CALL check_member(ensemble(ir),fluid)
         CALL GState_alloc(ensemble(ir)%field    ,num_components)
         CALL GState_alloc(ensemble(ir)%field_nxt,num_components)
         CALL GState_alloc(ensemble(ir)%force_nudged,num_components)
      END DO
! <<< SYNCHRO

! Initialization of the numerical domain
      CALL box_init('parameter.inp')

! Initialization of the particles
      if (dopart) then
         CALL particle%part_ctor('parameter.inp',fluid,workspace,part,part_nxt)
         icplist = init_icp_from_file('parameter.inp')
      endif

! Initializes the FFT library. This must be done at this
! stage as it requires status and benchmark initialization.
! Use FFTW_ESTIMATE or FFTW_MEASURE in short runs
! Use FFTW_PATIENT or FFTW_EXHAUSTIVE in very long runs
! FFTW_MEASURE is used below: at 128^3 on two MPI ranks it makes the FFTs
! about 23% faster than FFTW_ESTIMATE, 13 to 18% of the whole time step.
      nth = 1
!$    nth = omp_get_max_threads()
!$    CALL fftp3d_init_threads(ierr)
      CALL GTBenchInit(bench,ihcpu1,ihomp1,ihwtm1,ihcpu2,ihomp2,ihwtm2)
      IF (bench.eq.2) THEN
         CALL GTStart(ihcpu2); CALL GTStart(ihomp2); CALL GTStart(ihwtm2)
      ENDIF
      CALL fftp3d_create_plan(planrc,(/nx,ny,nz/),FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp3d_create_plan(plancr,(/nx,ny,nz/),FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)
      IF (bench.eq.2) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu2);  CALL GTStop(ihomp2);  CALL GTStop(ihwtm2)
      ENDIF

! Initial states
 IC : IF (stat.eq.0) THEN                 ! If stat=0 we start a new run
        ini  = 1
        sind = 0                          ! index for the spectrum
        tind = 0                          ! index for the binaries
        pind = 0                          ! index for the particles
        timet = tstep
        timep = pstep
        timec = cstep
        times = sstep
        timef = 0
      ELSE                    ! If stat.ne.0 a previous run is continued
        ini = int((stat-1)*tstep) + 1
        tind = int(stat)
        sind = int(real(ini,kind=GP)/real(sstep,kind=GP)+1)
        pind = int((stat-1)*lgmult+1)
        timet = 0
        timep = 0
        timec = int(modulo(float(ini-1),float(cstep)))
        times = int(modulo(float(ini-1),float(sstep)))
        timef = int(modulo(float(ini-1),float(fstep)))
      ENDIF IC
! >>> SYNCHRO
      seed = seed_reference               ! ICs of the reference: status seed
      CALL zero_state_host(field)         ! random ICs only fill their band
! <<< SYNCHRO
      CALL init_allstates(iclist,fluid,field)
      CALL GState_update_to(field)        ! Device copies (no-op in host builds).
      CALL GState_copy(field_nxt,field)   ! nxt is used by I/O and all steppers.
! >>> SYNCHRO
! The forcing is always built from seed_forcing, whatever the ICs drew
! from the RNG: a restart rebuilds the same (static) forcing.
      seed = seed_forcing
      CALL zero_state_host(force)
! <<< SYNCHRO
      CALL init_forcing(forcemethod,fluid,force)
! >>> SYNCHRO
      seed_after = seed
      seed       = seed_reference
! <<< SYNCHRO
      CALL GState_update_to(force)
      CALL fluid%sync_device()            ! Auxiliary arrays of the solver (if any).
      fstatic = forcing_is_static(forcemethod)
! >>> SYNCHRO
      CALL forcing_report(force,seed_forcing,seed_after,fstatic,kndg,    &
                          fluid%todir_)
! <<< SYNCHRO
      if (dopart) then
         CALL init_allpstates(icplist,fluid,field,particle,part)
         if (size(part(1)%rcomp) .ne. size(part_nxt(1)%rcomp)) then
            call GPState_resize(part_nxt,particle%partbuff_) ! We resize part_nxt
         endif
         CALL GPState_update_to(part)     ! Device copies of the particle
         CALL particle%sync_device()      ! states and of the class arrays.
         CALL GPState_copy(part_nxt,part) ! We also update part_nxt.
      endif

! >>> SYNCHRO
! Initial states of the members (built on the host copies, then copied
! to the device), and their time steppers.
!  new members, 'scramble': u_m = u_ref, random phases in shells > kndg
!  new members, 'iclist'  : member ICs (seed ensemble_seed_base+10007*ir)
!                           and shells <= kndg replaced by the reference
!  continued members      : member files (read_v) as written, with no
!                           replacement: the nudged shells are part of
!                           the state, and the run continues smoothly
      DO ir = 1,num_realizations
         IF (newmembers) THEN
            IF (TRIM(ensemble_ic).eq.'scramble') THEN
               CALL copy_state_host(ensemble(ir)%field,field)
               CALL scramble_phases_above(ensemble(ir)%field,kndg,       &
                                          ensemble_seed_base+10007*ir)
            ELSE
               seed = ensemble_seed_base+10007*ir
               CALL zero_state_host(ensemble(ir)%field)
               CALL init_allstates(ensemble(ir)%iclist,ensemble(ir)%pde, &
                                   ensemble(ir)%field)
               CALL replace_scales(field,ensemble(ir)%field,kndg)
            ENDIF
         ELSE
            CALL require_read_v(ensemble(ir))
            CALL zero_state_host(ensemble(ir)%field)
            CALL init_allstates(ensemble(ir)%iclist,ensemble(ir)%pde,    &
                                ensemble(ir)%field)
         ENDIF
         CALL GState_update_to(ensemble(ir)%field)
         CALL ensemble(ir)%pde%sync_device()
         CALL GState_copy(ensemble(ir)%field_nxt,ensemble(ir)%field)
         CALL build_stepper_from_file(TRIM(ensemble(ir)%infile),         &
                                      ensemble(ir)%stepper,workspace,    &
                                      ensemble(ir)%pde)
      END DO
      seed = seed_reference
! <<< SYNCHRO

! Sets up the time stepper
      if (dopart) then
        CALL build_stepper_from_file('parameter.inp',stepper,workspace,fluid,particle)
      else
        CALL build_stepper_from_file('parameter.inp',stepper,workspace,fluid)
      endif

! Time integration scheme starts here. In offload builds the fields are
! worked on the device inside the time step (gdev_active set) and on
! their host copies elsewhere (for I/O and diagnostics): the host copies
! of the fields are refreshed before any output, and the forcing, 
! which is computed on the host, is copied to the device in the steps
! in which an update method modifies it or its building blocks (e.g.,
! in shuffle schemes the blending is computed over the device copies).
! If we are doing a benchmark, we measure cputime before starting. We
! also re-inititialize the fftp timers.
      IF (bench.eq.1) THEN
         ffttime  = 0.D00; tratime  = 0.0D0; comtime  = 0.D00; tottime  = 0.0D0
         CALL GTStart(ihcpu1); CALL GTStart(ihomp1); CALL GTStart(ihwtm1)
      ENDIF

 RK : DO t = ini,step
         time = (t-1)*dt
         ! Refreshes the host copies of the fields if any output is due
         IF (((timet.eq.tstep).or.(timec.eq.cstep).or.(times.eq.sstep) &
              .or.(dopart.and.(timep.eq.pstep))).and.(bench.eq.0)) THEN
            CALL GState_update_from(field_nxt)
! >>> SYNCHRO
            DO ir = 1,num_realizations
               CALL GState_update_from(ensemble(ir)%field_nxt)
            END DO
! <<< SYNCHRO
         ENDIF
         ! Every 'tstep' steps, stores the fields in binary files
         IF ((timet.eq.tstep).and.(bench.eq.0)) THEN
            timet = 0
            tind = tind+1
            CALL fluid%write_states(field_nxt, planio)
! >>> SYNCHRO
            DO ir = 1,num_realizations
               CALL ensemble(ir)%pde%write_states(ensemble(ir)%field_nxt,planio)
            END DO
! <<< SYNCHRO
         ENDIF

         ! Every 'pstep' steps, stores the particle states
         if (dopart) then
            IF ((timep.eq.pstep).and.(bench.eq.0)) THEN
               timep = 0
               pind = pind+1
               CALL GPState_update_from(part_nxt) ! host copies for the I/O
               CALL particle%write_pstate(time,fluid,field_nxt,part_nxt)
            ENDIF
         endif

         ! Every 'cstep' steps writes global quantities (the host copy
         ! of a time dependent forcing is refreshed: its updates may
         ! have run on the device)
         IF ((timec.eq.cstep).and.(bench.eq.0)) THEN
            timec = 0
            IF (.not.fstatic) CALL GState_update_from(force)
            CALL fluid%global(field_nxt, force, t)
! >>> SYNCHRO
            DO ir = 1,num_realizations
               CALL member_global(ensemble(ir),field_nxt,force,t)
            END DO
! <<< SYNCHRO
         ENDIF

         ! Every 'sstep' steps writes spectra
         IF ((times.eq.sstep).and.(bench.eq.0)) THEN
            times = 0
            sind = sind+1
            CALL fluid%spectra(field_nxt)
! >>> SYNCHRO
            DO ir = 1,num_realizations
               CALL member_spectra(ensemble(ir),field_nxt)
            END DO
! <<< SYNCHRO
         ENDIF

         ! Time evolution
         CALL update_forcing(forcemethod,fluid,force)
         IF (forcing_changed(forcemethod)) CALL GState_update_to(force)
         gdev_active = .TRUE.
! >>> SYNCHRO
! Forcing of each member: the forcing of the reference plus the nudging
! term, evaluated with the states at the beginning of the step (u_ref
! and u_member are still in field_nxt), on the device copies in offload
! builds. Done before any step, so it does not depend on the stepper.
         DO ir = 1,num_realizations
            CALL GState_copy(ensemble(ir)%force_nudged,force)
            CALL add_nudge_forcing(field_nxt,ensemble(ir)%field_nxt,kndg, &
                                   nudge_lambda,ensemble(ir)%force_nudged)
         END DO
! <<< SYNCHRO
         if (dopart) then
            CALL GState_copy(field,field_nxt)
            CALL GPState_copy(part,part_nxt)
            CALL stepper%gstep(time, field, part, force, dt, field_nxt, part_nxt)
         else
            CALL GState_copy(field,field_nxt)
            CALL stepper%gstep(time, field, force, dt, field_nxt)
         endif
! >>> SYNCHRO
! Each member takes the same step with its nudged forcing. No modes are
! overwritten: the synchronization is only the relaxation term.
         DO ir = 1,num_realizations
            CALL GState_copy(ensemble(ir)%field,ensemble(ir)%field_nxt)
            CALL ensemble(ir)%stepper%gstep(time,ensemble(ir)%field,     &
                                            ensemble(ir)%force_nudged,dt,&
                                            ensemble(ir)%field_nxt)
         END DO
! <<< SYNCHRO
         gdev_active = .FALSE.
         timet = timet+1; timep = timep+1; timec = timec+1; times = times+1
      END DO RK

! Finishes and writes the benchmark results
      IF (bench.gt.0) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL GTStop(ihcpu1); CALL GTStop(ihomp1); CALL GTStop(ihwtm1)
      ENDIF
      CALL GTBenchReport(bench,myrank,ini,step,ihcpu1,ihomp1,ihwtm1,             &
               ihcpu2,ihomp2,ihwtm2)
      IF (dopart) THEN
        rbal = rbal + GetLoadBal(particle) ! Get load balancing
        CALL GTBenchReportParticles(bench, myrank, ini, step, maxparts, rbal,    &
             [ GetTime(particle,GPTIME_STEP),   GetTime(particle,GPTIME_COMM),   &
               GetTime(particle,GPTIME_SPLINE), GetTime(particle,GPTIME_TRANSP), &
               GetTime(particle,GPTIME_DATAEX), GetTime(particle,GPTIME_INTERP), &
               GetTime(particle,GPTIME_PUPDATE) ], 'gpbenchmark.txt')
      ENDIF

! End of main
      CALL workspace%report_peaks()
      CALL GTFree(ihcpu1)
      CALL GTFree(ihomp1)
      CALL GTFree(ihwtm1)
      CALL GTFree(ihcpu2)
      CALL GTFree(ihomp2)
      CALL GTFree(ihwtm2)
      CALL MPI_FINALIZE(ierr)
      CALL fftp3d_destroy_plan(plancr)
      CALL fftp3d_destroy_plan(planrc)
      
      END PROGRAM MAIN
