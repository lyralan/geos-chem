!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fullchem_lb_mod.F90
!
! !DESCRIPTION: Module fullchem\_LB\_Mod balances the KPP integration work of
!  the fullchem mechanism across MPI ranks.  Do\_FullChem then makes two
!  passes over the grid: pass 1 prepares each grid box as usual and stores
!  the KPP inputs (C, RCONST, ICNTRL, RCNTRL) here; LB\_Run moves part of the
!  integration work from ranks with more work than the mean to ranks with
!  less, integrates every box (including the retry after a failed first
!  attempt), and returns the results to the owning rank; pass 2 finishes
!  each box with these results.  Work is weighted by the wall time each box
!  took in the previous chemistry step.
!\\
!\\
!  Integrate depends only on C, RCONST, ICNTRL, RCNTRL, TIN, TOUT and the
!  tolerances (ICNTRL(15) = -1: no rate updates inside the integrator), so a
!  box integrated on another rank gives the same answer bit for bit.
!\\
!\\
!  Turned on by "operations:chemistry:load_balance: true" in
!  geoschem_config.yml (GCHP only).
!\\
!\\
! !INTERFACE:
!
MODULE FullChem_LB_Mod
!
! !USES:
!
  USE gckpp_Parameters, ONLY : NSPEC, NREACT
  USE gckpp_Precision,  ONLY : dp
  USE Precision_Mod
#if defined( MODEL_GCHP )
  USE mpi
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: LB_Begin
  PUBLIC  :: LB_Store
  PUBLIC  :: LB_IsActive
  PUBLIC  :: LB_Restore
  PUBLIC  :: LB_Result
  PUBLIC  :: LB_Neg0
  PUBLIC  :: LB_Time
  PUBLIC  :: LB_Run
!
! !REVISION HISTORY:
!  27 Oct 2026 - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Sizes of one box's record in the work and result messages [REAL(dp)]
  INTEGER,  PARAMETER   :: SZIN  = NSPEC + NREACT + 40
  INTEGER,  PARAMETER   :: SZOUT = NSPEC + 80 + 6

  ! Message tags
  INTEGER,  PARAMETER   :: TAG_WORK = 7701, TAG_RESULT = 7702

  ! Grid dimensions and number of boxes on this rank
  INTEGER               :: NX_ = 0, NY_ = 0, NZ_ = 0, NBOX = 0

  ! Per-box inputs (stored in pass 1)
  LOGICAL,  ALLOCATABLE :: Active(:)
  REAL(dp), ALLOCATABLE :: Cin(:,:), Rin(:,:), Rcn(:,:)
  INTEGER,  ALLOCATABLE :: Icn(:,:)

  ! Per-box results (read in pass 2)
  REAL(dp), ALLOCATABLE :: Cout(:,:), Rst1(:,:), Rst2(:,:)
  INTEGER,  ALLOCATABLE :: Ist1(:,:), Ist2(:,:)
  INTEGER,  ALLOCATABLE :: Ierr1(:), Ierr2(:), nNeg(:), rNV(:), cNZ(:)

  ! Wall time of each box's integration [s]: from the previous chemistry
  ! step it is the work estimate for the next one
  REAL(dp), ALLOCATABLE :: Cost(:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Begin
!
! !DESCRIPTION: Allocates the per-box arrays on first use and clears the
!  active flags at the start of a chemistry step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LB_Begin( NX, NY, NZ )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: NX, NY, NZ   ! Local grid dimensions
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( NBOX /= NX * NY * NZ ) THEN
       IF ( ALLOCATED( Active ) ) THEN
          DEALLOCATE( Active, Cin, Rin, Rcn, Icn, Cout, Rst1, Rst2, Ist1,    &
                      Ist2, Ierr1, Ierr2, nNeg, rNV, cNZ, Cost )
       ENDIF
       NX_  = NX
       NY_  = NY
       NZ_  = NZ
       NBOX = NX * NY * NZ
       ALLOCATE( Active(NBOX),          Cin(NSPEC,NBOX), Rin(NREACT,NBOX),   &
                 Rcn(20,NBOX),          Icn(20,NBOX),    Cout(NSPEC,NBOX),   &
                 Rst1(20,NBOX),         Rst2(20,NBOX),   Ist1(20,NBOX),      &
                 Ist2(20,NBOX),         Ierr1(NBOX),     Ierr2(NBOX),        &
                 nNeg(NBOX),            rNV(NBOX),       cNZ(NBOX),          &
                 Cost(NBOX)                                                 )
       Cost = 0.0_dp
    ENDIF
    Active = .FALSE.

  END SUBROUTINE LB_Begin
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Store
!
! !DESCRIPTION: Pass 1: stores the KPP inputs of box (I,J,L) as they are
!  just before the call to Integrate.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LB_Store( I, J, L, ICNTRL, RCNTRL )
!
! !USES:
!
    USE gckpp_Global, ONLY : C, RCONST
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I, J, L
    INTEGER,  INTENT(IN) :: ICNTRL(20)
    REAL(dp), INTENT(IN) :: RCNTRL(20)
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: K

    K           = BoxIndex( I, J, L )
    Active(K)   = .TRUE.
    Cin(:,K)    = C
    Rin(:,K)    = RCONST
    Icn(:,K)    = ICNTRL
    Rcn(:,K)    = RCNTRL

  END SUBROUTINE LB_Store
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_IsActive
!
! !DESCRIPTION: Returns true if box (I,J,L) was stored in pass 1.
!\\
!\\
! !INTERFACE:
!
  LOGICAL FUNCTION LB_IsActive( I, J, L )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I, J, L
!EOP
!------------------------------------------------------------------------------
!BOC
    LB_IsActive = Active( BoxIndex( I, J, L ) )

  END FUNCTION LB_IsActive
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Restore
!
! !DESCRIPTION: Pass 2: restores C, RCONST, ICNTRL and RCNTRL of box (I,J,L)
!  as they were just before the call to Integrate.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LB_Restore( I, J, L, ICNTRL, RCNTRL )
!
! !USES:
!
    USE gckpp_Global, ONLY : C, RCONST
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: I, J, L
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: ICNTRL(20)
    REAL(dp), INTENT(OUT) :: RCNTRL(20)
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: K

    K      = BoxIndex( I, J, L )
    C      = Cin(:,K)
    RCONST = Rin(:,K)
    ICNTRL = Icn(:,K)
    RCNTRL = Rcn(:,K)

  END SUBROUTINE LB_Restore
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Result
!
! !DESCRIPTION: Pass 2: returns what the first (Attempt = 1) or second
!  (Attempt = 2) call to Integrate gave for box (I,J,L), as if Integrate had
!  been called here.  C is set to the final concentrations (the first
!  attempt's when there was no second one).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LB_Result( Attempt, I, J, L, ISTATUS, RSTATE, IERR )
!
! !USES:
!
    USE gckpp_Global, ONLY : C, rNVAR, cNONZERO
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: Attempt, I, J, L
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: ISTATUS(20)
    REAL(dp), INTENT(OUT) :: RSTATE(20)
    INTEGER,  INTENT(OUT) :: IERR
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: K

    K = BoxIndex( I, J, L )
    C = Cout(:,K)
    IF ( Attempt == 1 ) THEN
       ISTATUS  = Ist1(:,K)
       RSTATE   = Rst1(:,K)
       IERR     = Ierr1(K)
       rNVAR    = rNV(K)
       cNONZERO = cNZ(K)
    ELSE
       ISTATUS  = Ist2(:,K)
       RSTATE   = Rst2(:,K)
       IERR     = Ierr2(K)
    ENDIF

  END SUBROUTINE LB_Result
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Neg0
!
! !DESCRIPTION: Number of negative concentrations after the first call to
!  Integrate for box (I,J,L) (KppNegatives0 diagnostic).
!\\
!\\
! !INTERFACE:
!
  INTEGER FUNCTION LB_Neg0( I, J, L )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I, J, L
!EOP
!------------------------------------------------------------------------------
!BOC
    LB_Neg0 = nNeg( BoxIndex( I, J, L ) )

  END FUNCTION LB_Neg0
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Time
!
! !DESCRIPTION: Wall time [s] of the integration of box (I,J,L), wherever it
!  ran (KppTime diagnostic).
!\\
!\\
! !INTERFACE:
!
  REAL(fp) FUNCTION LB_Time( I, J, L )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I, J, L
!EOP
!------------------------------------------------------------------------------
!BOC
    LB_Time = Cost( BoxIndex( I, J, L ) )

  END FUNCTION LB_Time
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: LB_Run
!
! !DESCRIPTION: Integrates all stored boxes, moving work between ranks so
!  that each rank's estimated integration time is close to the mean.  The
!  plan is computed identically on every rank from the gathered per-rank
!  work estimates: ranks above the mean (donors) are paired with ranks below
!  it (receivers), largest first; each donor sends its most expensive boxes
!  until the agreed amount of work is covered.  Every rank integrates the
!  boxes it kept first and the boxes it received next.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LB_Run( Comm, TIN, TOUT, UseAR, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: Comm        ! MPI communicator of the chem ranks
    REAL(dp), INTENT(IN)  :: TIN, TOUT   ! Integration interval [s]
    LOGICAL,  INTENT(IN)  :: UseAR       ! Auto-reduce solver on?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( MODEL_GCHP )
    ! Transfers smaller than this fraction of the mean load are not worth
    ! the messages
    REAL(dp), PARAMETER   :: MIN_FRAC = 0.005_dp

    INTEGER               :: nRank, me, ierr, K, B, R, S, T, nAct
    INTEGER               :: nTr, iD, iR, nDon, nRec, idx, nDone
    REAL(dp)              :: wDef, myLoad, Mean, Amt, Got
    INTEGER,  ALLOCATABLE :: Act(:), Ord(:), Dest(:), DonR(:), RecR(:)
    INTEGER,  ALLOCATABLE :: TrD(:), TrR(:), nSend(:), nRecv(:)
    INTEGER,  ALLOCATABLE :: SendList(:,:), Req(:), ReqW(:), ReqR(:), ReqS(:)
    REAL(dp), ALLOCATABLE :: W(:), Load(:), Ex(:), De(:), TrA(:)
    REAL(dp), ALLOCATABLE :: WorkS(:,:), WorkR(:,:), ResS(:,:), ResR(:,:)
    INTEGER,  ALLOCATABLE :: OffS(:), OffR(:)

    RC = GC_SUCCESS
    CALL MPI_Comm_size( Comm, nRank, ierr )
    CALL MPI_Comm_rank( Comm, me,    ierr )

    !------------------------------------------------------------------
    ! Work estimate of each active box: last step's time; boxes without
    ! one (first step) get the mean of those that have one, or 1
    !------------------------------------------------------------------
    nAct = COUNT( Active )
    ALLOCATE( Act(nAct), W(nAct), Dest(nAct) )
    Act  = PACK( [ ( K, K = 1, NBOX ) ], Active )
    Dest = me
    W    = Cost(Act)
    IF ( ANY( W > 0.0_dp ) ) THEN
       wDef = SUM( W, MASK = W > 0.0_dp ) / COUNT( W > 0.0_dp )
    ELSE
       wDef = 1.0_dp
    ENDIF
    WHERE( W <= 0.0_dp ) W = wDef
    myLoad = SUM( W )

    ALLOCATE( Load(0:nRank-1) )
    CALL MPI_Allgather( myLoad, 1, MPI_DOUBLE_PRECISION,                     &
                        Load,   1, MPI_DOUBLE_PRECISION, Comm, ierr )
    Mean = SUM( Load ) / nRank

    !------------------------------------------------------------------
    ! Transfer plan (same on every rank)
    !------------------------------------------------------------------
    nDon = COUNT( Load > Mean )
    nRec = COUNT( Load < Mean )
    ALLOCATE( DonR(nDon), RecR(nRec), Ex(nDon), De(nRec) )
    DonR = PACK( [ ( R, R = 0, nRank-1 ) ], Load > Mean )
    RecR = PACK( [ ( R, R = 0, nRank-1 ) ], Load < Mean )
    CALL SortDesc( Load(DonR) - Mean, DonR )
    CALL SortDesc( Mean - Load(RecR), RecR )
    Ex = Load(DonR) - Mean
    De = Mean - Load(RecR)

    ALLOCATE( TrD(nDon+nRec), TrR(nDon+nRec), TrA(nDon+nRec) )
    nTr = 0
    iD  = 1
    iR  = 1
    DO WHILE ( iD <= nDon .and. iR <= nRec )
       Amt = MIN( Ex(iD), De(iR) )
       IF ( Amt > MIN_FRAC * Mean ) THEN
          nTr      = nTr + 1
          TrD(nTr) = DonR(iD)
          TrR(nTr) = RecR(iR)
          TrA(nTr) = Amt
       ENDIF
       Ex(iD) = Ex(iD) - Amt
       De(iR) = De(iR) - Amt
       IF ( Ex(iD) <= MIN_FRAC * Mean ) iD = iD + 1
       IF ( De(iR) <= MIN_FRAC * Mean ) iR = iR + 1
    ENDDO

    !------------------------------------------------------------------
    ! As a donor: pick the boxes for each transfer, most expensive first
    !------------------------------------------------------------------
    ALLOCATE( Ord(nAct) )
    Ord = [ ( B, B = 1, nAct ) ]
    CALL SortDesc( W, Ord )
    idx = 1
    DO T = 1, nTr
       IF ( TrD(T) /= me ) CYCLE
       Got = 0.0_dp
       DO WHILE ( idx <= nAct .and. Got < TrA(T) - 0.5_dp * W(Ord(idx)) )
          Dest(Ord(idx)) = TrR(T)
          Got = Got + W(Ord(idx))
          idx = idx + 1
       ENDDO
    ENDDO

    !------------------------------------------------------------------
    ! Counts, then post the receives of work and results
    !------------------------------------------------------------------
    ALLOCATE( nSend(0:nRank-1), nRecv(0:nRank-1) )
    nSend = 0
    DO B = 1, nAct
       IF ( Dest(B) /= me ) nSend(Dest(B)) = nSend(Dest(B)) + 1
    ENDDO
    CALL MPI_Alltoall( nSend, 1, MPI_INTEGER, nRecv, 1, MPI_INTEGER,         &
                       Comm,  ierr )

    ALLOCATE( OffS(0:nRank), OffR(0:nRank) )
    OffS(0) = 0
    OffR(0) = 0
    DO R = 0, nRank-1
       OffS(R+1) = OffS(R) + nSend(R)
       OffR(R+1) = OffR(R) + nRecv(R)
    ENDDO
    ALLOCATE( SendList(MAX(OffS(nRank),1),1) )
    ALLOCATE( WorkS(SZIN, MAX(OffS(nRank),1)),  ResR(SZOUT, MAX(OffS(nRank),1)) )
    ALLOCATE( WorkR(SZIN, MAX(OffR(nRank),1)),  ResS(SZOUT, MAX(OffR(nRank),1)) )
    ALLOCATE( Req(0:nRank-1), ReqW(0:nRank-1), ReqR(0:nRank-1), ReqS(0:nRank-1) )
    Req  = MPI_REQUEST_NULL     ! work sent
    ReqW = MPI_REQUEST_NULL     ! work received
    ReqR = MPI_REQUEST_NULL     ! results received
    ReqS = MPI_REQUEST_NULL     ! results sent

    ! Pack and send work; post receives for its results
    S = 0
    DO R = 0, nRank-1
       IF ( nSend(R) == 0 ) CYCLE
       DO B = 1, nAct
          IF ( Dest(B) /= R ) CYCLE
          S = S + 1
          SendList(S,1) = Act(B)
          CALL PackWork( Act(B), WorkS(:,S) )
       ENDDO
       CALL MPI_Isend( WorkS(1,OffS(R)+1), SZIN*nSend(R),                    &
                       MPI_DOUBLE_PRECISION, R, TAG_WORK, Comm, Req(R), ierr )
       CALL MPI_Irecv( ResR(1,OffS(R)+1), SZOUT*nSend(R),                    &
                       MPI_DOUBLE_PRECISION, R, TAG_RESULT, Comm, ReqR(R),   &
                       ierr )
    ENDDO

    ! Post receives for incoming work
    DO R = 0, nRank-1
       IF ( nRecv(R) == 0 ) CYCLE
       CALL MPI_Irecv( WorkR(1,OffR(R)+1), SZIN*nRecv(R),                    &
                       MPI_DOUBLE_PRECISION, R, TAG_WORK, Comm, ReqW(R), ierr )
    ENDDO

    !------------------------------------------------------------------
    ! Integrate the boxes kept here
    !------------------------------------------------------------------
    DO B = 1, nAct
       IF ( Dest(B) /= me ) CYCLE
       K = Act(B)
       CALL SolveBox( Cin(:,K), Rin(:,K), Icn(:,K), Rcn(:,K), TIN, TOUT,     &
                      UseAR,    ResS(:,1) )
       CALL UnpackResult( K, ResS(:,1) )
    ENDDO

    !------------------------------------------------------------------
    ! Integrate received boxes as their messages arrive; send results back
    !------------------------------------------------------------------
    nDone = 0
    DO WHILE ( nDone < COUNT( nRecv > 0 ) )
       CALL MPI_Waitany( nRank, ReqW, idx, MPI_STATUS_IGNORE, ierr )
       R = idx - 1
       nDone = nDone + 1
       DO S = OffR(R)+1, OffR(R+1)
          CALL SolveBox( WorkR(1:NSPEC,S),                                   &
                         WorkR(NSPEC+1:NSPEC+NREACT,S),                      &
                         NINT( WorkR(NSPEC+NREACT+1:NSPEC+NREACT+20,S) ),    &
                         WorkR(NSPEC+NREACT+21:NSPEC+NREACT+40,S),           &
                         TIN, TOUT, UseAR, ResS(:,S) )
       ENDDO
       CALL MPI_Isend( ResS(1,OffR(R)+1), SZOUT*nRecv(R),                    &
                       MPI_DOUBLE_PRECISION, R, TAG_RESULT, Comm, ReqS(R),   &
                       ierr )
    ENDDO

    !------------------------------------------------------------------
    ! Collect the results of the boxes sent away
    !------------------------------------------------------------------
    CALL MPI_Waitall( nRank, ReqR, MPI_STATUSES_IGNORE, ierr )
    DO S = 1, OffS(nRank)
       CALL UnpackResult( SendList(S,1), ResR(:,S) )
    ENDDO
    CALL MPI_Waitall( nRank, Req,  MPI_STATUSES_IGNORE, ierr )
    CALL MPI_Waitall( nRank, ReqS, MPI_STATUSES_IGNORE, ierr )

    DEALLOCATE( Act, W, Dest, Load, DonR, RecR, Ex, De, TrD, TrR, TrA, Ord,  &
                nSend, nRecv, OffS, OffR, SendList, WorkS, ResR, WorkR,      &
                ResS, Req, ReqW, ReqR, ReqS )
#else
    RC = GC_FAILURE
    CALL GC_Error( 'Chemistry load balancing is only available in GCHP',     &
                   RC, 'LB_Run (in GeosCore/fullchem_lb_mod.F90)' )
#endif

  END SUBROUTINE LB_Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SolveBox
!
! !DESCRIPTION: Integrates one box exactly as Do\_FullChem does: one call to
!  Integrate and, if it fails, a second one from the same concentrations
!  with Hstart = 0 and auto-reduce off.  Writes the result record Res.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SolveBox( Ci, Ri, Ic, Rc, TIN, TOUT, UseAR, Res )
!
! !USES:
!
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetKeepActive
    USE gckpp_Global,             ONLY : C, RCONST, rNVAR, cNONZERO
    USE gckpp_Integrator,         ONLY : Integrate
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: Ci(NSPEC), Ri(NREACT), Rc(20), TIN, TOUT
    INTEGER,  INTENT(IN)  :: Ic(20)
    LOGICAL,  INTENT(IN)  :: UseAR
!
! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(OUT) :: Res(SZOUT)
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER  :: ICNTRL(20), ISTATUS(20), IERR
    REAL(dp) :: RCNTRL(20), RSTATE(20), T0
    INTEGER  :: O

    T0     = Wtime()
    C      = Ci
    RCONST = Ri
    ICNTRL = Ic
    RCNTRL = Rc
    CALL fullchem_AR_SetKeepActive( option=.TRUE. )

    ! Result record: C | ISTATUS 1 | RSTATE 1 | ISTATUS 2 | RSTATE 2 |
    !                IERR 1, IERR 2, nNeg, rNVAR, cNONZERO, time
    Res     = 0.0_dp
    ISTATUS = 0
    RSTATE  = 0.0_dp
    IERR    = 0
    CALL Integrate( TIN,    TOUT,    ICNTRL,                                 &
                    RCNTRL, ISTATUS, RSTATE, IERR                           )
    O = NSPEC
    Res(O+1 :O+20) = ISTATUS
    Res(O+21:O+40) = RSTATE
    Res(O+81)      = IERR
    Res(O+83)      = COUNT( C < 0.0_dp )
    Res(O+84)      = rNVAR
    Res(O+85)      = cNONZERO

    IF ( IERR < 0 ) THEN
       RCNTRL(3) = 0.0_dp
       C         = Ci
       RCONST    = Ri
       IF ( UseAR ) RCNTRL(12) = -1.0_dp
       CALL Integrate( TIN,    TOUT,    ICNTRL,                              &
                       RCNTRL, ISTATUS, RSTATE, IERR                        )
       Res(O+41:O+60) = ISTATUS
       Res(O+61:O+80) = RSTATE
       Res(O+82)      = IERR
    ENDIF

    Res(1:NSPEC) = C
    Res(O+86)    = Wtime() - T0

  END SUBROUTINE SolveBox
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PackWork
!
! !DESCRIPTION: Writes the work record of box K.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PackWork( K, Work )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: K
!
! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(OUT) :: Work(SZIN)
!EOP
!------------------------------------------------------------------------------
!BOC
    Work(1:NSPEC)                       = Cin(:,K)
    Work(NSPEC+1:NSPEC+NREACT)          = Rin(:,K)
    Work(NSPEC+NREACT+1:NSPEC+NREACT+20)  = Icn(:,K)
    Work(NSPEC+NREACT+21:NSPEC+NREACT+40) = Rcn(:,K)

  END SUBROUTINE PackWork
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: UnpackResult
!
! !DESCRIPTION: Stores the result record of box K and its time as the work
!  estimate for the next step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE UnpackResult( K, Res )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: K
    REAL(dp), INTENT(IN) :: Res(SZOUT)
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: O

    O         = NSPEC
    Cout(:,K) = Res(1:NSPEC)
    Ist1(:,K) = NINT( Res(O+1 :O+20) )
    Rst1(:,K) =       Res(O+21:O+40)
    Ist2(:,K) = NINT( Res(O+41:O+60) )
    Rst2(:,K) =       Res(O+61:O+80)
    Ierr1(K)  = NINT( Res(O+81) )
    Ierr2(K)  = NINT( Res(O+82) )
    nNeg(K)   = NINT( Res(O+83) )
    rNV(K)    = NINT( Res(O+84) )
    cNZ(K)    = NINT( Res(O+85) )
    Cost(K)   = Res(O+86)

  END SUBROUTINE UnpackResult
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SortDesc
!
! !DESCRIPTION: Reorders Idx so that the matching keys are in descending
!  order; ties keep their original order (insertion sort, stable).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SortDesc( Key, Idx )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)    :: Key(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(INOUT) :: Idx(:)
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER               :: A, B, Ti
    REAL(dp)              :: Tk
    REAL(dp), ALLOCATABLE :: Kc(:)

    Kc = Key
    DO A = 2, SIZE( Idx )
       Tk = Kc(A)
       Ti = Idx(A)
       B  = A - 1
       DO WHILE ( B >= 1 )
          IF ( Kc(B) >= Tk ) EXIT
          Kc(B+1)  = Kc(B)
          Idx(B+1) = Idx(B)
          B        = B - 1
       ENDDO
       Kc(B+1)  = Tk
       Idx(B+1) = Ti
    ENDDO

  END SUBROUTINE SortDesc
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BoxIndex
!
! !DESCRIPTION: Linear index of box (I,J,L).
!\\
!\\
! !INTERFACE:
!
  INTEGER FUNCTION BoxIndex( I, J, L )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: I, J, L
!EOP
!------------------------------------------------------------------------------
!BOC
    BoxIndex = I + ( J - 1 ) * NX_ + ( L - 1 ) * NX_ * NY_

  END FUNCTION BoxIndex
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Wtime
!
! !DESCRIPTION: Wall clock time [s].
!\\
!\\
! !INTERFACE:
!
  REAL(dp) FUNCTION Wtime()
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( MODEL_GCHP )
    Wtime = MPI_Wtime()
#else
    INTEGER(8) :: Cnt, Rate
    CALL SYSTEM_CLOCK( Cnt, Rate )
    Wtime = REAL( Cnt, dp ) / REAL( Rate, dp )
#endif

  END FUNCTION Wtime
!EOC
END MODULE FullChem_LB_Mod
