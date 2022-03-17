! 2,4,6 fold
!説明書
!３次元も

      module VAR_REC2D
        type REC2D_PARAMETER
          integer             :: NSYMM            ! symmetry 2/4/6
          integer             :: ND               ! number of data channels for calc
          integer             :: NP               ! number of Compton profiles
          integer             :: NORG             ! number of data channels (original)
          real(kind(1d0))     :: DCH              ! Bin size : (a.u./bin)
          real(kind(1d0))     :: DPZ              ! momentum resolution (a.u.)
          real(kind(1d0))     :: SMPAR            ! smoothing parameter
          integer             :: FLPAR            ! filtering parameter
          integer             :: KXORD            ! Order of B-spline for B(R) inperpolation (theta)
          integer             :: KYORD            ! Order of B-spline for B(R) inperpolation (radius)
          real(kind(1d0))     :: ZMPAR            ! Zoom ratio for output data
          integer             :: NOUT             ! Number of BINs of Output files
          character*256       :: JPZ_FILENAME     ! Input / Compton Profiles (table format)
          character*256       :: RHO_FILENAME     ! Output / 2D rho(pz) table
          character*256       :: ASM_FILENAME     ! Output / 2D rho(pz) asymmetry table
          character*256       :: LCW_FILENAME     ! Output / 2D n(k) table
          integer             :: OUT_FORMAT       ! Output 2D file format  0:table, 1,coulmn

          real(kind(1d0))     :: LTC_CONST_A      ! Lattice constant a-axis (Angstrom)
          real(kind(1d0))     :: LTC_CONST_B      ! Lattice constant b-axis (Angstrom)
          integer             :: NLCW_A           ! Number of BINs of LCW output 
          integer             :: NLCW_B           ! Number of BINs of LCW output
          integer             :: N_FOLDING_A      ! Number of Foldingd 
          integer             :: N_FOLDING_B      ! Number of Foldings
          real(kind(1d0))     :: STEP_FOLDING_A   ! k (ch) 1st brilloiun zone size (ch) 
          real(kind(1d0))     :: STEP_FOLDING_B   ! 
    	  end type
      end module

      program Rec2D
      use CPSEC_INT                               ! CPU time
      use VAR_REC2D
      implicit none

      integer ARGC                                ! commandline parameter
      character ARGV*256
      type(REC2D_PARAMETER) :: P
      real(kind(1d0)), allocatable    :: TH(:)    ! Direction of Pz (Degree)    (NPROF)
      real(kind(1d0)), allocatable    :: PZ(:)    ! pz (a.u.)                   (NDATA)
      real(kind(1d0)), allocatable    :: JPZ(:,:) ! Compton profile      (NDATA, NPROF)
      real(kind(1d0)), allocatable    :: BRQ(:,:) ! B(rq) = FFT(J(Pz))   (NDATA, NPROF)
      real(kind(1d0)), allocatable    :: B_R(:,:) ! B(R)                 (NDATA, NDATA)
      real(kind(1d0)), allocatable    :: R_R(:,:) ! rho(R)               (NDATA, NDATA)
      real(kind(1d0)), allocatable    :: ARR(:,:) ! rho(R) asymmetry     (NDATA, NDATA)
      real(kind(1d0)), allocatable    :: N_K(:,:) ! n(k)                 (NLCW,  NLCW)
!-----------------------------------------------------------------------
      print '("**************************************************************")'
      print '(" **** Start 2D reconstruction *************")'
! Input
      call INIT_PARAMETER(P)
      if (IARGC().EQ.1) then
        call getarg(1, ARGV)
        call READ_INIFILE(P, ARGV)
      else
        call INPUT_PARAMETER(P)
      endif
      call CHECK_PARAMETER(P)

	  
	  
      call WRITE_INIFILE(P, 'INI.tmp')
   	  allocate (TH(P%NP), PZ(P%ND))
      allocate (JPZ(P%ND, P%NP), BRQ(P%ND, P%NP))
      allocate (B_R(P%ND, P%ND), R_R(P%ND, P%ND), ARR(P%ND, P%ND))
      allocate (N_K(-P%NLCW_A:P%NLCW_A, -P%NLCW_B:P%NLCW_B))
      print '(" Read Files                                ")'
      call READ_COMPTON_PROFILE(P, TH, PZ, JPZ)
      call WRITE_PROFILE(P%NP, P%ND, TH, PZ, JPZ, 'JPZ.tmp')

      print '(" Smooth Compton profiles J(pz).            ",$)'
      call SMOOTH_JPZ(P, JPZ)
      call WRITE_PROFILE(P%NP, P%ND, TH, PZ, JPZ, 'JPZ_SM.tmp')
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Fourier Transform J(pz) -> B(rq).         ",$)'
      call FT_JPZ(P, JPZ, BRQ)
      call WRITE_PROFILE(P%NP, P%ND, TH, PZ, BRQ, 'BRQ.tmp')
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Low pass filter to reduce noise of B(rq). ",$)'
      call FILTER_BRQ(P, BRQ)
      call WRITE_PROFILE(P%NP, P%ND, TH, PZ, BRQ, 'BRQ_FL.tmp')
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Convert Axis B(rq) -> B(R).     ",$)'
      call TRANS_BP2BR(P, TH, PZ, BRQ, B_R)
      call WRITE_PROFILE_2D(P%ND, P%NOUT, PZ, B_R, 'B_R.tmp', P%OUT_FORMAT)
      call WRITE_PROFILE_BMP(P%ND, P%NOUT, B_R, 'B_R.bmp')
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Fourier Transform B(R) -> R(R). ",$)'
      call IFT_BR_2D(P, B_R, R_R)
      call WRITE_PROFILE_2D(P%ND, P%NOUT, PZ, R_R, P%RHO_FILENAME, P%OUT_FORMAT)
      call WRITE_PROFILE_BMP(P%ND, P%NOUT, R_R, P%RHO_FILENAME)
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Make asymmetry of R(r).         ",$)'
      call ASYMMETRY(P, R_R, ARR)
      call WRITE_PROFILE_2D(P%ND, P%NOUT, PZ, ARR, P%ASM_FILENAME, P%OUT_FORMAT)
      call WRITE_PROFILE_BMP(P%ND, P%NOUT, ARR, P%ASM_FILENAME)
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" LCW folding.                    ",$)'
      call LCW_FOLDING(P, R_R, N_K)
      call WRITE_N_K_2D(P%NLCW_A, P%NLCW_B, PZ, P%ND, N_K, P%LCW_FILENAME, P%OUT_FORMAT)
      call WRITE_N_K_BMP(P%NLCW_A, P%NLCW_B, N_K, P%LCW_FILENAME)
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

      print '(" Remake Compton J(pz) from R(R). ",$)'
      call REMAKE_JPZ(P, R_R, TH, JPZ)
      call WRITE_PROFILE(P%NP, P%ND, TH, PZ, JPZ, 'JPZ_NEW.tmp')
      print '(A, F8.2, A)', ' done.', CPSEC(), ' sec'

! terminate
      print '("**************************************************************")'
      deallocate (TH, PZ)
      deallocate (JPZ, BRQ, B_R, R_R)
      deallocate (N_K)
      end program Rec2D

      subroutine INIT_PARAMETER(P)
      use VAR_REC2D
      implicit none
      type(REC2D_PARAMETER) :: P
      P%NSYMM = 2
      P%ND = 2048
      P%NP = 0
      P%DCH = 0.1
      P%DPZ = 0.1
      P%SMPAR = -1
      P%FLPAR = 5
      P%KXORD = 3
      P%KYORD = 3
      P%ZMPAR = 1
      P%NOUT = 512
      P%JPZ_FILENAME='*'
      P%RHO_FILENAME = 'RHO.TXT'
      P%LCW_FILENAME = 'LCW.TXT'
      P%ASM_FILENAME = 'ASM.TXT'
      P%LTC_CONST_A = -1.
      P%LTC_CONST_B = -1.
      P%N_FOLDING_A = -1
      P%N_FOLDING_B = -1
      P%OUT_FORMAT = 0
      end subroutine
      
      subroutine READ_INIFILE(P, FILENAME)
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      integer         :: I,J,IOS
      character*256   :: FILENAME
      character*256   :: S, S1, S2

      open (11, FILE = trim(FILENAME), READONLY, STATUS='OLD')
      do I=1, 999
        read (11, *, iostat=IOS) S
        if (IOS.LT.0) exit
        if (S(1:1).EQ.'[') cycle  ! skip comments
        if (S(1:1).EQ.';') cycle  ! skip comments
        J=index(S,'=')
        if (J.EQ.0) cycle        ! skip error line
        S1=S(1:J-1)
        S2=S(J+1:255)
        if (trim(S1).EQ.'NSYMM') read (S2,*) P%NSYMM
        if (trim(S1).EQ.'NP') read (S2,*) P%NP
        if (trim(S1).EQ.'ND') read (S2,*) P%ND
        if (trim(S1).EQ.'DCH') read (S2,*) P%DCH
        if (trim(S1).EQ.'DPZ') read (S2,*) P%DPZ
        if (trim(S1).EQ.'SMPAR') read (S2,*) P%SMPAR
        if (trim(S1).EQ.'FLPAR') read (S2,*) P%FLPAR
        if (trim(S1).EQ.'ZMPAR') read (S2,*) P%ZMPAR
        if (trim(S1).EQ.'KXORD') read (S2,*) P%KXORD
        if (trim(S1).EQ.'KYORD') read (S2,*) P%KYORD
        if (trim(S1).EQ.'JPZ_FILENAME') read (S2,*) P%JPZ_FILENAME
        if (trim(S1).EQ.'RHO_FILENAME') read (S2,*) P%RHO_FILENAME
        if (trim(S1).EQ.'ASM_FILENAME') read (S2,*) P%ASM_FILENAME
        if (trim(S1).EQ.'LCW_FILENAME') read (S2,*) P%LCW_FILENAME
        if (trim(S1).EQ.'OUT_FORMAT') read (S2,*) P%OUT_FORMAT
        if (trim(S1).EQ.'NOUT') read (S2,*) P%NOUT
        if (trim(S1).EQ.'LTC_CONST_A') read (S2,*) P%LTC_CONST_A
        if (trim(S1).EQ.'LTC_CONST_B') read (S2,*) P%LTC_CONST_B
        if (trim(S1).EQ.'N_FOLDING_A') read (S2,*) P%N_FOLDING_A
        if (trim(S1).EQ.'N_FOLDING_B') read (S2,*) P%N_FOLDING_B
      end do
      close (11)
      end subroutine

      subroutine INPUT_PARAMETER(P)
      use VAR_REC2D
      implicit none
      type(REC2D_PARAMETER) :: P

      print *, 'symmetry 2/4/6 ?                          '
      read (*,*) P%NSYMM
      print *, 'Number of channels for calc (256-2048) ?  '
      read (*,*) P%ND
      print *, 'Number of Compton profiles ?              '
      read (*,*) P%NP
      print *, 'BIN size (a.u./channel) ?                 '
      read (*,*) P%DCH
      print *, 'Momentum resolution (a.u.) ?              '
      read (*,*) P%DPZ
      print *, 'Smoothing parameter for Compton profiles  '
      print *, ' >0:c-cpline factor, 0:off, -1:auto     ? '
      read (*,*) P%SMPAR
      print *, 'Low pass filter type for B(qr) profiles'
      print *, ' 0:off, 1:Parzon(triangle) 2:Hamming(Cos)'
      print *, ' 3:Welch(parabola) (fixed) 4:trapezoid'
      print *, ' 5:Gaussian 6:Tanaka (depends on dPz)     '
      read (*,*) P%FLPAR
      print *, 'Order of B-spline for inperpolation B(qr) q-axis? '
      read (*,*) P%KXORD
      print *, 'Order of B-spline for inperpolation B(qr) r-axis? '
      read (*,*) P%KYORD
      print *, 'Zoom ratio for output data ?              '
      read (*,*) P%ZMPAR
      print *, 'Number of BINs of output files            '
      read (*,*) P%NOUT
      print *, 'OUTPUT FileName (rho(X,Y))                '
      read (*,*) P%RHO_FILENAME
      print *, 'OUTPUT FileName (rho(X,Y) asymmetry       '
      read (*,*) P%ASM_FILENAME
      print *, 'OUTPUT FileName (n_k(X,Y))                '
      read (*,*) P%LCW_FILENAME
      print *, 'Output file format 0:table / 1:coulmn ?   '
      read (*,*) P%OUT_FORMAT
      print *, 'Lattice constant a (Angstrom)             '
      read (*,*) P%LTC_CONST_A
      print *, 'Lattice constant b (Angstrom)             '
      read (*,*) P%LTC_CONST_B
      print *, 'Number of BINs of output files            '
      read (*,*) P%NOUT
      print *, 'Number of holding times (a-axis)          '
      read (*,*) P%N_FOLDING_A
      print *, 'Number of holding times (b-axis)          '
      read (*,*) P%N_FOLDING_B
      end subroutine

      subroutine CHECK_PARAMETER(P)
      use VAR_REC2D
      implicit none
      type(REC2D_PARAMETER) :: P

      real(kind(1d0)), parameter      :: a0 = 0.52917721092 ! 1 a.u. (length)
      real(kind(1d0)), parameter      :: pi = 3.14159265358979324d0
      real(kind(1d0))      :: A, B

      if (P%LTC_CONST_B == -1) then
        P%LTC_CONST_B = P%LTC_CONST_A
      end if

      if (P%LTC_CONST_A < 0) then
        P%N_FOLDING_A = -1
        P%N_FOLDING_B = -1
      else
        P%NLCW_A = floor((pi * a0 / P%LTC_CONST_A ) / (P%DCH / P%ZMPAR))+1  ![ch]
        P%NLCW_B = floor((pi * a0 / P%LTC_CONST_B ) / (P%DCH / P%ZMPAR))+1
	  
        P%STEP_FOLDING_A = (2. * pi * a0 / P%LTC_CONST_A ) / (P%DCH / P%ZMPAR) ![ch]
        P%STEP_FOLDING_B = (2. * pi * a0 / P%LTC_CONST_B ) / (P%DCH / P%ZMPAR)
	  
        A = floor((P%ND - P%STEP_FOLDING_A/2. - 1.) / P%STEP_FOLDING_A)
        B = floor((P%ND - P%STEP_FOLDING_B/2. - 1.) / P%STEP_FOLDING_B)

        if ((A < P%N_FOLDING_A) .or. (P%N_FOLDING_A < 0)) then
            P%N_FOLDING_A = A
        end if
        if ((B < P%N_FOLDING_B) .or. (P%N_FOLDING_B < 0)) then
            P%N_FOLDING_B = B
        end if
      end if

      end subroutine
	  
      subroutine READ_COMPTON_PROFILE(P, TH, PZ, JPZ)
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: JPZ(P%ND, P%NP)
      real(kind(1d0)) :: TH(P%NP), PZ(P%ND)
      integer         :: I,J,K
      integer         :: I1,I2,I3
      character*256   :: S

      JPZ = 0.d0
      if (P%JPZ_FILENAME(1:1).EQ.'*') then ! single profile file
        do I=1, P%NP
          print *, ' Direction of Profile (Degree), FileName'
          read (*,*) TH(I), P%JPZ_FILENAME
          open (11, file = trim(P%JPZ_FILENAME), READONLY, status='OLD')
            do J=1, P%ND
              read (11, *, iostat=K) JPZ(J,I)
              if (K.LT.0) exit
            end do
            P%NORG = J
          close (11)
        end do
      else  ! multi profile file / table format
        open (11, file = trim(P%JPZ_FILENAME), READONLY, status='OLD')
          read (11, *, iostat=K) TH(1:P%NP)
          do J=1, P%ND
            read (11, *, iostat=K) PZ(J), JPZ(J,1:P%NP)
            if (K.LT.0) exit
          end do
          P%NORG = J
        close (11)
      end if
      do J=1, P%ND
        PZ(J) = float(J-1) * P%DCH
      end do
      end subroutine
      
      subroutine SMOOTH_JPZ(P, JPZ)
! Smooth Compton profiles J(pz) using C-SPLINE_FITTING
!   input  : P, JP(,)
!   output : JP(,)
      use CSSMH_INT       ! c-spline smoothing
      use CSVAL_INT       ! p761
      use CSSCV_INT       ! 
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: JPZ(P%ND, P%NP)
      
      integer         :: I,J
      real(kind(1d0)) :: X

      integer         :: IEQUAL
      real(kind(1d0)) :: XDATA(P%ND*2)
      real(kind(1d0)) :: FDATA(P%ND*2)
      real(kind(1d0)) :: BREAK(P%ND*2)
      real(kind(1d0)) :: CSCOEF(4,P%ND*2)
      real(kind(1d0)) :: WEIGHT(P%ND*2)

      if (P%SMPAR .EQ. 0.) then ! if 0 then skip smoothing
        return
      end if 
      do J=1, P%ND*2
        XDATA(J) = float(J)
      end do
      do I=1, P%NP
        do J=1, P%ND ! expand profiles from + to +- region 
          FDATA(J + P%ND) = JPZ(J,I)
          FDATA(P%ND-J+2) = JPZ(J,I)
        end do
        FDATA(1) = FDATA(2) ! fill end channel to next channel value
        if (P%SMPAR > 0.) then ! user input value
          call CSSMH (XDATA, FDATA, P%SMPAR, BREAK, CSCOEF)
        else ! default value (automatic)
          IEQUAL = 1
          call CSSCV (XDATA, FDATA, IEQUAL, BREAK, CSCOEF)
        end if
        do J=1, P%NORG ! collapse profiles from -+ to + region
          X = float(J + P%ND)
          JPZ(J,I) = CSVAL(X, BREAK, CSCOEF) 
        end do
        do J=P%NORG+1, P%ND ! fill 0 over data channel of original file
          JPZ(J,I) = 0d0
        end do
      end do
      end subroutine
      
      subroutine FT_JPZ(P, JPZ, BRQ)
! FFT (forward transform) J(pz) -> B(rq)
!   input  : P, JPZ(:)
!   output : BRQ(:)
      use FCOSI_INT       ! COS Transform p.1028
      use F2OST_INT
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: JPZ(P%ND, P%NP)
      real(kind(1d0)) :: BRQ(P%ND, P%NP)

      integer :: I
      real(kind(1d0)) :: WFCOS(3*P%ND+15)        ! 3 * NDATA + 15

      call FCOSI (P%ND, WFCOS)
      do I=1, P%NP
        call F2OST(P%ND, JPZ(:,I), BRQ(:,I), WFCOS)
      end do
      BRQ(1,:) = sum(BRQ(1,:)) / float(P%NP)  ! B(r=0) should be same value
      end subroutine
      
      subroutine FILTER_BRQ(P, BRQ)
! Filter (Low pass filter to reduce high frequency noise for B(pz))
!   input  : P, BRQ(:)
!   output : BRQ(:)
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: BRQ(P%ND, P%NP)    ! B(rq) profile

      integer :: I, J, D
      real(kind(1d0)) :: FL(P%ND)  ! Filter for FFT
      real(kind(1d0)) :: PI
      PI      = 4d0*datan(1d0)

      if (P%FLPAR .EQ. 0.) then ! if 0 then skip filtering
        return
      endif
      FL = 1d0
      do I=1, P%ND
        select case (P%FLPAR)
          case (1) ! Parzen : Triangle window
            FL(I) = float(P%ND-(I-1))/float(P%ND)
          case (2) ! Hamming : Cosine window
            FL(I) = (1.d0+ cos((float(I-1)*PI)/float(P%ND)) ) / 2.d0
          case (3) ! Welch : Parabolla
            FL(I) = 1. - (float(I-1) / float(P%ND) ) ** 2.
          case (4) !  : trapezoid
            D = float(P%ND) / P%DPZ * P%DCH
            if ((I-1) > int(D)) FL(I) = (float(P%ND) - float(I-1)) / (float(P%ND)-D)
          case (5) !  : Gaussian
            D = float(P%ND) / P%DPZ * P%DCH / sqrt(2. * log(2.))
            FL(I) = exp( - float(I-1) * float(I-1) / 2.d0 / D / D)
          case (6) ! Tanaka/Nagao
            D = float(P%ND) / P%DPZ * P%DCH * log(2.d0) / PI 
            if (I > int(D)) FL(I) = 0.5d0 ** ((float(I-1)-D)/ D)
        end select
      end do
      open (22, FILE = 'FILTER.tmp')
        do I=1, P%ND
          write (22, '(E13.5)') FL(I)
        end do
      close (22)
      do J=1, P%NP
        BRQ(:,J) = BRQ(:,J) * FL
      end do
      end subroutine

      subroutine TRANS_BP2BR(P, TH, PZ, BRQ, B_R)
! 2D interpolation / B(rq) -> B(R)
!   input  : P, TH(), pz(), BRQ(,)
!   output : pz(), B_R(,)
      use BS2IN_INT       ! B-Spline p.616
      use BSNAK_INT
      use BS2VL_INT
      use BS2DR_INT
      use VAR_REC2D
!     use IMSL_LIBRARIES
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: BRQ(P%ND, P%NP)
      real(kind(1d0)) :: B_R(P%ND, P%ND)
      real(kind(1d0)) :: TH(P%NP), PZ(P%ND)

      integer         :: NT            ! NTHET = NPROF*3-2
      real(kind(1d0)), allocatable    :: TDATA(:)      ! (NT)
      real(kind(1d0)), allocatable    :: RDATA(:)      ! (ND)
      real(kind(1d0)), allocatable    :: FDATA(:,:)    ! (ND, NT)             ! 
      real(kind(1d0)), allocatable    :: XKNOT(:)      ! (ND+KXORD)
      real(kind(1d0)), allocatable    :: YKNOT(:)      ! (NP+KYORD)
      real(kind(1d0)), allocatable    :: BSCOEF(:,:)   ! (NP, ND)

      integer         :: I, J
      real(kind(1d0))                 :: THETA,RADIUS  ! polar coodinate
      real(kind(1d0))                 :: X,Y           ! rectangular coodinate
      real(kind(1d0))                 :: PI
      PI      = 4d0*datan(1d0)

      NT = P%NP * 3 - 2 ! number of profiles for interpolation
      allocate (FDATA(P%ND, NT),      TDATA(NT)        , RDATA(P%ND) )
      allocate (XKNOT(P%ND +P%KXORD), YKNOT(NT+P%KYORD), BSCOEF(NT,P%ND))
      do I=1, P%ND
        RDATA(I) = float (I-1)
      end do

      ! copy profiles to symmetric position
      do I=1, NT	! J(original) to I theta
        if (I < P%NP) then ! - region
          J = P%NP - I + 1
          TDATA(I) = - TH(J) * PI / 180.d0
        else if (I < (P%NP) * 2) then ! real position
          J = I - P%NP + 1
          TDATA(I) = TH(J) * PI / 180.d0
        else ! + region
          J = 3 * P%NP - I - 1
          if (P%NSYMM == 2) then
            TDATA(I) = (180.d0 - TH(J)) * PI / 180.d0
          else if (P%NSYMM == 4) then
            TDATA(I) = (90.d0 - TH(J)) * PI / 180.d0
          else if (P%NSYMM == 6) then
            TDATA(I) = (60.d0 - TH(J)) * PI / 180.d0
          end if
        end if
!        print *, I, J, TDATA(I)
        FDATA(:,I) = BRQ(:,J)                               ! 
      end do

      ! prepare for interpolation
      ! BSNAK : Knot sequence given interpolation data
      call BSNAK (P%ND , RDATA, P%KXORD, XKNOT) ! X
      call BSNAK (NT   , TDATA, P%KYORD, YKNOT) ! Y
      call BS2IN (RDATA, TDATA, FDATA, P%KXORD, P%KYORD, XKNOT, YKNOT, BSCOEF)

      ! interpolation  ! To chage Bin size, tume R factor
      do I=1, P%ND      ! B(X,Y), calc (theta 0 : PI/2)
        if (mod(I, (P%ND/10)) == 0) print '("*",$)' ! disp *
        do J=1, P%ND
          X = float(I-1)
          Y = float(J-1)
          RADIUS = sqrt(X**2. + Y**2. ) * P%ZMPAR
          if (X == 0) then
            THETA = PI/2.
          else
            THETA = atan( Y/X )
          end if
          if (P%NSYMM == 4) then
            if (THETA > PI/4.) then
              THETA = PI/2. - THETA
            end if
          end if
          if (P%NSYMM == 6) then
            if (THETA > PI/3.) then
              THETA = THETA - PI/3.
            else if (THETA > PI/6.) then
              THETA = PI/3. - THETA
            end if
          end if
          if ( RADIUS > float(P%ND-1) ) then
            B_R(I,J) = 0.d0
          else
            B_R(I,J) = BS2VL(RADIUS, THETA, P%KXORD, P%KYORD, XKNOT, YKNOT, P%ND, NT, BSCOEF)
!           B_R(I,J) = QD2VL(RADIUS, THETA, RDATA, TDATA, FDATA)
          end if
        end do
      end do

      do I=1, P%ND
        PZ(I) = float(I-1)*P%DCH/P%ZMPAR
      end do

      deallocate (FDATA, TDATA, RDATA)
      deallocate (XKNOT, YKNOT, BSCOEF)
      end subroutine
      
      subroutine IFT_BR_2D(P, B_R, R_R)
! iFFT 2D   B(R) -> rho(R)
!   input  : P, B_R(,)
!   output : R_R(,)
      use FCOSI_INT       !   FCOST ()		IMSL manual p.1028
      use F2OST_INT
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: B_R(P%ND, P%ND)
      real(kind(1d0)) :: R_R(P%ND, P%ND)
      real(kind(1d0)) :: T_XY(P%ND, P%ND)       ! (NDATA, NDATA) temp for 2dFFT
      real(kind(1d0)) :: WFCOS(3*P%ND + 15)            ! 3 * NDATA + 15                   Cos変換用変数
      integer :: I
      real(kind(1d0)) :: X
      
      R_R(:,:)=0
      call FCOSI (P%ND, WFCOS)
      do I=1, P%ND
        if (mod(I, (P%ND/5)) == 0) print '("*",$)'
        call F2OST(P%ND, B_R(I,:), T_XY(I,:), WFCOS)
      end do
      do I=1, P%ND
        if (mod(I, (P%ND/5)) == 0) print '("*",$)'
        call F2OST(P%ND, T_XY(:,I), R_R(:,I), WFCOS)
      end do
      X = (2. / P%ZMPAR * float(P%ND))**2. * P%DCH
      R_R = R_R / X
      end subroutine

      subroutine ASYMMETRY(P, R_R, ARR)
      use CSINT_INT       ! c-spline smoothing
      use CSVAL_INT       ! p761
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: R_R(P%ND, P%ND)
      real(kind(1d0)) :: ARR(P%ND, P%ND)
      
      integer         :: I,J,K
      real(kind(1d0)) :: R, DR, X, Y

      real(kind(1d0)) :: XDT(P%ND)
      real(kind(1d0)) :: FDT(P%ND)
      real(kind(1d0)) :: BREAK(P%ND)
      real(kind(1d0)) :: CSCOEF(4,P%ND)
      real(kind(1d0)) :: NDT(P%ND)

      NDT(:)=0
      FDT(:)=0
      do I=1, P%ND ! make average data (1d)
        do J=1, P%ND
          R = sqrt(float((I-1)*(I-1) + (J-1)*(J-1))) ! radius
          DR = R - floor(R) 
          K = floor(R)+1 ! ch
          if ((K>0).and.(K<P%ND))  then
            FDT(K) = FDT(K) + R_R(I,J) * (1.-DR)
            NDT(K) = NDT(K) + (1.-DR)
            FDT(K + 1) = FDT(K + 1) + R_R(I,J) * DR
            NDT(K + 1) = NDT(K + 1) + DR
          end if
        end do
      end do

      do I=1, P%ND
        XDT(I) = float(I-1)
        if (NDT(I)>0) then
          FDT(I) = FDT(I) / NDT(I)
        end if
      end do

	   ! subtract average data (2d)
      ARR(:,:)=0
      call CSINT (XDT, FDT, BREAK, CSCOEF)

      do I=1, P%ND      ! R(X,Y), calc (theta 0 : PI/2)
        if (mod(I, (P%ND/10)) == 0) print '("*",$)'
        do J=1, P%ND
          R = sqrt(float((I-1)*(I-1) + (J-1)*(J-1)))
          if ( R > float(P%ND-1) ) then
            ARR(I,J) = 0.d0
          else
            ARR(I,J) = R_R(I,J) - CSVAL(R, BREAK, CSCOEF) 
          end if
        end do
      end do

      end subroutine
	  
      subroutine LCW_FOLDING(P, R_R, N_K)
!LCW
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: R_R(P%ND,        P%ND)
      real(kind(1d0)) :: N_K(-P%NLCW_A:P%NLCW_A, -P%NLCW_B:P%NLCW_B)
      real(kind(1d0)) :: R(-P%ND+1:P%ND, -P%ND+1:P%ND)

      integer         :: I,J
      integer         :: IA, IB, JA, JB, KA, KB
      integer         :: SA, SB
      real(kind(1d0)) :: STEP_A, STEP_B
      real(kind(1d0)) :: DA, DB
      real(kind(1d0)) :: F1, F2, F3, F4
      real(kind(1d0)) :: X, S

      R(:,:) = 0.
      N_K(:,:) = 0.
      X = 0.
      do I=-P%ND+1, P%ND-1 ! ４象限 に 展開
        do J=-P%ND+1, P%ND-1
          R(J, I) = R_R(abs(J)+1, abs(I)+1) 
          X = X + R(J, I)
        end do
      end do

      print *, "N foldings ", P%N_FOLDING_A, P%N_FOLDING_B
      do IA = -P%N_FOLDING_A, P%N_FOLDING_A  ! 重ねる
        SA = floor(P%STEP_FOLDING_A * float(IA)) !境目のマイナス側の ch
        DA = (P%STEP_FOLDING_A * float(IA)) - SA
        do IB = -P%N_FOLDING_B, P%N_FOLDING_B
          SB = floor(P%STEP_FOLDING_B * float(IB))
          DB = (P%STEP_FOLDING_B * float(IB)) - SB

          F1 = (1. - DA)*(1. - DB)
          F2 = DA * (1. - DB)
          F3 = DB * (1. - DA)
          F4 = DA * DB
          
          do JA = -P%NLCW_A, P%NLCW_A
            do JB = -P%NLCW_B, P%NLCW_B
               KA = JA + SA
               KB = JB + SB
				S = R(KA, KB)*F1 + R(KA+1, KB)*F2 + R(KA, KB+1)*F3 + R(KA+1, KB+1)*F4
               N_K(JA, JB) = N_K(JA, JB) + S
            end do
          end do
        end do
      end do
      end subroutine
  
      subroutine REMAKE_JPZ(P, R_R, TH, JPZ)
! integrate rho(R) to make J(PZ)
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: R_R(P%ND, P%ND)
      real(kind(1d0)) :: JPZ(P%ND, P%NP)
      real(kind(1d0)) :: TH(P%NP), PZ(P%ND)
      integer :: I, J, K, IX, IY
      real(kind(1d0)) :: X, Y, Z
      real :: cosT, sinT
      real(kind(1d0))                 :: PI
      PI      = 4d0*datan(1d0)

      do I=1, P%NP
        print '("*",$)'
        cosT = cos(-TH(I)*PI/180.)
        sinT = sin(-TH(I)*PI/180.)
        do J=1, P%ND
          Z = 0.
          do K=-P%ND, P%ND
            IX = abs(anint(float(J-1)*cosT - float(K-1)*sinT))
            IY = abs(anint(float(J-1)*sinT + float(K-1)*cosT))
            if ((IX < P%ND-1) .and. (IY < P%ND-1)) then
            Z = Z + R_R(IX+1, IY+1)
            end if
          end do
          JPZ(J,I) = Z
        end do
      end do
      JPZ = JPZ * P%DCH / P%ZMPAR
      end subroutine

      subroutine REMAKE_JPZ_FINE(P, R_R, TH, JPZ)
! integrate rho(R) to make J(PZ)
      use BS2IN_INT       ! B-Spline p.616
      use BSNAK_INT
      use BS2VL_INT
      use BS2DR_INT
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      real(kind(1d0)) :: R_R(P%ND, P%ND)
      real(kind(1d0)) :: JPZ(P%ND, P%NP)
      real(kind(1d0)) :: TH(P%NP), PZ(P%ND)
      integer :: I, J, K
      real(kind(1d0)) :: X, Y, Z
      real :: cosT, sinT
      real(kind(1d0))                 :: PI
      real(kind(1d0)), allocatable    :: XDATA(:)      ! (ND)
      real(kind(1d0)), allocatable    :: YDATA(:)      ! (ND)
      real(kind(1d0)), allocatable    :: XKNOT(:)      ! (ND+KXORD)
      real(kind(1d0)), allocatable    :: YKNOT(:)      ! (ND+KYORD)
      real(kind(1d0)), allocatable    :: BSCOEF(:,:)   ! (ND, ND)

      PI      = 4d0*datan(1d0)
      allocate (XDATA(P%ND)    , YDATA(P%ND) )
      allocate (XKNOT(P%ND+P%KXORD),YKNOT(P%ND+P%KYORD),BSCOEF(P%ND,P%ND))
      do I=1, P%ND
        XDATA(I) = float (I-1)
        YDATA(I) = float (I-1)
      end do

      print '("*",$)'
      call BSNAK (P%ND , XDATA, P%KXORD, XKNOT) ! X
      call BSNAK (P%ND , YDATA, P%KYORD, YKNOT) ! Y
      call BS2IN (XDATA, YDATA, R_R, P%KXORD, P%KYORD, XKNOT, YKNOT, BSCOEF)

      do I=1, P%NP
        print '("*",$)'
        cosT = cos(-TH(I)*PI/180.)
        sinT = sin(-TH(I)*PI/180.)
        do J=1, P%ND
          Z = 0.
          do K=-P%ND, P%ND
            X = (float(J-1)*cosT - float(K-1)*sinT)
            Y = (float(J-1)*sinT + float(K-1)*cosT)
            X = abs(X)
            Y = abs(Y)
            if ((X < P%ND-1) .and. (Y < P%ND-1)) then
             Z = Z + BS2VL(X, Y, P%KXORD, P%KYORD, XKNOT, YKNOT, P%ND, P%ND, BSCOEF)
            end if
          end do
          JPZ(J,I) = Z
        end do
      end do
      JPZ = JPZ * P%DCH
      deallocate (XDATA, YDATA, XKNOT, YKNOT, BSCOEF)
      end subroutine
      
      subroutine WRITE_PROFILE(NX, NY, X, Y, D, FILENAME)
      implicit none

      integer         :: NX, NY
      real(kind(1d0)) :: D(NY,NX) ! Compton profile
      real(kind(1d0)) :: X(NX), Y(NY)
      integer         :: I
      character*256   :: FILENAME 
      
      open (22, FILE = trim(FILENAME))
        write (22, '(99E13.5)') X(1:NX)
        do I=1, NY
          write (22, '(99E13.5)') Y(I), D(I,1:NX)
        end do
      close (22)
      end subroutine
      
      subroutine WRITE_PROFILE_2D(NX, NOUT, X, D, FILENAME, F)
      implicit none

      integer         :: NX, NOUT
      real(kind(1d0)) :: X(NX)
      real(kind(1d0)) :: D(NX,NX)
      integer         :: I,J,F
      character*256   :: FILENAME 

      open (22, FILE = trim(FILENAME))
        if (F==1) then
          do I=1, NOUT
            do J=1, NOUT
              write (22, '(3E13.5)') X(I), X(J), D(I,J)
            end do
          end do
        else
          write (22, '(9999E13.5)') X(1:NOUT)
          do I=1, NOUT
            write (22, '(9999E13.5)') X(I), D(I,1:NOUT)
          end do
        end if
      close (22)
      end subroutine

      subroutine WRITE_N_K_2D(NA, NB, X, NX, D, FILENAME, F)
      implicit none

      integer         :: NX, NA, NB
      real(kind(1d0)) :: X(NX)
      real(kind(1d0)) :: D(-NA:NA,-NB:NB)
      integer         :: I,J,F
      character*256   :: FILENAME 

      open (22, FILE = trim(FILENAME))
        if (F==1) then
          do I=-NA, NA
            do J=-NB, NB
              write (22, '(9999E13.5)') sign(X(abs(I)+1), float(I)), sign(X(abs(J)+1), float(J)), D(I,J)
            end do
          end do
        else
!          do I=-NA, NA
!            write (22, '(9999E13.5)') -X(abs(NA:1)), 0 ,-X(abs(1:NA)
!          end do
          do I=-NB, NB
!            write (22, '(9999E13.5)') sign(X(abs(J)+1), D(I,-NA:NA)
            write (22, '(9999E13.5)') D(I,-NA:NA)
          end do
        end if
      close (22)
      end subroutine

      subroutine WRITE_PROFILE_BMP(N, NOUT, D, FILENAME)
      implicit none

      integer         :: N, NOUT, N1
      real(kind(1d0)) :: D(N,N) ! Compton profile
      character*256   :: FILENAME 
      integer             :: I,J
      real(kind(1d0))     :: DMAX, DMIN, DH
      integer(1)          :: B                        ! byte
      integer(2)          :: W                        ! word
      integer(4)          :: L                        ! double word
      integer H, H1

!BitMap
      H=1078 !header size
      N1=NOUT*2-1
      H1=mod((NOUT*2-1),4)
      open (22, FILE = (trim(FILENAME) // ".BMP"), form='binary') ! Bitmap 作成 固定の ヘッダに
        write(22) 'BM'
        L=(N1+H1)*N1+H; write(22) L ! file size
        L=0; write(22) L; L=H; write(22) L
        
        L=40; write(22) L
        L=N1; write(22) L; write(22) L ! Bitmap size
        W=1; write(22) W;  W=8; write(22) W ! 256 color bitmap
        L=0; write(22) L
        L=(N1+H1)*N1; write(22) L ! Image size (byte)
        L=0; write(22) L; write(22) L; L=256; write(22) L; L=0; write(22) L
        do I=0, 255 ! color palette
          B = I
          write (22) B; write (22) B; write (22) B; B=0; write (22) B
        end do
        
        DMAX = maxval(D)
        DMIN = minval(D)
        DH=255./(DMAX-DMIN)
        do J=-NOUT+1, NOUT-1
          do I=-NOUT+1, NOUT-1
            B = int((D(abs(J)+1, abs(I)+1)-DMIN)*DH)
            write (22) B
          end do
          if (H1==3) write (22) B ! padding
          if (H1==1) write (22) B,B,B
        end do
      close(22)
      end subroutine

      subroutine WRITE_N_K_BMP(NA, NB, D, FILENAME)
      implicit none

      integer         :: NX, NA, NB
      real(kind(1d0)) :: D(-NA:NA,-NB:NB)
      integer         :: I,J,F
      character*256   :: FILENAME 

      integer             :: IA, IB
      real(kind(1d0))     :: DMAX, DMIN, DH
      integer(1)          :: B                        ! byte
      integer(2)          :: W                        ! word
      integer(4)          :: L                        ! double word
      integer H, H1

!BitMap
      H=1078 !header size
      H1=mod((NB*2+1),4)

      open (22, FILE = (trim(FILENAME) // ".BMP"), form='binary') ! Bitmap 作成 固定の ヘッダに
        write(22) 'BM'
        L=((NA*2+1)+H1)*(NB*2+1)+H; write(22) L ! file size

        L=0; write(22) L; L=H; write(22) L
        L=40; write(22) L
        L=(NB*2+1); write(22) L ! Bitmap size
		L=(NA*2+1); write(22) L ! Bitmap size
        W=1; write(22) W;  W=8; write(22) W ! 256 color bitmap
        L=0; write(22) L
        L=((NA*2+1)+H1)*(NB*2+1); write(22) L ! Image size (byte)
        L=0; write(22) L; write(22) L; L=256; write(22) L; L=0; write(22) L
        do I=0, 255 ! color palette
          B = I
          write (22) B; write (22) B; write (22) B; B=0; write (22) B
        end do
        
        DMAX = maxval(D)
        DMIN = minval(D)
        DH=255./(DMAX-DMIN)
        do IA = -NA, NA
          do IB= -NB, NB
            B = int((D(IA, IB)-DMIN)*DH)
            write (22) B
          end do
          if (H1==3) write (22) B ! padding
          if (H1==1) write (22) B,B,B
        end do
      close(22)
      end subroutine

      subroutine WRITE_INIFILE(P, FILENAME)
      use VAR_REC2D
      implicit none

      type(REC2D_PARAMETER) :: P
      character*256   :: FILENAME
      character*256   :: S

      open (22, FILE = trim(FILENAME))
        write (S, *) P%NSYMM
        write (22, '(A)') 'NSYMM='//trim(adjustl(S))
        write (S, *) P%NP
        write (22, '(A)') 'NP='//trim(adjustl(S))
        write (S, *) P%ND
        write (22, '(A)') 'ND='//trim(adjustl(S))
        write (S, *) P%NOUT
        write (22, '(A)') 'NOUT='//trim(adjustl(S))
        write (S, *) P%DCH
        write (22, '(A)') 'DCH='//trim(adjustl(S))
        write (S, *) P%DPZ
        write (22, '(A)') 'DPZ='//trim(adjustl(S))
        write (S, *) P%SMPAR
        write (22, '(A)') 'SMPAR='//trim(adjustl(S))
        write (S, *) P%FLPAR
        write (22, '(A)') 'FLPAR='//trim(adjustl(S))
        write (S, *) P%ZMPAR
        write (22, '(A)') 'ZMPAR='//trim(adjustl(S))
        write (S, *) P%KXORD
        write (22, '(A)') 'KXORD='//trim(adjustl(S))
        write (S, *) P%KYORD
        write (22, '(A)') 'KYORD='//trim(adjustl(S))
        write (22, '(A)') 'JPZ_FILENAME='//trim(P%JPZ_FILENAME)
        write (22, '(A)') 'RHO_FILENAME='//trim(P%RHO_FILENAME)
        write (22, '(A)') 'LCW_FILENAME='//trim(P%LCW_FILENAME)
        write (S, *) P%LTC_CONST_A
        write (22, '(A)') 'LTC_CONST_A='//trim(adjustl(S))
        write (S, *) P%LTC_CONST_B
        write (22, '(A)') 'LTC_CONST_B='//trim(adjustl(S))
        write (S, *) P%OUT_FORMAT
        write (22, '(A)') 'OUT_FORMAT='//trim(adjustl(S))
        write (S, *) P%N_FOLDING_A
        write (22, '(A)') 'N_FOLDING_A='//trim(adjustl(S))
        write (S, *) P%N_FOLDING_B
        write (22, '(A)') 'N_FOLDING_B='//trim(adjustl(S))
      end subroutine
