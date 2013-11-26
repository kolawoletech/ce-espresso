!
! Copyright (C) 2013 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Nov 2013, D. Ceresoli: Z2 calculation according Soluyanov and Vanderbilt,
!                        PRB 83, 235401 (2011)
!
! This routine is largerly inspired from bp_c_phase.f90
!
!##############################################################################!
!#                                                                            #!
!#                                                                            #!
!#   This is the main one of a set of Fortran 90 files designed to compute    #!
!#   the Z2 invarinat without inversion symmetry.                             #!
!#                                                                            #!
!#   This routine is largerly inspired from bp_c_phase.f90                    #!
!#                                                                            #!
!#                                                                            #!
!#   BRIEF SUMMARY OF THE METHODOLOGY  (REWRITE)                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                         #!
!#   The standard procedure would be for the user to first perform a          #!
!#   self-consistent (sc) calculation to obtain a converged charge density.   #!
!#   With well-converged sc charge density, the user would then run one       #!
!#   or more non-self consistent (or "band structure") calculations,          #!
!#   using the same main code, but with a flag to ask for the polarization.   #!
!#   Each such run would calculate the projection of the polarization         #!
!#   onto one of the three primitive reciprocal lattice vectors. In           #!
!#   cases of high symmetry (e.g. a tetragonal ferroelectric phase), one      #!
!#   such run would suffice. In the general case of low symmetry, the         #!
!#   user would have to submit up to three jobs to compute the three          #!
!#   components of polarization, and would have to obtain the total           #!
!#   polarization "by hand" by summing these contributions.                   #!
!#                                                                            #!
!#   Accurate calculation of the electronic or "Berry-phase" polarization     #!
!#   requires overlaps between wavefunctions along fairly dense lines (or     #!
!#   "strings") in k-space in the direction of the primitive G-vector for     #!
!#   which one is calculating the projection of the polarization. The         #!
!#   code would use a higher-density k-mesh in this direction, and a          #!
!#   standard-density mesh in the two other directions. See below for         #!
!#   details.                                                                 #!
!#                                                                            #!
!#                                                                            #!
!#   FUNCTIONALITY/COMPATIBILITY                                              #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   * Spin-polarized systems supported.                                      #!
!#                                                                            #!
!#   * Ultrasoft and norm-conserving pseudopotentials supported.              #!
!#                                                                            #!
!#   * Calculation must be non collinear and including spin orbit             #!
!#                                                                            #!
!#                                                                            #!
!#   NEW INPUT PARAMETERS                                                     #!
!#   ~~~~~~~~~~~~~~~~~~~~                                                     #!
!#   * lcalc_z2 (.TRUE. or .FALSE.)                                           #!
!#     Tells PWSCF that a Berry phase calcultion is desired.                  #!
!#                                                                            #!
!#   * gdir (1, 2, or 3)                                                      #!
!#     Specifies the direction of the k-point strings in reciprocal space.    #!
!#     '1' refers to the first reciprocal lattice vector, '2' to the          #!
!#     second, and '3' to the third.                                          #!
!#                                                                            #!
!#   * nppstr (integer)                                                       #!
!#     Specifies the number of k-points to be calculated along each           #!
!#     symmetry-reduced string.                                               #!
!#                                                                            #!
!#                                                                            #!
!#   EXPLANATION OF K-POINT MESH  (REWRITE)                                   #!
!#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~                                              #!
!#   If gdir=1, the program takes the standard input specification of the     #!
!#   k-point mesh (nk1 x nk2 x nk3) and stops if the k-points in dimension    #!
!#   1 are not equally spaced or if its number is not equal to nppstr,        #!
!#   working with a mesh of dimensions (nppstr x nk2 x nk3).  That is, for    #!
!#   each point of the (nk2 x nk3) two-dimensional mesh, it works with a      #!
!#   string of nppstr k-points extending in the third direction.  Symmetry    #!
!#   will be used to reduce the number of strings (and assign them weights)   #!
!#   if possible.  Of course, if gdir=2 or 3, the variables nk2 or nk3 will   #!
!#   be overridden instead, and the strings constructed in those              #!
!#   directions, respectively.                                                #!
!#                                                                            #!
!#                                                                            #!
!#   BIBLIOGRAPHY                                                             #!
!#   ~~~~~~~~~~~~                                                             #!
!#                                                                            #1
!#   [1] A. A. Soluyanov and D. Vanderbilt, "Computing topological            #!
!#       invariants without inversion symmetry",                              #!
!#       Phys Rev B 83, 235401 (1993).                                        #!
!#   [2] http://www.physics.rutgers.edu/z2pack/                               #!
!#                                                                            #!
!#                                                                            #!
!##############################################################################!


!======================================================================
SUBROUTINE c_phase_z2
!======================================================================
   !
   !   Geometric phase calculation along a strip of nppstr k-points
   !   averaged over a 1D grid of nkort k-points ortogonal to nppstr 
   !  --- Make use of the module with common information ---
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : iunwfc, nwordwfc
   USE buffers,              ONLY : get_buffer
   USE cell_base,            ONLY : at, tpiba, omega, tpiba2
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm, g, gcutm, ngm_g, ig_l2g
   USE fft_base,             ONLY : dfftp
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, lmaxq, nbetam, nh, nhm
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nks, xk, wk
   USE wvfct,                ONLY : npwx, nbnd, ecutwfc
   USE wavefunctions_module, ONLY : evc
   USE bp,                   ONLY : gdir, nppstr, mapgm_global
   USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                    deallocate_bec_type
   USE noncollin_module,     ONLY : noncolin, npol, nspin_lsda
   USE spin_orb,             ONLY : lspinorb
   USE mp_bands,             ONLY : intra_bgrp_comm, nproc_bgrp
   USE mp,                   ONLY : mp_sum

!  --- Avoid implicit definitions ---
   IMPLICIT NONE

!  --- Internal definitions ---
   INTEGER :: i
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: info
   INTEGER :: is
   INTEGER :: istring
   INTEGER :: iv
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: jv
   INTEGER :: kindex
   INTEGER :: kort
   INTEGER :: kpar
   INTEGER :: kpoint
   INTEGER :: kstart
   INTEGER :: mb
   INTEGER :: mk1
   INTEGER :: mk2
   INTEGER :: mk3
   INTEGER , ALLOCATABLE :: ln(:,:,:)
   INTEGER :: n1
   INTEGER :: n2
   INTEGER :: n3
   INTEGER :: na
   INTEGER :: nb
   INTEGER :: ng
   INTEGER :: nhjkb
   INTEGER :: nhjkbm
   INTEGER :: nkbtona(nkb)
   INTEGER :: nkbtonh(nkb)
   INTEGER :: nkort
   INTEGER :: np
   INTEGER :: npw1
   INTEGER :: npw0
   INTEGER :: nstring
   INTEGER :: nbnd_occ
   INTEGER :: nt
   INTEGER, ALLOCATABLE :: map_g(:)
   LOGICAL :: l_para
   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for occupied/empty states
   REAL(DP) :: dk(3)
   REAL(DP) :: dkmod
   REAL(DP) :: el_loc
   REAL(DP), parameter :: eps = 1d-6
   REAL(DP) :: fac
   REAL(DP) :: g2kin_bp(npwx)
   REAL(DP) :: gpar(3)
   REAL(DP) :: gtr(3)
   REAL(DP) :: gvec
   REAL(DP), ALLOCATABLE :: loc_k(:)
   REAL(DP), ALLOCATABLE :: phik(:)
   REAL(DP) :: phik_ave
   REAL(DP) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(DP) :: weight
   REAL(DP) :: phidw
   REAL(DP) :: phiup
   REAL(DP), ALLOCATABLE :: wstring(:)
   REAL(DP) :: ylm_dk(lmaxq*lmaxq)
   COMPLEX(DP), ALLOCATABLE :: aux(:)
   COMPLEX(DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(DP), ALLOCATABLE :: aux0(:)
   TYPE (bec_type) :: becp0
   TYPE (bec_type) :: becp_bp
   COMPLEX(DP) :: cave
   COMPLEX(DP) , ALLOCATABLE :: cphik(:)
   COMPLEX(DP) :: dtheta
   COMPLEX(DP) :: mat(nbnd,nbnd)
   COMPLEX(DP) :: pref
   COMPLEX(DP), ALLOCATABLE :: psi(:,:)
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)
   COMPLEX(DP) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(DP) :: struc(nat)
   COMPLEX(DP) :: theta0
   COMPLEX(DP), external :: zdotc

!  -------------------------------------------------------------------------   !
!                               Z2 variables
!  -------------------------------------------------------------------------   !
   real(dp), parameter :: m_threshold = 0.7d0
   complex(dp), allocatable :: UU(:,:), VT(:,:), work(:), lambda(:,:), eig(:)
   real(dp), allocatable :: SV(:), rwork(:), zz(:), gaps(:)
   integer, allocatable :: ind(:)
   real(dp) :: point, max_gap
   integer :: lwork, kk

!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !
   ALLOCATE (psi(npwx*npol,nbnd))
   ALLOCATE (aux(ngm*npol))
   ALLOCATE (aux0(ngm*npol))
   IF (okvan) THEN
      CALL allocate_bec_type ( nkb, nbnd, becp0 )
      CALL allocate_bec_type ( nkb, nbnd, becp_bp )
      IF (lspinorb) ALLOCATE(q_dk_so(nhm,nhm,4,ntyp))
   END IF

   l_para= (nproc_bgrp > 1 .AND. gdir /= 3)
   IF (l_para) THEN
      ALLOCATE ( aux_g(ngm_g*npol) )
   ELSE
      ALLOCATE ( map_g(ngm) )
   ENDIF

!  --- Write header ---
   WRITE( stdout,"(/,/,/,15X,50('='))")
   WRITE( stdout,"(27X,'Z2 INVARIANT CALCULATION')")
   WRITE( stdout,"(25X,'!!! NOT THOROUGHLY TESTED !!!')")
   WRITE( stdout,"(15X,50('='),/)")
   if (.not. lspinorb) call errore('z2_c_phase', 'Z2 needs spin orbit', 1)
   if (nspin_lsda /= 1) call errore('z2_c_phase', 'internal error: nspin_lsda=', nspin_lsda)

!  --- Recalculate FFT correspondence (see ggen.f90) ---
   ALLOCATE (ln (-dfftp%nr1:dfftp%nr1, -dfftp%nr2:dfftp%nr2, -dfftp%nr3:dfftp%nr3) )
   DO ng=1,ngm
      mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   END DO

   if(okvan) then
!  --- Initialize arrays ---
      jkb_bp=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na).eq.nt) THEN
               DO i=1, nh(nt)
                  jkb_bp=jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i
               END DO
            END IF
         END DO
      END DO
   endif
!  --- Get the number of strings ---
   nstring=nks/nppstr
   nkort=nstring/nspin_lsda

!  --- Allocate memory for arrays ---
   ALLOCATE(phik(nstring))
   ALLOCATE(loc_k(nstring))
   ALLOCATE(cphik(nstring))
   ALLOCATE(wstring(nstring))

!  -------------------------------------------------------------------------   !
!           electronic polarization: set values for k-points strings           !
!  -------------------------------------------------------------------------   !

!  --- Find vector along strings ---
   gpar(1)=xk(1,nppstr)-xk(1,1)
   gpar(2)=xk(2,nppstr)-xk(2,1)
   gpar(3)=xk(3,nppstr)-xk(3,1)
   gvec=dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba

!  --- Find vector between consecutive points in strings ---
   dk(1)=xk(1,2)-xk(1,1)
   dk(2)=xk(2,2)-xk(2,1) 
   dk(3)=xk(3,2)-xk(3,1)
   dkmod=SQRT(dk(1)**2+dk(2)**2+dk(3)**2)*tpiba
   IF (ABS(dkmod-gvec/(nppstr-1)) > eps) & 
     CALL errore('c_phase_z2','Wrong k-strings?',1)

!  --- Check that k-points form strings ---
   DO i=1,nspin_lsda*nkort
      DO j=2,nppstr
         kindex=j+(i-1)*nppstr
         IF (ABS(xk(1,kindex)-xk(1,kindex-1)-dk(1)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(xk(2,kindex)-xk(2,kindex-1)-dk(2)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(xk(3,kindex)-xk(3,kindex-1)-dk(3)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings?',1)
         IF (ABS(wk(kindex)-wk(kindex-1)) > eps) &
            CALL errore('c_phase_z2','Wrong k-strings weights?',1)
      END DO
   END DO

   ! TODO: check k-points are on the edge of BZ

!  -------------------------------------------------------------------------   !
!                   electronic polarization: weight strings                    !
!  -------------------------------------------------------------------------   !

!  --- Calculate string weights, normalizing to 1 (no spin or noncollinear)
!       or 1+1 (spin) ---
   DO is=1,nspin_lsda
      weight=0.0_dp
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wk(nppstr*istring)
         weight=weight+wstring(istring)
      END DO
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wstring(istring)/weight
      END DO
   END DO

!  -------------------------------------------------------------------------   !
!                  electronic polarization: structure factor                   !
!  -------------------------------------------------------------------------   !

!  --- Calculate structure factor e^{-i dk*R} ---
   DO na=1,nat
      fac=(dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
      struc(na)=CMPLX(cos(fac),-sin(fac),kind=DP)
   END DO

!  -------------------------------------------------------------------------   !
!                     electronic polarization: form factor                     !
!  -------------------------------------------------------------------------   !
   if(okvan) then
!  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
      CALL calc_btq(dkmod,qrad_dk,0)

!  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
      dkmod=dk(1)**2+dk(2)**2+dk(3)**2
      CALL ylmr2(lmaxq*lmaxq, 1, dk, dkmod, ylm_dk)
!  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
      q_dk = (0.d0, 0.d0)
      DO np =1, ntyp
         if( upf(np)%tvanp ) then
            DO iv = 1, nh(np)
               DO jv = iv, nh(np)
                  call qvan3(iv,jv,np,pref,ylm_dk,qrad_dk)
                  q_dk(iv,jv,np) = omega*pref
                  q_dk(jv,iv,np) = omega*pref
               ENDDO
            ENDDO
         endif
      ENDDO
      IF (lspinorb) CALL transform_qq_so(q_dk,q_dk_so)
   endif

!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !

   el_loc=0.d0
   kpoint=0

   allocate (SV(nbnd), UU(nbnd,nbnd), VT(nbnd,nbnd), lambda(nbnd,nbnd), eig(nbnd), zz(nbnd))
   allocate (ind(nbnd), gaps(nbnd))
   allocate (l_cal(nbnd))  ! l_cal(n) = .true./.false. if n-th state is occupied/empty
   nbnd_occ = nbnd
   l_cal(1:nbnd) = .true.
   call weights()

!  --- Start loop over spin ---
   DO is=1,nspin_lsda

!     --- Start loop over orthogonal k-points ---
      DO kort=1,nkort

!        --- Index for this string ---
         istring=kort+(is-1)*nkort

   lambda = (0.d0, 0.d0)
   do nb = 1, nbnd
       lambda(nb,nb) = (1.d0, 0.d0)
   enddo

!        --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr

!           --- Set index of k-point ---
            kpoint = kpoint + 1

!           --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1) THEN

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
               CALL gk_sort(xk(1,kpoint-1),ngm,g,ecutwfc/tpiba2, &
                            npw0,igk0,g2kin_bp) 
               CALL get_buffer (psi,nwordwfc,iunwfc,kpoint-1)
               if (okvan) then
                  CALL init_us_2 (npw0,igk0,xk(1,kpoint-1),vkb)
                  CALL calbec (npw0, vkb, psi, becp0)
               endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= nppstr) THEN
                  CALL gk_sort(xk(1,kpoint),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)        
                  CALL get_buffer(evc,nwordwfc,iunwfc,kpoint)
                  if (okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kpoint),vkb)
                     CALL calbec (npw1, vkb, evc, becp_bp)
                  endif
               ELSE
                  kstart = kpoint-nppstr+1
                  CALL gk_sort(xk(1,kstart),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)  
                  CALL get_buffer(evc,nwordwfc,iunwfc,kstart)
                  if (okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kstart),vkb)
                     CALL calbec(npw1, vkb, evc, becp_bp)
                  endif
               ENDIF

               IF (kpar == nppstr .AND. .NOT. l_para) THEN
                  map_g(:) = 0
                  DO ig=1,npw1
!                          --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
!                          --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---

                     gtr(1)=g(1,igk1(ig)) - gpar(1)
                     gtr(2)=g(2,igk1(ig)) - gpar(2)
                     gtr(3)=g(3,igk1(ig)) - gpar(3)
!                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
!                          --- and the position ng in the ngm array ---

                     IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                        n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                             +gtr(3)*at(3,1))
                        n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                             +gtr(3)*at(3,2))
                        n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                             +gtr(3)*at(3,3))
                        ng=ln(n1,n2,n3)

                        IF ( (ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                             (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                             (ABS(g(3,ng)-gtr(3)) > eps) ) THEN
                           WRITE(6,*) ' error: translated G=', &
                                gtr(1),gtr(2),gtr(3), &
                                &     ' with crystal coordinates',n1,n2,n3, &
                                &     ' corresponds to ng=',ng,' but G(ng)=', &
                                &     g(1,ng),g(2,ng),g(3,ng)
                           WRITE(6,*) ' probably because G_par is NOT', &
                                &    ' a reciprocal lattice vector '
                           WRITE(6,*) ' Possible choices as smallest ', &
                                ' G_par:'
                           DO i=1,50
                              WRITE(6,*) ' i=',i,'   G=', &
                                   g(1,i),g(2,i),g(3,i)
                           ENDDO
                           CALL errore('c_phase_z2','wrong g',1)
                        ENDIF
                     ELSE
                        WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                             gtr(1),gtr(2),gtr(3)
                        CALL errore('c_phase_z2','wrong gtr',1)
                     END IF
                     map_g(ig)=ng
                  END DO
               END IF

!              --- Matrix elements calculation ---

               mat(:,:) = (0.d0, 0.d0)
               DO mb=1,nbnd
                  IF ( .NOT. l_cal(mb) ) THEN
                      mat(mb,mb)=(1.d0, 0.d0)
                  ELSE
                     aux(:) = (0.d0, 0.d0)
                     IF (kpar /= nppstr) THEN
                        DO ig=1,npw1
                           aux(igk1(ig))=evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux(igk1(ig)+ngm)=evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                     ELSEIF (.NOT. l_para) THEN
                        DO ig=1,npw1
                           aux(map_g(ig))=evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux(map_g(ig)+ngm)=evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                     ELSE
!
!   In this case this processor might not have the G-G_0
!
                        aux_g=(0.d0,0.d0)
                        DO ig=1,npw1
                           aux_g(mapgm_global(ig_l2g(igk1(ig)),gdir)) &
                                                =evc(ig,mb)
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,npw1
                              aux_g(mapgm_global(ig_l2g(igk1(ig)),gdir) &
                                                + ngm_g) =evc(ig+npwx,mb)
                           ENDDO
                        ENDIF
                        CALL mp_sum(aux_g(:), intra_bgrp_comm )
                        DO ig=1,ngm
                           aux(ig) = aux_g(ig_l2g(ig))
                        ENDDO
                        IF (noncolin) THEN
                           DO ig=1,ngm
                              aux(ig+ngm) = aux_g(ig_l2g(ig)+ngm_g)
                           ENDDO
                        ENDIF
                     ENDIF
!
                     DO nb=1,nbnd
                        IF ( l_cal(nb) ) THEN
                           aux0(:)= (0.d0, 0.d0)
                           DO ig=1,npw0
                              aux0(igk0(ig))=psi(ig,nb)
                           END DO
                           IF (noncolin) THEN
                              DO ig=1,npw0
                                aux0(igk0(ig)+ngm)=psi(ig+npwx,nb)
                              END DO
                           ENDIF
                           mat(nb,mb) = zdotc (ngm*npol,aux0,1,aux,1)
                        END IF
                     END DO
                  END IF
               END DO
               !
               call mp_sum( mat, intra_bgrp_comm )
               !
               DO nb=1,nbnd
                  DO mb=1,nbnd
!                    --- Calculate the augmented part: ij=KB projectors, ---
!                    --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
!                    --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
!                    --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                     IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                        if (okvan) then
                           pref = (0.d0,0.d0)
                           DO jkb=1,nkb
                              nhjkb = nkbtonh(jkb)
                              na = nkbtona(jkb)
                              np = ityp(na)
                              nhjkbm = nh(np)
                              jkb1 = jkb - nhjkb
                              DO j = 1,nhjkbm
                                 IF (noncolin) THEN
                                    IF (lspinorb) THEN
                                       pref = pref+(CONJG(becp0%nc(jkb,1,nb))* &
                                                  becp_bp%nc(jkb1+j,1,mb)  &
                                            *q_dk_so(nhjkb,j,1,np)   &
                                            +CONJG(becp0%nc(jkb,1,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb)  &
                                            *q_dk_so(nhjkb,j,2,np) &
                                            +CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,1,mb)  &
                                            *q_dk_so(nhjkb,j,3,np) &
                                            +CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb)   &
                                            *q_dk_so(nhjkb,j,4,np))*struc(na)
                                    ELSE
                                       pref = pref+(CONJG(becp0%nc(jkb,1,nb))* &
                                                   becp_bp%nc(jkb1+j,1,mb) + &
                                             CONJG(becp0%nc(jkb,2,nb))* &
                                                   becp_bp%nc(jkb1+j,2,mb))  &
                                             *q_dk(nhjkb,j,np)*struc(na)
                                    END IF
                                 ELSE
                                    pref = pref+CONJG(becp0%k(jkb,nb))* &
                                           becp_bp%k(jkb1+j,mb) &
                                      *q_dk(nhjkb,j,np)*struc(na)
                                 END IF
                              ENDDO
                           ENDDO
                           mat(nb,mb) = mat(nb,mb) + pref
                        endif
                     endif
                  ENDDO
               ENDDO

!  -------------------------------------------------------------------------   !
!  THE FOLLOWING CODE IS TAKEN FROM z2_main.f90 BY ALEXEY SOLUYANOV
!  -------------------------------------------------------------------------   !
!              --- Calculate SVD of M matrix: M = UU * SV * VV^* ---
               allocate (work(1))
               lwork = -1
               call ZGESVD('A','A',nbnd,nbnd,aux,nbnd,SV,UU,nbnd,VT,nbnd,work,lwork,rwork,info)
               lwork = int(dble(work(1)))
               deallocate (work)

               allocate (work(lwork), rwork(5*nbnd))
               call ZGESVD('A','A',nbnd,nbnd,mat,nbnd,SV,UU,nbnd,VT,nbnd,work,lwork,rwork,info)
               deallocate (work, rwork)
               if (info > 0) call errore('c_phase_z2','error in SVD factorization, info > 0', info)
               if (info < 0) call errore('c_phase_z2','error in SVD factorization, info < 0', -info)

               write(stdout,*)
               write(stdout,'(5X,''kort,kpar='',2I4,10X,''k,kp='',2I4)') kort, kpar, kpoint-1, kpoint
               write(stdout,'(5X,''sigmas'')')
               write(stdout,'(''  '',8F9.4)') (SV(nb), nb=1,nbnd)
               ! test for Mmn threshold
               do nb = 1, nbnd
                  if (SV(nb) < m_threshold) call errore('c_phase_z2', 'k-point string too coarse', -1)
               enddo
               mat = matmul(UU, VT)
               lambda = matmul(lambda, mat)
!           --- End of dot products between wavefunctions and betas ---
            ENDIF

!        --- End loop over parallel k-points ---
         END DO  ! kpar

         ! diagonalize lamdba matrix
         lwork = -1
         allocate(work(1))
         call ZGEEV('N','N',nbnd,aux,nbnd,eig,VT,nbnd,UU,nbnd,work,lwork,rwork,info)
         lwork = int(dble(work(1)))
         deallocate (work)

         allocate (work(lwork), rwork(2*nbnd))
         call ZGEEV('N','N',nbnd,lambda,nbnd,eig,VT,nbnd,UU,nbnd,work,lwork,rwork,info)
         deallocate (work, rwork)
         if (info > 0) call errore('c_phase_z2','error in diagonalization, info > 0', info)
         if (info < 0) call errore('c_phase_z2','error in diagonalization, info < 0', -info)
         write(stdout,*)
         write(stdout,*)
         write(stdout,'(5X,''========== kort='',I4,'' =========='')') kort
         write(stdout,'(5X,''eigenvalues (real part)'')')
         write(stdout,'(''  '',8F9.4)') (dble(eig(nb)), nb=1,nbnd)
         write(stdout,'(5X,''eigenvalues (imag part)'')')
         write(stdout,'(''  '',8F9.4)') (dimag(eig(nb)), nb=1,nbnd)
         do nb = 1, nbnd
            if (dreal(cdlog(eig(nb)))>1d-8) print*, 'eigenvalue', nb, 'has real part'
            zz(nb) = dimag(cdlog(eig(nb))) / (2.d0*PI)
            !print*, zz(nb)
            print*, nb, cdlog(eig(nb)), zz(nb)
         enddo
         write(70,*) kort, (zz(nb), nb=1,nbnd)
         ! putting eigenvalues within the [-0.5,0.5) window
         print*, 'wrap'
         do nb = 1, nbnd
             i = int(zz(nb))
             zz(nb) = zz(nb) - i
             if (zz(nb) <= -0.5d0) zz(nb) = zz(nb) + 1.d0
             if (zz(nb) > 0.5d0) zz(nb) = zz(nb) - 1.d0
            print*, nb, zz(nb)
         enddo
         write(71,*) kort, (zz(nb), nb=1,nbnd)
         call hpsort(nbnd, zz, ind)
         print*, 'wrap+sort'
         do nb = 1, nbnd
            print*, nb, zz(nb)
         enddo
         write(72,*) kort, (zz(nb), nb=1,nbnd)
         ! gaps
         do nb = 1, nbnd
           if (nb < nbnd) then
              gaps(nb) = zz(nb+1)-zz(nb)
           else
              gaps(nb) = zz(1)+1.d0-zz(nb)
           endif
         enddo
         !print*, 'gaps'
         !do nb = 1, nbnd
         !   print*, nb, gaps(nb)
         !enddo
         ! find largest gap
         max_gap = maxval(gaps)
         do nb = 1, nbnd
            if (gaps(nb) == max_gap) kk = nb
         enddo
         if (kk /= nbnd) then
            point = (zz(kk+1) + zz(kk))/2.d0
         else
            point = (zz(1) + 1.d0 + zz(kk))/2.d0
         endif
         print*, nb+1, zz(1) + 1
         print*, nb+2, point

!     --- End loop over orthogonal k-points ---
      END DO  ! kort

!  --- End loop over spin ---
   END DO
   DEALLOCATE ( l_cal ) 

!  -------------------------------------------------------------------------   !
!                    electronic polarization: phase average                    !
!  -------------------------------------------------------------------------   !

!  --- Start loop over spins ---
   DO is=1,nspin_lsda

!  --- Initialize average of phases as complex numbers ---
      cave=(0.0_dp,0.0_dp)
      phik_ave=(0.0_dp,0.0_dp)

!     --- Start loop over strings with same spin ---
      DO kort=1,nkort

!        --- Calculate string index ---
         istring=kort+(is-1)*nkort

!        --- Average phases as complex numbers ---
         cave=cave+wstring(istring)*cphik(istring)

!     --- End loop over strings with same spin ---
      END DO

!     --- Get the angle corresponding to the complex numbers average ---
      theta0=atan2(AIMAG(cave), DBLE(cave))
!     --- Put the phases in an around theta0 ---
      DO kort=1,nkort
        istring=kort+(is-1)*nkort
        cphik(istring)=cphik(istring)/cave
        dtheta=atan2(AIMAG(cphik(istring)), DBLE(cphik(istring)))
        phik(istring)=theta0+dtheta
        phik_ave=phik_ave+wstring(istring)*phik(istring)
      END DO

!     --- Assign this angle to the corresponding spin phase average ---
      IF (nspin == 1) THEN
         phiup=phik_ave !theta0+dtheta
         phidw=phik_ave !theta0+dtheta
      ELSE IF (nspin == 2) THEN
         IF (is == 1) THEN
            phiup=phik_ave !theta0+dtheta
         ELSE IF (is == 2) THEN
            phidw=phik_ave !theta0+dtheta
         END IF
      ELSE IF (nspin==4 ) THEN
         phiup=phik_ave
         phidw=0.0_DP
      END IF

!  --- End loop over spins
   END DO


!  -------------------------------------------------------------------------   !
!                           write output information                           !
!  -------------------------------------------------------------------------   !

!  --- Information about the k-points string used ---
   WRITE( stdout,"(/,21X,'K-POINTS STRINGS USED IN CALCULATIONS')")
   WRITE( stdout,"(21X,37('~'),/)")
   WRITE( stdout,"(7X,'G-vector along string (2 pi/a):',3F9.5)") &
           gpar(1),gpar(2),gpar(3)
   WRITE( stdout,"(7X,'Modulus of the vector (1/bohr):',F9.5)") &
           gvec
   WRITE( stdout,"(7X,'Number of k-points per string:',I4)") nppstr
   WRITE( stdout,"(7X,'Number of different strings  :',I4)") nkort

!  --- End of information relative to polarization calculation ---
   WRITE( stdout,"(/,/,15X,50('=')/,/)")

!  -------------------------------------------------------------------------   !
!                                  finalization                                !
!  -------------------------------------------------------------------------   !

!  --- Free memory ---
   DEALLOCATE(wstring)
   DEALLOCATE(cphik)
   DEALLOCATE(loc_k)
   DEALLOCATE(phik)
   DEALLOCATE(ln)
   DEALLOCATE(aux)
   DEALLOCATE(aux0)
   DEALLOCATE(psi)
   IF (l_para) THEN
      DEALLOCATE ( aux_g )
   ELSE
      DEALLOCATE ( map_g )
   ENDIF
   IF (okvan) THEN
      CALL deallocate_bec_type ( becp0 )
      CALL deallocate_bec_type ( becp_bp )
      IF (lspinorb) DEALLOCATE(q_dk_so)
   END IF



END SUBROUTINE c_phase_z2

