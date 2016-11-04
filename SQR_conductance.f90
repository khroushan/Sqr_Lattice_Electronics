!=============================================================
!- Program to calculate the conductivity of a square lattice
!- using the recerieve Green function. 
! =========          ===========              ==============
! A. Ahmadi           
! A. Ahmadi           2/12/2015        trace is defined as a function
! A. Ahmadi           3/11/2016        more comments and module base 
!-============================================================
module numeric_kind
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: sp = kind(1.0), &
                        dp = selected_real_kind(2 * precision(1.0_sp) ), &
                        qp = selected_real_kind(2 * precision(1.0_dp) )
  real(kind=dp), parameter :: pi = 3.14159265358979323
contains
  !==================================
  subroutine init_random_seed()
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed
    
    call random_seed(size = n)
    allocate(seed(n))
    
    call system_clock(count=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    
    deallocate(seed)
    
  end subroutine init_random_seed
end module numeric_kind

! ==========================================================
! ===                Linear Algebra Routines             ===
! ==========================================================
module linalg_routines
  use numeric_kind

contains
  
  subroutine sqr_invers(n, a_mtrx, inv_mtrx)

    implicit none

    integer, intent(in) :: n
    integer ::  info, lda, lwork, nb
    integer, dimension(n) :: ipiv
    complex(kind=dp), intent(in), dimension(n,n) :: a_mtrx
    complex(kind=dp), intent(out), dimension(n,n) :: inv_mtrx
    complex(kind=dp), allocatable, dimension(:) :: work
    integer :: ilaenv

    lda = n !- leading dimension of array
    nb = ILAENV( 1, 'ZGETRI', '', n, n, -1, -1)
    lwork = n * nb
    !print *, "lwork is :", lwork
    allocate(work(lwork))
    !- ipiv(output) integer array, dimension(min(m,n))
    !- the pivot indicies; for 1 <= i <= min(m,n), row
    !- i of the matrix was interchanged with row ipiv(i).
    inv_mtrx = a_mtrx
    !===================
    !- LU factorization
    !===================
    call zgetrf(n, n, inv_mtrx, lda, ipiv, info)
    
    if (info > 0 ) then
       print *, "factor is singular :-("
    elseif (info < 0) then
       print *, "the ", info, &
            "th argument had an illegal value :-("
    end if
    !===================
    !- matrix inversion
    !===================
    call zgetri(n, inv_mtrx, lda, ipiv, work, lwork, info)
    
    if (info > 0 ) then
       print *, "the U(i,i), i =  ", info, &
            "is zero, singular matrix :-("
    elseif (info < 0) then
       print *, "the ", info, &
            "th argument had an illegal value :-("
    end if

    !print *, "for optimal performance Lwork = ", work(1)

    deallocate(work)
  end subroutine sqr_invers
  
end module linalg_routines

! ==========================================================
! ===                Lead Green's Function               ===
! ==========================================================
module lead_routines
  !====================
  !- The subroutine to calculate the Green's 
  !- function of the left and right leads
  !====================
  use numeric_kind
  use linalg_routines
contains
  
  subroutine lead(e, g_l, h, mlat)
    implicit none
    
    integer :: mlat ! the number of vertical sites in one slice  
    
    complex(kind=dp), dimension(mlat,mlat) :: h, u 
    complex(kind=dp), dimension(mlat,mlat) :: aux1,aux2,aux3, g_l
    
    real(kind=dp) :: e
    
    complex(kind=dp), parameter  :: eta = (0.d0,1.d0)*1.0d-8
    complex(kind=dp), dimension(mlat,mlat) :: Gamma ,G   
    
    integer :: i, j, l1,l2
  
    h = (0.d0,0.d0) !Initializing
    
    u = (0.d0,0.d0) !Initializing
  
    !====================
    !- Constructing the Hamiltonian and Coupling matrices
    !====================
    do i = 1, mlat - 1
       h(i,i+1) = (-1.d0,0.d0)
       h(i+1,i) = conjg(h(i,i+1))
    end do
  
    !===================
    do i = 1, Mlat
       u(i,i) = (-1.d0,0.d0)
    end do
    !============
  
    g_l = (0.d0,0.d0)
  
    !============================================
    !- Loop to calculate the lead Green's function
    !=============================================
 
    call glead_dec(e,h,u,g_l,mlat)
  
  end subroutine lead
  ! ============================================================

  ! ============================================================
  ! Standard subroutine to compute the contact (surface) Green's function
  ! of a semi-infinite lead using the decimation method and a 
  ! tight-binding (nearest-neighbor) representation.
  !================================
  !   DATE            PROGRAMER                  DESCRIPTION OF CHANGE
  !  =======         ===========                =======================
  ! 3/18/2012        A. Ahmadi                    
  !
  subroutine glead_dec(e_l,h,u, g, mlat)
    
    implicit none
    
    ! * input/output variables 
    integer :: mlat ! the number of sites in vertical slice
    complex(kind=dp), intent(out), dimension(mlat,mlat) :: g
    complex(kind=dp), intent(in), dimension(mlat,mlat) :: u,h
    
    complex(kind=dp), parameter :: eta = (0.d0,1.d0)*1.0d-8
    complex(kind=dp), dimension(mlat,mlat) :: aux1, aux2, aux3
    complex(kind=dp), dimension(mlat,mlat) :: alpha, beta, eps, eps_s


    != local variables 
    integer :: l1, l2, k
    complex(kind=dp) :: ee
    real(kind=dp) :: e_l

    !to pass true dimension matrix to inversion routine
    complex(kind=dp), dimension(mlat,mlat) :: mtx1, mtx2 
    
    !========================
    ! - matrix initialization
    !========================
    ee = (1.d0,0.d0)*e_l + eta

    beta = conjg(transpose(u))
    alpha = u
    eps = h
    eps_s = h

    !=============================
    ! - starts the decimation loop
    !=============================
    do k = 1, 40
       !===========================================
       ! - construct: aux3 = (E - eps_{k-1})**(-1)
       !===========================================

       do l2 = 1, mlat
          do l1 = 1, mlat
             aux1(l1,l2) = (0.d0,0.d0)
             aux2(l1,l2) = - eps(l1,l2)
          end do
          aux1(l2,l2) = (1.d0,0.d0)
          aux2(l2,l2) = ee - eps(l2,l2)
       end do

       call sqr_invers(mlat, aux2, aux3)

       !================================================
       ! - construct: aux2 = alpha_{k-1}*aux3*beta_{k-1}
       !================================================

       aux2 = matmul(matmul(alpha,aux3),beta)

       !======================================
       ! updates decimation matrix: eps_s & eps
       !======================================

       eps_s = eps_s + aux2
       eps = eps + aux2

       !===============================================
       ! - construct: aux2 = beta_{k-1}*aux3*alpha_{k-1}
       !===============================================

       aux2 = matmul(matmul(beta,aux3),alpha)

       !===============================
       ! updates decimation matrix: EPS
       !===============================

       eps = eps + aux2

       !=================================================
       ! - construct: aux2 = alpha_{k-1}*aux3*alpha_{k-1}
       !=================================================

       aux2 = matmul(matmul(alpha,aux3),alpha)

       !==============================
       ! update decimation matrix ALPHA
       !==============================

       alpha = aux2

       !===============================================
       ! - construct: aux2 = beta_{k-1}*aux3*beta_{k-1}
       !===============================================

       aux2 = matmul(matmul(beta,aux3),beta)

       !==============================
       ! update decimation matrix BETA
       !==============================

       beta = aux2     

       !=========================
       ! - end of decimation loop
       !=========================
    end do

    !============
    ! - compute G
    !============
    
    do l2 = 1, mlat
       do l1 = 1, mlat
          aux1(l1,l2) = - eps_s(l1,l2)
          aux2(l1,l2) = (0.d0,0.d0)
       end do
       aux1(l2,l2) = ee + aux1(l2,l2)
       aux2(l2,l2) = (1.d0,0.d0)
    end do
    
    call sqr_invers(mlat, aux1, g) 
    
    
  end subroutine glead_dec
end module lead_routines

! ==========================================================
! ===                Constructing Hamiltonian            ===
! ==========================================================
module hamil
  use numeric_kind
  
contains
  subroutine hamilt(mlat, nlat, b_ext, imp_perc, v_imp, xi, h_d, tau_l, tau_r)
    implicit none
    
    integer, intent(in) :: mlat, nlat
    integer :: n_imp
    real(kind=sp), intent(in)  :: b_ext, imp_perc, v_imp, xi ! b_ext in tesla unit
    real(kind=sp) :: b_cof = 0.0
    complex(kind=dp), intent(out), dimension(mlat*nlat, mlat*nlat) :: h_d
    complex(kind=dp), dimension(mlat*nlat, mlat*nlat) :: v 
    real(kind=sp) , allocatable, dimension(:) :: aux_x, aux_y
    integer, allocatable, dimension(:) :: x_imp, y_imp
    complex(kind=dp), intent(out) , dimension(mlat, mlat*nlat):: tau_l, tau_r
    integer :: m1, m2, i !dummy counter
    real(kind=dp), parameter :: a_0 = 0.142e-9, hbar = 1.05e-34, e = 1.6e-19
    ! a_0 lattice constant all in si unit
    
    
    h_d = (0.d0, 0.d0); tau_l = (0.d0, 0.d0); tau_r = (0.d0, 0.d0); v = (0.d0,0.d0)
    b_cof = (b_ext * e * a_0**2)/(2.0 * hbar )
    !=========================
    !- constructing hamiltonian
    !=========================
    vertical :do m1 = 0, nlat - 1
       do m2 = 1, mlat - 1
          h_d( m1*mlat + m2, m1*mlat + m2 +1 ) = (-1.d0, 0.d0) * &
               exp( (-1.d0, 0.d0)*m1*b_cof  )
          h_d( m1*mlat + m2 +1, m1*mlat + m2 ) = &
               conjg(h_d( m1*mlat + m2, m1*mlat + m2 +1 ) )
       end do
    end do vertical

    horizintal: do m1 = 0, nlat - 2
       do m2 = 1, mlat
          h_d( m1*mlat + m2, (m1 + 1)*mlat  + m2 ) = (-1.d0, 0.d0) * &
               exp( (1.d0, 0.d0)*m2*b_cof  )
          h_d( (m1 + 1)*mlat  + m2 , m1*mlat + m2) = &
               conjg( h_d( m1*mlat + m2, (m1 + 1)*mlat  + m2 ) )
       end do
    end do horizintal

    !====================================
    !- potential due to gassian impurities
    !====================================

    n_imp = int(imp_perc*mlat*nlat / 100)

    allocate( aux_x(n_imp), aux_y(n_imp), x_imp(n_imp), y_imp(n_imp) )
    call random_number(aux_x)
    x_imp = mod(int(aux_x*nlat), nlat)
    call random_number(aux_y)
    y_imp = mod(int(aux_y*mlat), mlat) + 1

    do m1 = 0, nlat - 1
       do m2 = 1, mlat

          do i = 1, n_imp !sum over all impurities

             v(m1*mlat+m2, m1*mlat+m2) =  v(m1*mlat+m2, m1*mlat+m2) + v_imp * &
                  exp( ((m1 - x_imp(i) )*mlat + ( m2 - y_imp(i) ))**2 * a_0**2 &
                  / xi**2 ) 

          end do
       end do
    end do

    deallocate(aux_x, aux_y, x_imp, y_imp )

    h_d = h_d + v ! update the hamiltonian with impurities potential

    !====================================
    !- constructing left nd right contacts
    !====================================
    do m1 = 1, mlat
       tau_l(m1, m1) = (-1.d0, 0.d0)
       tau_r(m1, m1 + (nlat - 1) * mlat) = (-1.d0, 0.d0)
    end do

  end subroutine hamilt
end module hamil
! ============================================================

program square_conductivity
  
  use numeric_kind
  use linalg_routines
  use lead_routines
  use hamil
  
  implicit none
  
  
  integer :: mlat, nlat, dim,  n_energies ! the number of sites in vertical slice and the width of device
  real(kind=sp) :: b_ext, imp_perc, v_imp, xi
  integer :: ie , m1, m2, i_avg
  real(kind=dp) :: e
  complex(kind=dp), allocatable, dimension(:,:) :: gamma, g, g_l, h, sigma_l, sigma_r, &
       csigma_l, csigma_r, gamma_l, gamma_r
  complex(kind=dp), allocatable, dimension(:,:) :: h_d, tau_l, tau_r, aux1
  real(kind=dp), allocatable, dimension(:) :: dos
  real(kind=dp) :: g1, f1, g1_avg = 0.d0
  real(kind=dp) :: e_l
  real(kind=dp) :: cond  
  complex(kind=dp), parameter :: eta = (0.d0,1.d0)*1.0d-8

  !==============================
  !- initializing the input values
  !==============================
  open(unit=1, file='input.dat', action ='read')
  read(1,*) mlat
  read(1,*) nlat
  read(1,*) n_energies
  read(1,*) b_ext
  read(1,*) imp_perc
  read(1,*) v_imp
  read(1,*) xi
  close(unit=1)

  call init_random_seed()
  dim = nlat*mlat

  allocate( h(mlat,mlat), g_l(mlat,mlat), dos(dim) )
  allocate(h_d(dim, dim), tau_l(mlat, dim), tau_r(mlat, dim), gamma_l(dim,dim),  &
       gamma_r(dim,dim), g(dim, dim) )
  allocate(sigma_l(dim,dim), sigma_r(dim,dim), csigma_l(dim,dim), &
       csigma_r(dim,dim),  aux1(dim, dim) )
 
  open(unit=14, file="gf", action='write')
  open(unit=16, file="dos", action='write')
  energy :do ie = 0, n_energies
     g1_avg = 0.d0; dos = 0.d0
     e = -2.0d0 + (ie)*((4.d0)/n_energies) 
     
     avg :  do i_avg = 0, 100
        
        h_d = (0.d0, 0.d0) ; tau_l = (0.d0, 0.d0); tau_r = (0.d0, 0.d0)
     
        call hamilt(mlat, nlat, b_ext, imp_perc, v_imp, xi,  h_d, tau_l, tau_r)
                     
        call lead(e, g_l, h, mlat)
        
        sigma_r = matmul(transpose(tau_l),matmul(g_l,tau_l))
        sigma_l = matmul(transpose(tau_r),matmul(g_l,tau_r))
        
        
        csigma_l = transpose(dconjg(sigma_l))
        csigma_r = transpose(dconjg(sigma_r))
        
        
        gamma_l = (0.d0,1.d0)*(sigma_l - csigma_l)
        gamma_r = (0.d0,1.d0)*(sigma_r - csigma_r)
        
        !============================
        !- calculate the g.f. of bulk
        !============================
        do m1 = 1, dim
           do m2 = 1, dim
              aux1(m1,m2) = -h_d(m1,m2) - sigma_l(m1,m2) - sigma_r(m1,m2)
           end do
           aux1(m1,m1) = e + eta + aux1(m1,m1)
        end do
        
        
        call sqr_invers(dim, aux1, g) 
        
        g1 = cond(gamma_l, gamma_r, g, dim)
     
        g1_avg = g1_avg + g1
     
        do m1 = 0, nlat-1
           do m2 = 1, mlat
              dos(m1*mlat+m2) = dos(m1*mlat+m2) +  (-1.0/pi) * &
                   aimag(g(m1*mlat + m2, m1*mlat + m2 ))
           end do
        end do

        
     end do avg
     
     write(14,1)ie, e, g1/100.0
1    format(i3, 2x, f6.3, 2x, f6.3 ) 

     do m1 = 0, nlat-1
        do m2 = 1, mlat
           write(16,3) m1, m2, (1.0/100.0)*dos(m1*mlat +m2)
        end do
     end do
     
3    format(i2,2x,i2,2x, f6.3)
     
  end do energy
     


  close(unit=14)
  close(unit=16)
  
  deallocate( h, g_l, h_d, tau_l, tau_r, gamma_l, gamma_r, g, sigma_l, sigma_r, &
       csigma_l, csigma_r, aux1 )
end program square_conductivity

!==============================
pure function cond(gammal, gammar, gl_ext, dim)
  !================================
  !- return the conductance as the 
  !- trace of gammal*g*gammar*g_dg
  !================================
  implicit none
  !==============
  ! * input/output
  !==============
  integer, intent(in) :: dim
  complex*16, intent(in), dimension(dim, dim) :: gammal, gammar, gl_ext
  double precision :: cond
  complex*16, dimension(dim,dim) :: aux1, aux2
  complex*16 ::  trg
  integer :: m1, m2, m3
  intrinsic dconjg

  aux1 = dconjg(transpose(gl_ext))
  
  aux2 = matmul(gammal, matmul(gl_ext, matmul(gammar, aux1) ) )
  
  trg = (0.d0, 0.d0)

  do m1=1, dim

     trg = trg + aux2(m1, m1)
  end do
  cond = real(trg)
end function
!===========================








