program PowerReconstr
    ! Program for the one-dimensional power reconstruction calculation
    ! based on the Analitical Function Expansion Nidal method 
    ! Last Update                 :: Vlad, 23.04.2021 
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    real                          :: keff ! Effective multiplication factor
    real                          :: S12
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(NG)        :: Fav ! Average flux in node
    real,    dimension(NG)        :: Fcorn ! Corner flux
    real,    dimension(NG)        :: Fs ! Surface flux
    real,    dimension(NG)        :: Fa ! Surface flux
    real,    dimension(NG)        :: Fb ! Surface flux
    real,    dimension(NG)        :: D ! Diffusion coefficient
    real,    dimension(NG)        :: NuSf ! Nuf * Fission cross-section
    real,    dimension(NG)        :: Sa ! Absorbtion cross-section
    real,    dimension(NG)        :: k0 ! Flux vector
    real,    dimension(NG)        :: k1 ! Flux vector
    real,    dimension(NG)        :: k2 ! Flux vector
    real,    dimension(NG)        :: flux_result ! Flux vector
    real,    dimension(NG*3)      :: f ! Flux vector
    real,    dimension(NG,NG)     :: B2
    real,    dimension(NG,NG)     :: test1
    real,    dimension(NG,NG)     :: test11
    real,    dimension(NG)        :: test2
    real,    dimension(NG,NG)     :: xappa ! alpha from the analytical expression, xappa = SQRT(B2) 
    real,    dimension(NG,NG)     :: Matr_Func_xappa
    real,    dimension(NG,NG)     :: MatFunc_Sin
    real,    dimension(NG,NG)     :: MatFunc_Cos
    real,    dimension(3*NG,3*NG) :: SensitivityMatrix
    real,    dimension(3*NG,3*NG) :: SensitivityMatrix1
    real,    dimension(NG,NG, 9)  :: Matr_Func
    integer                       :: i,j ! counters
    real :: a,x
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    ! lambda - generic function for sqrt calculation 
    ! lambda = sqrt(x) if x>0, lambda = sqrt(-x) if x<0  
    
    ! B2 matrix definition using the results of Fenics calculations
    ! and mathematica B2 matrix expression
    !data B2 / 3.30358, -14.2291, 0., -6.91217/ ! NG = 2
    
    ! Funtionals definition obtained from Fenics calculations
    Fav(1) = 6.03E+16 !6.80E+16
    Fav(2) = 1.13E+16 !1.43E+16
    !Fcorn(1) = 0.010870036379003406
    !Fcorn(2) = 0.0023402839954890418
    !Fs(1) = 0.010686852064653892
    !Fs(2) = 0.0022345909901370685
    Fa(1) = 6.79E+16 !7.98E+16 !5.65E+16
    Fa(2) = 1.41E+16 !1.84E+16 !1.04E+16
    Fb(1) = 6.79E+16 !7.98E+16
    Fb(2) = 1.41E+16 !1.84E+16
    keff = 1.110327
    
    D(1) = 1.25478 
    D(2) = 0.233905
    nusf(1) = 0.00719674
    nusf(2) = 0.157254
    s12 = 0.0239031
    sa(1) = 0.0123492
    sa(2) = 0.130017
    
    f(1) = Fav(1)
    f(2) = Fav(2)
    f(3) = Fa(1)
    f(4) = Fa(2)
    f(5) = Fb(1)
    f(6) = Fb(2)
    
    ! calculation of the eigenvalues of B2 matrix (just for test)
    !call MAT_Compute_Eigenv_2x2(B2, eigenv)
    !print*, eigenv
    
    data xappa / 3.30358, -14.2291, 0., -6.91217/
    !data xappa / 13.2143, -56.9166, 0., -27.6487/
    
    ! calculation of the xappa matrix (xappa = SQRT(B2))
    !call MAT_Set_MatFunc_B2_ANM(B2, Matr_Func_xappa)
    !xappa = Matr_Func_xappa
    !print*, xappa
 
    ! calculation of matrix functions included in the Sensitivity matrix 
    call MAT_Set_MatFunc_ANM_1D(xappa, Matr_Func)
    !print*, Matr_Func
    
    ! forming the Sensitivity matrix using the results of the matrix functions calculation
    SensitivityMatrix = 0.   
    do i = 1, 3
        do j = 1, 3
            SensitivityMatrix((i*2)-1,(j*2)-1) = Matr_Func(1,1,(i-1)*3+j)
            SensitivityMatrix((i*2)-1,(j*2)) = Matr_Func(1,2,(i-1)*3+j)
            SensitivityMatrix((i*2),(j*2)-1) = Matr_Func(2,1,(i-1)*3+j)
            SensitivityMatrix((i*2),(j*2)) = Matr_Func(2,2,(i-1)*3+j)
        enddo
    enddo
    
    !print the result
    do i = 1,3*NG
        print'(9(1x,f5.2))', (SensitivityMatrix(i,j), j=1,3*NG)
    enddo

    !end
    
    SensitivityMatrix1 = SensitivityMatrix
    call MSC_LU_Solve(SensitivityMatrix1, NG*3, f)
    print*, f
    print*, matmul(SensitivityMatrix, f)
    
    
    k0(1)=f(1)
    k0(2)=f(2)
    k1(1)=f(3)
    k1(2)=f(4)
    k2(1)=f(5)
    k2(2)=f(6)
    
    
    data Unity / 1., 0., 0., 1./ 
    x = 0.0
    call MAT_Set_MatFunc_Sin_ANM(xappa, x, MatFunc_Sin)
    call MAT_Set_MatFunc_Cos_ANM(xappa, x, MatFunc_Cos)
    !print*,xappa*x
    flux_result = k0 + matmul(MatFunc_Sin, k1) + matmul(MatFunc_Cos, k2)
    
    print*,'Result:'
    print*,flux_result
    
    !print*, matmul(SensitivityMatrix, f)
    
    read*, a

end program PowerReconstr

    
subroutine MAT_Compute_Eigenv_2x2(xappa, eigenv_xappa)
!*=====================================================================*
!* Compute Eigenvalue of 2x2 matrix                                    *
!* is used to compute eigenvalue of buckling matrix in the case        *
!*  of 2 group ANM                                                     *
!* Last Update              Slava (c) 5.06.1998                        *
!*=====================================================================*
      implicit none
      integer NG
      parameter(NG=2)
!* Input: xappa(2,2) - matrix
      real xappa(NG,NG)
!* Output: eigenv_xappa(2) - eigenvalues
      real eigenv_xappa(NG)
!* Local Variables:
!*      determ - Determinant of the matrix
      real determ, d

      determ = xappa(1,1)*xappa(2,2) - xappa(1,2)*xappa(2,1)
      d = (xappa(1,1)+xappa(2,2))**2 - 4.*determ
      eigenv_xappa(2)  = 0.5*(xappa(1,1) + xappa(2,2) + sqrt(d))
      eigenv_xappa(1)  = 0.5*(xappa(1,1) + xappa(2,2) - sqrt(d))

      return
end    

    
subroutine MAT_Set_MatFunc_B2_ANM(xappa, Matr_Func)
    ! Calculation of the Matrix Function SQRT(B2)=xappa  
    ! Last Update                 :: Vlad, 23.04.2021
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    
    ! Input variables
    real,    dimension(NG,NG)     :: xappa ! buckling matrix
    
    !Local variables
    real                          :: a1 ! temporary value
    real                          :: a2 ! temporary value
    real                          :: denom ! temporary value
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(NG)        :: gamma ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    integer                       :: n ! counter
    integer                       :: m ! counter
    
    ! Output variables
    real,    dimension(NG,NG)     :: Matr_Func ! Matrix Function SQRT(B2) = xappa
    
    data Unity / 1., 0., 0., 1./ 
 
    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            gamma(n) = SQRT(-eigenv(n))
        else
            gamma(n) = SQRT(eigenv(n))
        end if
    end do

    ! calculation of the matrix function 
    denom = Eigenv(1) - Eigenv(2)
    a1 = (MAT_F1B2_ANM(gamma(1),eigenv(1))-MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom                        
    a2 = (-Eigenv(2)*MAT_F1B2_ANM(gamma(1),eigenv(1)) &
        + Eigenv(1)*MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom

    do n = 1, NG
        do m = 1, NG
            Matr_Func(n,m) = xappa(n,m)*a1 + a2*Unity(n,m)
        end do
    end do
    
    contains

        real function MAT_F1B2_ANM(gamma,eigenv)
            ! Response function for the calculation of the xappa matrix  
            ! Last Update                 :: Vlad, 23.04.2021
            real                          :: eigenv ! eigenvalue, eigenv = x
            real                          :: gamma ! gamma = SQRT(x)

            MAT_F1B2_ANM = gamma
    
        end function MAT_F1B2_ANM
end
    
    
subroutine MAT_Set_MatFunc_ANM_1D(xappa, Matr_Func)
    ! Calculation of the Matrix Functions included in the Sensitivity matrix  
    ! Last Update                 :: Vlad, 23.04.2021
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    
    ! Input variables
    real,    dimension(NG,NG)     :: xappa ! xappa matrix (alpha from the analitical expression)
    
    !Local variables
    real                          :: denom ! temporary value
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(9)         :: a1 ! temporary values
    real,    dimension(9)         :: a2 ! temporary values
    real,    dimension(NG)        :: lambda ! generic function 
                                           ! cosh(x) if x > 0; cos(x) if x < 0;
    real,    dimension(NG)        :: gamma ! generic function 
                                           ! sinh(x) if x > 0; sin(x) if x < 0;
    real,    dimension(NG)        :: beta  ! generic function 
                                           ! sinh(sqrt(x))/sqrt(x) if x > 0; sinh(sqrt(-x))/sqrt(-x) if x < 0;
    real,    dimension(NG)        :: alpha ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    integer                       :: n ! counter
    integer                       :: m ! counter
    integer                       :: i ! counter
    
    ! Output variables
    real,    dimension(NG,NG,9)   :: Matr_Func ! Array for Matrix Functions values storage
    
    data Unity / 1., 0., 0., 1./ 

    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up lambda and gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            alpha(n) = sqrt(-eigenv(n))
            lambda(n) = cos(alpha(n))
            gamma(n) = sin(alpha(n))
            beta(n) = sin(alpha(n))/alpha(n)
        else
            alpha(n) = sqrt(eigenv(n))
            lambda(n) = cosh(alpha(n))
            gamma(n) = sinh(alpha(n))
            beta(n) = sinh(alpha(n))/alpha(n)
        end if
        !print*,alpha(n)
    end do

    ! calculation of the matrix function 
    denom = Eigenv(1) - Eigenv(2)
    
    a1(1) = (F_01_01(gamma(1),eigenv(1)) &
        - F_01_01(gamma(2),eigenv(2)))/denom
    a1(2) = (F_01_02(gamma(1),eigenv(1)) &
        - F_01_02(gamma(2),eigenv(2)))/denom
    a1(3) = (F_01_03(beta(1),eigenv(1)) &
        - F_01_03(beta(2),eigenv(2)))/denom
      
    a1(4) = (F_02_01(gamma(1),eigenv(1)) &
        - F_02_01(gamma(2),eigenv(2)))/denom
    a1(5) = (F_02_02(gamma(1),eigenv(1)) &
        - F_02_02(gamma(2),eigenv(2)))/denom
    a1(6) = (F_02_03(lambda(1),eigenv(1)) &
        - F_02_03(lambda(2),eigenv(2)))/denom
      
    a1(7) = (F_03_01(gamma(1),eigenv(1)) &
        - F_03_01(gamma(2),eigenv(2)))/denom
    a1(8) = (F_03_02(gamma(1),eigenv(1)) &
        - F_03_02(gamma(2),eigenv(2)))/denom
    a1(9) = (F_03_03(lambda(1),eigenv(1)) &
        - F_03_03(lambda(2),eigenv(2)))/denom
                     
    a2(1) = (-Eigenv(2)*F_01_01(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_01_01(gamma(2),eigenv(2)))/denom
    a2(2) = (-Eigenv(2)*F_01_02(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_01_02(gamma(2),eigenv(2)))/denom
    a2(3) = (-Eigenv(2)*F_01_03(beta(1),eigenv(1)) &
        + Eigenv(1)*F_01_03(beta(2),eigenv(2)))/denom
      
    a2(4) = (-Eigenv(2)*F_02_01(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_02_01(gamma(2),eigenv(2)))/denom
    a2(5) = (-Eigenv(2)*F_02_02(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_02_02(gamma(2),eigenv(2)))/denom
    a2(6) = (-Eigenv(2)*F_02_03(lambda(1),eigenv(1)) &
        + Eigenv(1)*F_02_03(lambda(2),eigenv(2)))/denom
      
    a2(7) = (-Eigenv(2)*F_03_01(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_03_01(gamma(2),eigenv(2)))/denom
    a2(8) = (-Eigenv(2)*F_03_02(gamma(1),eigenv(1)) &
        + Eigenv(1)*F_03_02(gamma(2),eigenv(2)))/denom
    a2(9) = (-Eigenv(2)*F_03_03(lambda(1),eigenv(1)) &
        + Eigenv(1)*F_03_03(lambda(2),eigenv(2)))/denom
      
    do i = 1, 9
        do n = 1, NG
            do m = 1, NG
                Matr_Func(n,m,i) = xappa(n,m)*a1(i) + a2(i)*Unity(n,m)
            end do
        end do
    end do


    contains
        ! Matrix functions expressions via  generic (lambda and gamma) functions
        ! first number - number of row in the Sensitivity matrix
        ! second number - number of column in the Sensitivity matrix
        ! Sensitivity matrix has the following form:
        ! 1     0                Sin(xappa)/(xappa) !
        ! 1     -Sin(xappa)      Cos(xappa)         !
        ! 1     Sin(xappa)       Cos(xappa)         !
        real function F_01_01(gamma,eigenv)
            real gamma, eigenv

            F_01_01 = 1.

        end function
        
        real function F_01_02(gamma,eigenv)
            real gamma, eigenv

            F_01_02 = 0.

        end function
        
        real function F_01_03(beta,eigenv)
            real beta, eigenv

            F_01_03 = beta
            !print*,F_01_03, eigenv

        end function
        
        real function F_02_01(gamma,eigenv)
            real gamma, eigenv

            F_02_01 = 1.

        end function
        
        real function F_02_02(gamma,eigenv)
            real gamma, eigenv

            F_02_02 = - gamma

        end function
        
        real function F_02_03(lambda,eigenv)
            real lambda, eigenv

            F_02_03 = lambda
            !print*,F_02_03, eigenv

        end function
        
        real function F_03_01(gamma,eigenv)
            real gamma, eigenv

            F_03_01 = 1.

        end function
        
        real function F_03_02(gamma,eigenv)
            real gamma, eigenv

            F_03_02 = gamma

        end function
        
        real function F_03_03(lambda,eigenv)
            real lambda, eigenv

            F_03_03 = lambda

        end function

    end subroutine
    

subroutine MAT_Set_MatFunc_Sin_ANM(xappa, x, Matr_Func_Sin)
    ! Calculation of the Matrix Function SQRT(B2)=xappa  
    ! Last Update                 :: Vlad, 23.04.2021
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    
    ! Input variables
    real,    dimension(NG,NG)     :: xappa ! buckling matrix
    
    !Local variables
    real                          :: a1 ! temporary value
    real                          :: a2 ! temporary value
    real                          :: denom ! temporary value
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(NG)        :: gamma ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG)        :: alpha ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    integer                       :: n ! counter
    integer                       :: m ! counter
    
    ! Output variables
    real,    dimension(NG,NG)     :: Matr_Func_Sin ! Matrix Function SQRT(B2) = xappa
    real :: x
    
    data Unity / 1., 0., 0., 1./ 
 
    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            alpha(n) = sqrt(-eigenv(n))*x
            gamma(n) = sin(alpha(n))
        else
            alpha(n) = sqrt(eigenv(n))*x
            gamma(n) = sinh(alpha(n))
        end if
    end do

    ! calculation of the matrix function 
    denom = Eigenv(1) - Eigenv(2)
    a1 = (MAT_F1B2_ANM(gamma(1),eigenv(1))-MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom                        
    a2 = (-Eigenv(2)*MAT_F1B2_ANM(gamma(1),eigenv(1)) &
        + Eigenv(1)*MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom

    do n = 1, NG
        do m = 1, NG
            Matr_Func_Sin(n,m) = xappa(n,m)*a1 + a2*Unity(n,m)
        end do
    end do
    
    contains

        real function MAT_F1B2_ANM(gamma,eigenv)
            ! Response function for the calculation of the xappa matrix  
            ! Last Update                 :: Vlad, 23.04.2021
            real                          :: eigenv ! eigenvalue, eigenv = x
            real                          :: gamma ! gamma = SQRT(x)

            MAT_F1B2_ANM = gamma
    
        end function MAT_F1B2_ANM
    end
    
subroutine MAT_Set_MatFunc_Cos_ANM(xappa, x, Matr_Func_Cos)
    ! Calculation of the Matrix Function SQRT(B2)=xappa  
    ! Last Update                 :: Vlad, 23.04.2021
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    
    ! Input variables
    real,    dimension(NG,NG)     :: xappa ! buckling matrix
    
    !Local variables
    real                          :: a1 ! temporary value
    real                          :: a2 ! temporary value
    real                          :: denom ! temporary value
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(NG)        :: gamma ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG)        :: alpha ! generic function 
                                           ! sqrt(x) if x > 0; sqrt(-x) if x < 0;
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    integer                       :: n ! counter
    integer                       :: m ! counter
    
    ! Output variables
    real,    dimension(NG,NG)     :: Matr_Func_Cos ! Matrix Function SQRT(B2) = xappa
    real :: x
    
    data Unity / 1., 0., 0., 1./ 
 
    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            alpha(n) = sqrt(-eigenv(n))*x
            gamma(n) = cos(alpha(n))
        else
            alpha(n) = sqrt(eigenv(n))*x
            gamma(n) = cosh(alpha(n))
        end if
    end do

    ! calculation of the matrix function 
    denom = Eigenv(1) - Eigenv(2)
    a1 = (MAT_F1B2_ANM(gamma(1),eigenv(1))-MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom                        
    a2 = (-Eigenv(2)*MAT_F1B2_ANM(gamma(1),eigenv(1)) &
        + Eigenv(1)*MAT_F1B2_ANM(gamma(2),eigenv(2)))/denom

    do n = 1, NG
        do m = 1, NG
            Matr_Func_Cos(n,m) = xappa(n,m)*a1 + a2*Unity(n,m)
        end do
    end do
    
    contains

        real function MAT_F1B2_ANM(gamma,eigenv)
            ! Response function for the calculation of the xappa matrix  
            ! Last Update                 :: Vlad, 23.04.2021
            real                          :: eigenv ! eigenvalue, eigenv = x
            real                          :: gamma ! gamma = SQRT(x)

            MAT_F1B2_ANM = gamma
    
        end function MAT_F1B2_ANM
    end
    

subroutine MSC_LU_Solve(A1, N, B1)
      implicit none 
!*=====================================================================*
!* Solution of the equation A1 * x = B1 Using LU Decomposition         * 
!* Last Update              Slava (c) April 1998                       *
!*=====================================================================*
!* Input: A1 - Matrix, N - Dimension, B1 - Right Side
      integer N
      real A1(N,N),B1(N)
!* Output: A1 - LU Decomposition of the MAtrix A1, B - Solution
!* Local Variables: 
      integer k, i, j
      real sum 

!* LU Decomposition
      do j = 1,N
!* first part beta_ij expression (2.3.12)
         do i= 1, j
            sum = a1(i,j)
            do k = 1, i-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
            a1(i,j) = sum
         end do
!* second part alfa_ij expression (2.3.13)
       do i = j+1,N
            sum = a1(i,j)
            do k = 1, j-1
               sum = sum - a1(i,k)*a1(k,j)
            end do
         a1(i,j) = sum/a1(j,j)
       end do
      end do

!* Solution LU x = b

!* forward substitution
      do i = 1,N
         sum = b1(i)
         do j = 1, i-1
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i) = sum
       end do

!* backsubstitution

      do i = N, 1,-1
         sum = b1(i)
         do j = i+1,N
            sum = sum - a1(i,j)*b1(j)
         end do
         b1(i)=sum/a1(i,i)
      end do

      return
      end
     