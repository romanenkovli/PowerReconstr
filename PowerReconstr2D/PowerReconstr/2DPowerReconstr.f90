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
    real,    dimension(NG)        :: k3 ! Flux vector
    real,    dimension(NG)        :: k4 ! Flux vector
    real,    dimension(NG)        :: k5 ! Flux vector
    real,    dimension(NG)        :: k6 ! Flux vector
    real,    dimension(NG)        :: k7 ! Flux vector
    real,    dimension(NG)        :: k8 ! Flux vector
    real,    dimension(NG)        :: k9 ! Flux vector
    real,    dimension(NG)        :: k10 ! Flux vector
    real,    dimension(NG)        :: k11 ! Flux vector
    real,    dimension(NG)        :: k12 ! Flux vector
    real,    dimension(NG)        :: flux_result ! Flux vector
    real,    dimension(NG*13)      :: f ! Flux vector
    real,    dimension(NG,NG)     :: B2
    real,    dimension(NG,NG)     :: test1
    real,    dimension(NG,NG)     :: test11
    real,    dimension(NG)        :: test2
    real,    dimension(NG,NG)     :: xappa ! alpha from the analytical expression, xappa = SQRT(B2) 
    real,    dimension(NG,NG)     :: Matr_Func_xappa
    real,    dimension(NG,NG)     :: MatFunc_Sin1
    real,    dimension(NG,NG)     :: MatFunc_Sin2
    real,    dimension(NG,NG)     :: MatFunc_Sin3
    real,    dimension(NG,NG)     :: MatFunc_Sin4
    real,    dimension(NG,NG)     :: MatFunc_Sin5
    real,    dimension(NG,NG)     :: MatFunc_Sin6
    real,    dimension(NG,NG)     :: MatFunc_Cos1
    real,    dimension(NG,NG)     :: MatFunc_Cos2
    real,    dimension(NG,NG)     :: MatFunc_Cos3
    real,    dimension(NG,NG)     :: MatFunc_Cos4
    real,    dimension(NG,NG)     :: MatFunc_Cos5
    real,    dimension(NG,NG)     :: MatFunc_Cos6
    real,    dimension(13*NG,13*NG) :: SensitivityMatrix
    real,    dimension(13*NG,13*NG) :: SensitivityMatrix1
    real,    dimension(NG,NG,169)  :: Matr_Func
    integer                       :: i,j ! counters
    real :: a,x,y
    ! lambda - generic function for sqrt calculation 
    ! lambda = sqrt(x) if x>0, lambda = sqrt(-x) if x<0  
    
    ! B2 matrix definition using the results of Fenics calculations
    ! and mathematica B2 matrix expression
    data B2 / 3.21931, -13.7509, 0., -3.54266/ ! NG = 2
    
    ! Funtionals definition obtained from Fenics calculations
    Fav(1) = 0.010062599566301773 !6.03E+16 !6.80E+16
    Fav(2) = 0.0019200392896564895 !1.13E+16 !1.43E+16
    Fcorn(1) = 0.010870036379003406
    Fcorn(2) = 0.0023402839954890418
    Fs(1) = 0.010686852064653892
    Fs(2) = 0.0022345909901370685
    !Fa(1) = 6.79E+16 !7.98E+16 !5.65E+16
    !Fa(2) = 1.41E+16 !1.84E+16 !1.04E+16
    !Fb(1) = 6.79E+16 !7.98E+16
    !Fb(2) = 1.41E+16 !1.84E+16
    keff = 1.1547918497089225 !1.110327
    
    D(1) = 1.25478 
    D(2) = 0.233905
    nusf(1) = 0.00719674
    nusf(2) = 0.157254
    s12 = 0.0239031
    sa(1) = 0.0123492
    sa(2) = 0.130017
    
    f(1) = Fav(1)
    f(2) = Fav(2)
    f(3) = Fs(1)
    f(4) = Fs(2)
    f(5) = Fs(1)
    f(6) = Fs(2)
    f(7) = Fs(1)
    f(8) = Fs(2)
    f(9) = Fs(1)
    f(10) = Fs(2)
    f(11) = Fs(1)
    f(12) = Fs(2)
    f(13) = Fs(1)
    f(14) = Fs(2)
    
    f(15) = Fcorn(1)
    f(16) = Fcorn(2)
    f(17) = Fcorn(1)
    f(18) = Fcorn(2)
    f(19) = Fcorn(1)
    f(20) = Fcorn(2)
    f(21) = Fcorn(1)
    f(22) = Fcorn(2)
    f(23) = Fcorn(1)
    f(24) = Fcorn(2)
    f(25) = Fcorn(1)
    f(26) = Fcorn(2)
    ! calculation of the eigenvalues of B2 matrix (just for test)
    !call MAT_Compute_Eigenv_2x2(B2, eigenv)
    !print*, eigenv
    
    ! calculation of the xappa matrix (xappa = SQRT(B2))
    call MAT_Set_MatFunc_B2_ANM(B2, Matr_Func_xappa)
    xappa = Matr_Func_xappa
    print*, xappa
 
    ! calculation of matrix functions included in the Sensitivity matrix 
    call MAT_Set_MatFunc_ANM_1D(xappa, B2, Matr_Func)
    !print*, Matr_Func
    
    ! forming the Sensitivity matrix using the results of the matrix functions calculation
    SensitivityMatrix = 0.   
    do i = 1, 13
        do j = 1, 13
            SensitivityMatrix((i*2)-1,(j*2)-1) = Matr_Func(1,1,(i-1)*13+j)
            SensitivityMatrix((i*2)-1,(j*2)) = Matr_Func(1,2,(i-1)*13+j)
            SensitivityMatrix((i*2),(j*2)-1) = Matr_Func(2,1,(i-1)*13+j)
            SensitivityMatrix((i*2),(j*2)) = Matr_Func(2,2,(i-1)*13+j)
        enddo
    enddo
    
    !print the result
    do i = 1,13*NG
        print'(26(1x,f5.2))', (SensitivityMatrix(i,j), j=1,13*NG)
    enddo

    !end
    
    print*,f
    SensitivityMatrix1 = SensitivityMatrix
    call INVERT1(SensitivityMatrix, SensitivityMatrix1, 13*NG)
    call MSC_LU_Solve(SensitivityMatrix1, NG*13, f)
    !print*, 'f', f
    print*, matmul(SensitivityMatrix, f)
    
    
    k0(1)=f(1)
    k0(2)=f(2)
    k1(1)=f(3)
    k1(2)=f(4)
    k2(1)=f(5)
    k2(2)=f(6)
    k3(1)=f(7)
    k3(2)=f(8)
    k4(1)=f(9)
    k4(2)=f(10)
    k5(1)=f(11)
    k5(2)=f(12)
    k6(1)=f(13)
    k6(2)=f(14)
    
    k7(1)=f(15)
    k7(2)=f(16)
    k8(1)=f(17)
    k8(2)=f(18)
    k9(1)=f(19)
    k9(2)=f(20)
    k10(1)=f(21)
    k10(2)=f(22)
    k11(1)=f(23)
    k11(2)=f(24)
    k12(1)=f(25)
    k12(2)=f(26)
    
    
    !data Unity / 1., 0., 0., 1./ 
    !do i=0,10
    x = 0.!0.+i*0.1
    call MAT_Set_MatFunc_Sin_ANM(xappa, 1., x, 0., y, MatFunc_Sin1)
    call MAT_Set_MatFunc_Sin_ANM(xappa, sqrt(3.)/2., x, 0.5, y, MatFunc_Sin2)
    call MAT_Set_MatFunc_Sin_ANM(xappa, 0.5, x, sqrt(3.)/2., y, MatFunc_Sin3)
    call MAT_Set_MatFunc_Sin_ANM(xappa, 0., x, 1., y, MatFunc_Sin4)
    call MAT_Set_MatFunc_Sin_ANM(xappa, -0.5, x, sqrt(3.)/2., y, MatFunc_Sin5)
    call MAT_Set_MatFunc_Sin_ANM(xappa, -sqrt(3.)/2., x, 0.5, y, MatFunc_Sin6)
    call MAT_Set_MatFunc_Cos_ANM(xappa, 1., x, 0., y, MatFunc_Cos1)
    call MAT_Set_MatFunc_Cos_ANM(xappa, sqrt(3.)/2., x, 0.5, y, MatFunc_Cos2)
    call MAT_Set_MatFunc_Cos_ANM(xappa, 0.5, x, sqrt(3.)/2., y, MatFunc_Cos3)
    call MAT_Set_MatFunc_Cos_ANM(xappa, 0., x, 1., y, MatFunc_Cos4)
    call MAT_Set_MatFunc_Cos_ANM(xappa, -0.5, x, sqrt(3.)/2., y, MatFunc_Cos5)
    call MAT_Set_MatFunc_Cos_ANM(xappa, -sqrt(3.)/2., x, 0.5, y, MatFunc_Cos6)
    !print*,xappa*x
    
    print*, ''
    flux_result = k0 + matmul(MatFunc_Sin1, k1) + matmul(MatFunc_Cos1, k7) &
        +  matmul(MatFunc_Sin2, k2) + matmul(MatFunc_Cos2, k8) &
        +  matmul(MatFunc_Sin3, k3) + matmul(MatFunc_Cos3, k9) &
        +  matmul(MatFunc_Sin4, k4) + matmul(MatFunc_Cos4, k10) &
        +  matmul(MatFunc_Sin5, k5) + matmul(MatFunc_Cos5, k11) &
        +  matmul(MatFunc_Sin6, k6) + matmul(MatFunc_Cos6, k12)
    
    print*,x, 'Result:'
    print*,flux_result
    !enddo
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
end subroutine   

    
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
    
    
subroutine MAT_Set_MatFunc_ANM_1D(xappa, B2, Matr_Func)
    ! Calculation of the Matrix Functions included in the Sensitivity matrix  
    ! Last Update                 :: Vlad, 23.04.2021
    implicit none
    integer, parameter            :: NG = 2 ! Number of neutron groups
    
    ! Input variables
    real,    dimension(NG,NG)     :: xappa ! xappa matrix (alpha from the analitical expression)
    real,    dimension(NG,NG)     :: B2 ! xappa matrix (alpha from the analitical expression)
    
    !Local variables
    real                          :: denom,a ! temporary value
    real,    dimension(NG)        :: eigenv !eigenvalue array
    real,    dimension(NG)        :: eigenvb2 !eigenvalue array
    real,    dimension(169)         :: a1 ! temporary values
    real,    dimension(169)         :: a2 ! temporary values
    real,    dimension(NG)        :: lambda ! generic function 
                                           ! cos(x) if x > 0; cosh(x) if x < 0;
    real,    dimension(NG)        :: gamma ! generic function 
                                           ! sin(x) if x > 0; sinh(x) if x < 0;
    real,    dimension(NG)        :: f1    ! generic function 
                                           ! cos(x) if x > 0; cosh(x) if x < 0;
    real,    dimension(NG)        :: f2    ! generic function 
                                           ! sin(x) if x > 0; sinh(x) if x < 0;
    real,    dimension(NG)        :: f3    ! generic function 
                                           ! cos(x/sqrt(3)) if x > 0; cosh(x/sqrt(3)) if x < 0;
    real,    dimension(NG)        :: f4    ! generic function 
                                           ! cos(2*x/sqrt(3)) if x > 0; cosh(2*x/sqrt(3)) if x < 0;
    real,    dimension(NG)        :: f5    ! generic function 
                                           ! sin(x/sqrt(3)) if x > 0; sinh(x/sqrt(3)) if x < 0;
    real,    dimension(NG)        :: f6    ! generic function 
                                           ! sin(x/sqrt(3)) if x > 0; sinh(x/sqrt(3)) if x < 0;
    real,    dimension(NG,NG)     :: Unity ! Unity Marix
    integer                       :: n ! counter
    integer                       :: m ! counter
    integer                       :: i ! counter
    
    ! Output variables
    real,    dimension(NG,NG,169)   :: Matr_Func ! Array for Matrix Functions values storage
    
    data Unity / 1., 0., 0., 1./ 

    call MAT_Compute_Eigenv_2x2(xappa, eigenv)
    call MAT_Compute_Eigenv_2x2(b2, eigenvb2)
    

    ! setting up lambda and gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenvb2(n).LT.0) then
            f1(n) = cos(eigenv(n))
            f2(n) = sin(eigenv(n))
            f3(n) = cos(eigenv(n)/sqrt(3.))
            f4(n) = cos(2.*eigenv(n)/sqrt(3.))
            f5(n) = sin(eigenv(n)/sqrt(3.))
            f6(n) = sin(2.*eigenv(n)/sqrt(3.))
            !lambda(n) = cos(eigenv(n))
            !gamma(n) = sin(eigenv(n))
        else
            !lambda(n) = cosh(eigenv(n))
            !gamma(n) = sinh(eigenv(n))
            f1(n) = cosh(eigenv(n))
            f2(n) = sinh(eigenv(n))
            f3(n) = cosh(eigenv(n)/sqrt(3.))
            f4(n) = cosh(2.*eigenv(n)/sqrt(3.))
            f5(n) = sinh(eigenv(n)/sqrt(3.))
            f6(n) = sinh(2.*eigenv(n)/sqrt(3.))
        end if
    end do

    ! calculation of the matrix function 
    denom = Eigenv(1) - Eigenv(2)
    
    a1(1) = (F_01_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(2) = (F_01_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(3) = (F_01_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(4) = (F_01_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(5) = (F_01_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(6) = (F_01_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(7) = (F_01_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(8) = (F_01_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(9) = (F_01_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(10) = (F_01_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(11) = (F_01_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(12) = (F_01_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(13) = (F_01_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_01_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(14) = (F_02_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(15) = (F_02_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(16) = (F_02_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(17) = (F_02_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(18) = (F_02_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(19) = (F_02_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(20) = (F_02_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(21) = (F_02_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(22) = (F_02_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(23) = (F_02_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(24) = (F_02_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(25) = (F_02_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(26) = (F_02_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_02_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(27) = (F_03_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(28) = (F_03_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(29) = (F_03_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(30) = (F_03_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(31) = (F_03_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(32) = (F_03_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(33) = (F_03_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(34) = (F_03_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(35) = (F_03_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(36) = (F_03_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(37) = (F_03_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(38) = (F_03_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(39) = (F_03_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_03_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(40) = (F_04_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(41) = (F_04_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(42) = (F_04_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(43) = (F_04_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(44) = (F_04_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(45) = (F_04_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(46) = (F_04_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(47) = (F_04_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(48) = (F_04_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(49) = (F_04_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(50) = (F_04_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(51) = (F_04_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(52) = (F_04_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_04_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(53) = (F_05_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(54) = (F_05_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(55) = (F_05_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(56) = (F_05_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(57) = (F_05_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(58) = (F_05_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(59) = (F_05_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(60) = (F_05_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(61) = (F_05_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(62) = (F_05_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(63) = (F_05_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(64) = (F_05_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(65) = (F_05_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_05_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(66) = (F_06_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(67) = (F_06_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(68) = (F_06_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(69) = (F_06_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(70) = (F_06_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(71) = (F_06_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(72) = (F_06_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(73) = (F_06_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(74) = (F_06_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(75) = (F_06_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(76) = (F_06_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(77) = (F_06_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(78) = (F_06_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_06_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(79) = (F_07_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(80) = (F_07_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(81) = (F_07_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(82) = (F_07_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(83) = (F_07_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(84) = (F_07_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(85) = (F_07_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(86) = (F_07_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(87) = (F_07_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(88) = (F_07_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(89) = (F_07_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(90) = (F_07_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(91) = (F_07_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_07_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(92) = (F_08_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(93) = (F_08_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(94) = (F_08_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(95) = (F_08_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(96) = (F_08_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(97) = (F_08_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(98) = (F_08_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(99) = (F_08_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(100) = (F_08_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(101) = (F_08_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(102) = (F_08_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(103) = (F_08_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(104) = (F_08_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_08_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(105) = (F_09_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(106) = (F_09_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(107) = (F_09_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(108) = (F_09_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(109) = (F_09_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(110) = (F_09_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(111) = (F_09_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(112) = (F_09_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(113) = (F_09_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(114) = (F_09_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(115) = (F_09_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(116) = (F_09_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(117) = (F_09_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_09_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(118) = (F_10_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(119) = (F_10_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(120) = (F_10_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(121) = (F_10_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(122) = (F_10_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(123) = (F_10_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(124) = (F_10_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(125) = (F_10_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(126) = (F_10_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(127) = (F_10_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(128) = (F_10_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(129) = (F_10_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(130) = (F_10_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_10_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(131) = (F_11_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(132) = (F_11_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(133) = (F_11_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(134) = (F_11_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(135) = (F_11_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(136) = (F_11_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(137) = (F_11_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(138) = (F_11_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(139) = (F_11_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(140) = (F_11_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(141) = (F_11_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(142) = (F_11_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(143) = (F_11_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_11_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(144) = (F_12_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(145) = (F_12_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(146) = (F_12_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(147) = (F_12_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(148) = (F_12_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(149) = (F_12_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(150) = (F_12_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(151) = (F_12_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(152) = (F_12_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(153) = (F_12_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(154) = (F_12_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(155) = (F_12_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(156) = (F_12_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_12_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a1(157) = (F_13_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(158) = (F_13_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(159) = (F_13_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(160) = (F_13_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(161) = (F_13_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(162) = (F_13_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(163) = (F_13_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(164) = (F_13_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(165) = (F_13_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(166) = (F_13_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(167) = (F_13_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(168) = (F_13_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a1(169) = (F_13_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        - F_13_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    
      
    
                     
    a2(1) = (-Eigenv(2)*F_01_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(2) = (-Eigenv(2)*F_01_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(3) = (-Eigenv(2)*F_01_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(4) = (-Eigenv(2)*F_01_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(5) = (-Eigenv(2)*F_01_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(6) = (-Eigenv(2)*F_01_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(7) = (-Eigenv(2)*F_01_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(8) = (-Eigenv(2)*F_01_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(9) = (-Eigenv(2)*F_01_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(10) = (-Eigenv(2)*F_01_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(11) = (-Eigenv(2)*F_01_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(12) = (-Eigenv(2)*F_01_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(13) = (-Eigenv(2)*F_01_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_01_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(14) = (-Eigenv(2)*F_02_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(15) = (-Eigenv(2)*F_02_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(16) = (-Eigenv(2)*F_02_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(17) = (-Eigenv(2)*F_02_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(18) = (-Eigenv(2)*F_02_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(19) = (-Eigenv(2)*F_02_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(20) = (-Eigenv(2)*F_02_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(21) = (-Eigenv(2)*F_02_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(22) = (-Eigenv(2)*F_02_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(23) = (-Eigenv(2)*F_02_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(24) = (-Eigenv(2)*F_02_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(25) = (-Eigenv(2)*F_02_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(26) = (-Eigenv(2)*F_02_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_02_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(27) = (-Eigenv(2)*F_03_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(28) = (-Eigenv(2)*F_03_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(29) = (-Eigenv(2)*F_03_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(30) = (-Eigenv(2)*F_03_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(31) = (-Eigenv(2)*F_03_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(32) = (-Eigenv(2)*F_03_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(33) = (-Eigenv(2)*F_03_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(34) = (-Eigenv(2)*F_03_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(35) = (-Eigenv(2)*F_03_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(36) = (-Eigenv(2)*F_03_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(37) = (-Eigenv(2)*F_03_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(38) = (-Eigenv(2)*F_03_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(39) = (-Eigenv(2)*F_03_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_03_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(40) = (-Eigenv(2)*F_04_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(41) = (-Eigenv(2)*F_04_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(42) = (-Eigenv(2)*F_04_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(43) = (-Eigenv(2)*F_04_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(44) = (-Eigenv(2)*F_04_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(45) = (-Eigenv(2)*F_04_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(46) = (-Eigenv(2)*F_04_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(47) = (-Eigenv(2)*F_04_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(48) = (-Eigenv(2)*F_04_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(49) = (-Eigenv(2)*F_04_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(50) = (-Eigenv(2)*F_04_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(51) = (-Eigenv(2)*F_04_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(52) = (-Eigenv(2)*F_04_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_04_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(53) = (-Eigenv(2)*F_05_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(54) = (-Eigenv(2)*F_05_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(55) = (-Eigenv(2)*F_05_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(56) = (-Eigenv(2)*F_05_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(57) = (-Eigenv(2)*F_05_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(58) = (-Eigenv(2)*F_05_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(59) = (-Eigenv(2)*F_05_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(60) = (-Eigenv(2)*F_05_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(61) = (-Eigenv(2)*F_05_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(62) = (-Eigenv(2)*F_05_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(63) = (-Eigenv(2)*F_05_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(64) = (-Eigenv(2)*F_05_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(65) = (-Eigenv(2)*F_05_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_05_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(66) = (-Eigenv(2)*F_06_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(67) = (-Eigenv(2)*F_06_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(68) = (-Eigenv(2)*F_06_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(69) = (-Eigenv(2)*F_06_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(70) = (-Eigenv(2)*F_06_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(71) = (-Eigenv(2)*F_06_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(72) = (-Eigenv(2)*F_06_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(73) = (-Eigenv(2)*F_06_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(74) = (-Eigenv(2)*F_06_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(75) = (-Eigenv(2)*F_06_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(76) = (-Eigenv(2)*F_06_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(77) = (-Eigenv(2)*F_06_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(78) = (-Eigenv(2)*F_06_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_06_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(79) = (-Eigenv(2)*F_07_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(80) = (-Eigenv(2)*F_07_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(81) = (-Eigenv(2)*F_07_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(82) = (-Eigenv(2)*F_07_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(83) = (-Eigenv(2)*F_07_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(84) = (-Eigenv(2)*F_07_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(85) = (-Eigenv(2)*F_07_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(86) = (-Eigenv(2)*F_07_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(87) = (-Eigenv(2)*F_07_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(88) = (-Eigenv(2)*F_07_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(89) = (-Eigenv(2)*F_07_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(90) = (-Eigenv(2)*F_07_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(91) = (-Eigenv(2)*F_07_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_07_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(92) = (-Eigenv(2)*F_08_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(93) = (-Eigenv(2)*F_08_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(94) = (-Eigenv(2)*F_08_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(95) = (-Eigenv(2)*F_08_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(96) = (-Eigenv(2)*F_08_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(97) = (-Eigenv(2)*F_08_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(98) = (-Eigenv(2)*F_08_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(99) = (-Eigenv(2)*F_08_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(100) = (-Eigenv(2)*F_08_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(101) = (-Eigenv(2)*F_08_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(102) = (-Eigenv(2)*F_08_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(103) = (-Eigenv(2)*F_08_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(104) = (-Eigenv(2)*F_08_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_08_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(105) = (-Eigenv(2)*F_09_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(106) = (-Eigenv(2)*F_09_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(107) = (-Eigenv(2)*F_09_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(108) = (-Eigenv(2)*F_09_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(109) = (-Eigenv(2)*F_09_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(110) = (-Eigenv(2)*F_09_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(111) = (-Eigenv(2)*F_09_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(112) = (-Eigenv(2)*F_09_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(113) = (-Eigenv(2)*F_09_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(114) = (-Eigenv(2)*F_09_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(115) = (-Eigenv(2)*F_09_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(116) = (-Eigenv(2)*F_09_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(117) = (-Eigenv(2)*F_09_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_09_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(118) = (-Eigenv(2)*F_10_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(119) = (-Eigenv(2)*F_10_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(120) = (-Eigenv(2)*F_10_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(121) = (-Eigenv(2)*F_10_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(122) = (-Eigenv(2)*F_10_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(123) = (-Eigenv(2)*F_10_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(124) = (-Eigenv(2)*F_10_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(125) = (-Eigenv(2)*F_10_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(126) = (-Eigenv(2)*F_10_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(127) = (-Eigenv(2)*F_10_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(128) = (-Eigenv(2)*F_10_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(129) = (-Eigenv(2)*F_10_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(130) = (-Eigenv(2)*F_10_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_10_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(131) = (-Eigenv(2)*F_11_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(132) = (-Eigenv(2)*F_11_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(133) = (-Eigenv(2)*F_11_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(134) = (-Eigenv(2)*F_11_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(135) = (-Eigenv(2)*F_11_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(136) = (-Eigenv(2)*F_11_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(137) = (-Eigenv(2)*F_11_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(138) = (-Eigenv(2)*F_11_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(139) = (-Eigenv(2)*F_11_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(140) = (-Eigenv(2)*F_11_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(141) = (-Eigenv(2)*F_11_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(142) = (-Eigenv(2)*F_11_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(143) = (-Eigenv(2)*F_11_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_11_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(144) = (-Eigenv(2)*F_12_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(145) = (-Eigenv(2)*F_12_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(146) = (-Eigenv(2)*F_12_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(147) = (-Eigenv(2)*F_12_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(148) = (-Eigenv(2)*F_12_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(149) = (-Eigenv(2)*F_12_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(150) = (-Eigenv(2)*F_12_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(151) = (-Eigenv(2)*F_12_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(152) = (-Eigenv(2)*F_12_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(153) = (-Eigenv(2)*F_12_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(154) = (-Eigenv(2)*F_12_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(155) = (-Eigenv(2)*F_12_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(156) = (-Eigenv(2)*F_12_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_12_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
    a2(157) = (-Eigenv(2)*F_13_01(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_01(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(158) = (-Eigenv(2)*F_13_02(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_02(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(159) = (-Eigenv(2)*F_13_03(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_03(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(160) = (-Eigenv(2)*F_13_04(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_04(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(161) = (-Eigenv(2)*F_13_05(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_05(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(162) = (-Eigenv(2)*F_13_06(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_06(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(163) = (-Eigenv(2)*F_13_07(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_07(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(164) = (-Eigenv(2)*F_13_08(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_08(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(165) = (-Eigenv(2)*F_13_09(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_09(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom    
    a2(166) = (-Eigenv(2)*F_13_10(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_10(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(167) = (-Eigenv(2)*F_13_11(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_11(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    a2(168) = (-Eigenv(2)*F_13_12(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_12(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom   
    a2(169) = (-Eigenv(2)*F_13_13(f1(1), f2(1), f3(1), f4(1), f5(1), f6(1), eigenv(1)) &
        + Eigenv(1)*F_13_13(f1(2), f2(2), f3(2), f4(2), f5(2), f6(2), eigenv(2)))/denom
    
      
    
      
    do i = 1, 169
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
        
        !Start of the 13th row
        ! ||
        ! \/       
        real function F_13_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_01 = 1.

        end function
        
        real function F_13_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_02 = f2

        end function
        
        real function F_13_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_03 = f5

        end function
        
        real function F_13_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_04 = 0.

        end function
        
        real function F_13_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_05 = -f5

        end function
        
        real function F_13_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_06 = -f2
            
        end function
        
        real function F_13_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_07 = -f6

        end function
        
        real function F_13_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_08 = f1

        end function
        
        real function F_13_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_09 = f3

        end function
        
        real function F_13_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_10 = 1.

        end function
        
        real function F_13_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_11 = f3

        end function
        
        real function F_13_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_12 = f1

        end function
        
        real function F_13_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_13_13 = f4

        end function
    
    
    
        !Start of the 12th row
        ! ||
        ! \/       
        real function F_12_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_01 = 1.

        end function
        
        real function F_12_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_02 = 0.

        end function
        
        real function F_12_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_03 = -f5

        end function
        
        real function F_12_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_04 = -f2

        end function
        
        real function F_12_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_05 = -f6

        end function
        
        real function F_12_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_06 = -f2
            
        end function
        
        real function F_12_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_07 = -f5

        end function
        
        real function F_12_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_08 = 1.

        end function
        
        real function F_12_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_09 = f3

        end function
        
        real function F_12_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_10 = f1

        end function
        
        real function F_12_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_11 = f4

        end function
        
        real function F_12_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_12 = f1

        end function
        
        real function F_12_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_12_13 = f3

        end function
    
    
    
        !Start of the 11th row
        ! ||
        ! \/       
        real function F_11_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_01 = 1.

        end function
        
        real function F_11_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_02 = -f2

        end function
        
        real function F_11_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_03 = -f6

        end function
        
        real function F_11_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_04 = -f2

        end function
        
        real function F_11_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_05 = -f5

        end function
        
        real function F_11_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_06 = 0.
            
        end function
        
        real function F_11_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_07 = f5

        end function
        
        real function F_11_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_08 = f1

        end function
        
        real function F_11_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_09 = f4

        end function
        
        real function F_11_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_10 = f1

        end function
        
        real function F_11_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_11 = f3

        end function
        
        real function F_11_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_12 = 1.

        end function
        
        real function F_11_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_11_13 = f3

        end function
    
    
    
        !Start of the 10th row
        ! ||
        ! \/       
        real function F_10_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_01 = 1.

        end function
        
        real function F_10_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_02 = -f2

        end function
        
        real function F_10_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_03 = -f5

        end function
        
        real function F_10_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_04 = 0.

        end function
        
        real function F_10_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_05 = f5

        end function
        
        real function F_10_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_06 = f2
            
        end function
        
        real function F_10_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_07 = f6

        end function
        
        real function F_10_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_08 = f1

        end function
        
        real function F_10_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_09 = f3

        end function
        
        real function F_10_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_10 = 1.

        end function
        
        real function F_10_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_11 = f3

        end function
        
        real function F_10_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_12 = f1

        end function
        
        real function F_10_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_10_13 = f4

        end function
    
    
    
        !Start of the 9th row
        ! ||
        ! \/       
        real function F_09_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_01 = 1.

        end function
        
        real function F_09_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_02 = 0.

        end function
        
        real function F_09_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_03 = f5

        end function
        
        real function F_09_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_04 = f2

        end function
        
        real function F_09_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_05 = f6

        end function
        
        real function F_09_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_06 = f2
            
        end function
        
        real function F_09_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_07 = f5

        end function
        
        real function F_09_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_08 = 1.

        end function
        
        real function F_09_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_09 = f3

        end function
        
        real function F_09_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_10 = f1

        end function
        
        real function F_09_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_11 = f4

        end function
        
        real function F_09_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_12 = f1

        end function
        
        real function F_09_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_09_13 = f3

        end function
    
    
    
        !Start of the 8th row
        ! ||
        ! \/       
        real function F_08_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_01 = 1.

        end function
        
        real function F_08_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_02 = f2

        end function
        
        real function F_08_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_03 = f6

        end function
        
        real function F_08_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_04 = f2

        end function
        
        real function F_08_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_05 = f5

        end function
        
        real function F_08_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_06 = 0.
            
        end function
        
        real function F_08_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_07 = -f5

        end function
        
        real function F_08_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_08 = f1

        end function
        
        real function F_08_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_09 = f4

        end function
        
        real function F_08_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_10 = f1

        end function
        
        real function F_08_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_11 = f3

        end function
        
        real function F_08_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_12 = 1.

        end function
        
        real function F_08_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_08_13 = f3

        end function
    
    
    
        !Start of the 7th row
        ! ||
        ! \/       
        real function F_07_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_01 = 2./sqrt(3.)

        end function
        
        real function F_07_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_02 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_07_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_03 = 0.

        end function
        
        real function F_07_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_04 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_07_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_05 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_07_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_06 = -2.*f2/sqrt(3.)
            
        end function
        
        real function F_07_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_07 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_07_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_08 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_07_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_09 = 2.*f5/eigenv

        end function
        
        real function F_07_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_10 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_07_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_11 = 2.*(-f5 + f6)/eigenv

        end function
        
        real function F_07_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_12 = 2.*f1/sqrt(3.)

        end function
        
        real function F_07_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_07_13 = 2.*(-f5 + f6)/eigenv

        end function
    
    
    
        !Start of the 6th row
        ! ||
        ! \/       
        real function F_06_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_01 = 2./sqrt(3.)

        end function
        
        real function F_06_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_02 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_06_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_03 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_06_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_04 = -2.*f2/sqrt(3.)

        end function
        
        real function F_06_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_05 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_06_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_06 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)
            
        end function
        
        real function F_06_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_07 = 0.

        end function
        
        real function F_06_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_08 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_06_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_09 = 2.*(-f5 + f6)/(eigenv)

        end function
        
        real function F_06_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_10 = 2.*f1/sqrt(3.)

        end function
        
        real function F_06_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_11 = 2.*(-f5 + f6)/eigenv

        end function
        
        real function F_06_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_12 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_06_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_06_13 = 2.*f5/eigenv

        end function
    
    
    
        !Start of the 5th row
        ! ||
        ! \/       
        real function F_05_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_01 = 2./sqrt(3.)

        end function
        
        real function F_05_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_02 = -2.*f2/sqrt(3.)

        end function
        
        real function F_05_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_03 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_05_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_04 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_05_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_05 = 0.

        end function
        
        real function F_05_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_06 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_05_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_07 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_05_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_08 = 2.*f1/sqrt(3.)

        end function
        
        real function F_05_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_09 = 2.*(-f5 + f6)/(eigenv)

        end function
        
        real function F_05_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_10 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_05_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_11 = 2.*f5/eigenv

        end function
        
        real function F_05_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_12 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_05_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_05_13 = 2.*(-f5 + f6)/eigenv

        end function
    
    
        !Start of the 4th row
        ! ||
        ! \/       
        real function F_04_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_01 = 2./sqrt(3.)

        end function
        
        real function F_04_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_02 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_04_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_03 = 0.

        end function
        
        real function F_04_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_04 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_04_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_05 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_04_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_06 = 2.*f2/sqrt(3.)

        end function
        
        real function F_04_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_07 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_04_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_08 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_04_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_09 = 2.*f5/eigenv

        end function
        
        real function F_04_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_10 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_04_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_11 = 2.*(-f5 + f6)/eigenv

        end function
        
        real function F_04_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_12 = 2.*f1/sqrt(3.)

        end function
        
        real function F_04_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_04_13 = 2.*(-f5 + f6)/eigenv

        end function
    
    
    
        !Start of the 3rd row
        ! ||
        ! \/       
        real function F_03_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_01 = 2./sqrt(3.)

        end function
        
        real function F_03_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_02 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_03_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_03 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_03_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_04 = 2.*f2/sqrt(3.)

        end function
        
        real function F_03_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_05 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_03_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_06 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_03_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_07 = 0.

        end function
        
        real function F_03_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_08 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_03_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_09 = 2.*(-f5 + f6)/(eigenv)

        end function
        
        real function F_03_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_10 = 2.*f1/sqrt(3.)

        end function
        
        real function F_03_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_11 = 2.*(-f5 + f6)/(eigenv)

        end function
        
        real function F_03_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_12 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_03_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_03_13 = 2.*f5/eigenv

        end function
    
    
        
        !Start of the 2nd row
        ! ||
        ! \/
        
        real function F_02_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_01 = 2./sqrt(3.)

        end function
        
        real function F_02_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_02 = 2.*f2/sqrt(3.)

        end function
        
        real function F_02_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_03 = 2.*(f3 - f4)/eigenv

        end function
        
        real function F_02_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_04 = -2.*(f1 - 1)/(sqrt(3.)*eigenv)

        end function
        
        real function F_02_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_05 = 0.

        end function
        
        real function F_02_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_06 = 2.*(f1 - 1)/(sqrt(3.)*eigenv)
            
        end function
        
        real function F_02_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_07 = 2.*(-f3 + f4)/eigenv

        end function
        
        real function F_02_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_08 = 2.*f1/sqrt(3.)

        end function
        
        real function F_02_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_09 = 2.*(-f5 + f6)/(eigenv)

        end function
        
        real function F_02_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_10 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_02_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_11 = 2.*f5/eigenv

        end function
        
        real function F_02_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_12 = 2.*f2/(sqrt(3.)*eigenv)

        end function
        
        real function F_02_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_02_13 = 2.*(-f5 + f6)/eigenv

        end function
    
    
        !Start of the 1st row
        ! ||
        ! \/
        real function F_01_01(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_01 = 1.

        end function
        
        real function F_01_02(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_02 = 0.

        end function
        
        real function F_01_03(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_03 = 0.

        end function
        
        real function F_01_04(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_04 = 0.

        end function
        
        real function F_01_05(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_05 = 0.

        end function
        
        real function F_01_06(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_06 = 0.

        end function
        
        real function F_01_07(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_07 = 0.

        end function
        
        real function F_01_08(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_08 = 2.*(1. - f1 + eigenv*f2)/(3.*eigenv*eigenv)

        end function
        
        real function F_01_09(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_09 = 2.*(f3 - f4)/(eigenv*eigenv)

        end function
        
        real function F_01_10(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_10 = 2.*(1. - f1 + eigenv*f2)/(3.*eigenv*eigenv)

        end function
        
        real function F_01_11(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_11 = 2.*(f3 - f4)/(eigenv*eigenv)

        end function
        
        real function F_01_12(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_12 = 2.*(1. - f1 + eigenv*f2)/(3.*eigenv*eigenv)

        end function
        
        real function F_01_13(f1, f2, f3, f4, f5, f6, eigenv)
            real f1, f2, f3, f4, f5, f6, eigenv

            F_01_13 = 2.*(f3 - f4)/(eigenv*eigenv)

        end function

    end subroutine
    

subroutine MAT_Set_MatFunc_Sin_ANM(xappa, kx, x, ky, y, Matr_Func_Sin)
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
    real,    dimension(NG,NG)     :: Matr_Func_Sin ! Matrix Function SQRT(B2) = xappa
    real :: kx, x, ky, y
    data Unity / 1., 0., 0., 1./ 
 
    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            gamma(n) = sin((kx*x+ky*y)*eigenv(n))
        else
            gamma(n) = sinh((kx*x+ky*y)*eigenv(n))
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


subroutine MAT_Set_MatFunc_Cos_ANM(xappa, kx, x, ky, y, Matr_Func_Cos)
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
    real,    dimension(NG,NG)     :: Matr_Func_Cos ! Matrix Function SQRT(B2) = xappa
    real :: kx, x, ky, y
    data Unity / 1., 0., 0., 1./ 
 
    call MAT_Compute_Eigenv_2x2(xappa, eigenv)

    ! setting up gamma function values depending on the eigenvalue
    do n = 1, NG
        if(Eigenv(n).LT.0) then
            gamma(n) = cos((kx*x+ky*y)*eigenv(n))
        else
            gamma(n) = cosh((kx*x+ky*y)*eigenv(n))
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
      print*, 'a1', a1
      print*, 'b1', b1

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
     
    
          SUBROUTINE INVERT1(A,A1,N)
!         A1=INV(A)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION A(N,N),A1(N,N)
      A1=A
      EPS=1E-20
      DO 21 I=1,N
      Q=A1(I,I)
      Q1=Q
      IF(Q.LT.0.) Q1=-Q1
      IF(Q1.GE.EPS) GO TO 3
	PRINT*,'matritsa A virojdennaia'  
       READ*; STOP
    3 A1(I,I)=1.
      DO 22 K=1,N
   22 A1(I,K)=A1(I,K)/Q
      DO 25 J=1,N
      IF(I.EQ.J) GO TO 25
      Q=A1(J,I)
      A1(J,I)=0.
      DO 26 K=1,N
      A1(J,K)=A1(J,K)-A1(I,K)*Q
   26 CONTINUE
   25 CONTINUE
   21 CONTINUE
      RETURN
      END