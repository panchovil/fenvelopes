program main
    use dtypes, only: envelope
    use constants, only: pr
    use legacy_ar_models, only: nc, z
    use legacy_thermo_properties, only: TERMO
    use io_nml, only: read_system
    
    implicit none
    character(len=*), parameter :: infile="test/input_prueba_jaco.nml"
    real(pr) :: pt_bub_t0 = 180
    real(pr) :: pt_dew_t0 = 180

    type(envelope) :: pt_bub, pt_dew, pt_hpl

    call read_system(infile)
    
    
    
    call pt_envelopes



contains
    subroutine pt_envelopes
        use envelopes, only: envelope2, max_points, k_wilson_bubble, &
                            max_points, p_wilson, k_wilson, find_hpl, get_case    
        real(pr), allocatable :: tv(:) ! Temperatures [K]
        real(pr), allocatable :: pv(:) ! Pressures [bar]
        real(pr), allocatable :: dv(:) ! Densities [mol/L]

        real(pr) :: tcri(4)            ! Critical points temperatures
        real(pr) :: pcri(4)            ! Critical points pressures
        real(pr) :: dcri(4)            ! Critical points densities

        real(pr) :: t, p               ! Temperature and pressure
        real(pr), allocatable :: k(:)  ! K factors
        integer :: n_points, icri(4), ncri, i

        character(len=:), allocatable :: pt_case
  
        integer :: n

        ! Capilar Preasure variables
        real(pr) :: IFT, Pcap, r_poro, ang_cont
        real(pr) :: Pliq, Pvap
        real(pr) :: Vx, Vy

        allocate (tv(max_points), pv(max_points), dv(max_points))
        allocate (k(size(z)))
  
        print *, "Prueba"
        !print *, "P y T: ", p,t

        
        ! ========================================================================
        !  Bubble envel
        ! ------------------------------------------------------------------------
        call k_wilson_bubble(z, t_0=pt_bub_t0, p_end=0.5_pr, t=t, p=p, k=k)
        !print *, "log P y T: ", log(p),log(t)
        !print *, "z: ", z
        !print *, "log k", log(k)
        
        print*, "-----------------------BULK-------------------------------"
        
        !iniciador bulk y Laplace para sistema capilar
        call iniciador_bulk(1, nc, z, T, P, k, Vx, Vy)
        
        call Laplace(r_poro, ang_cont, Vx, Vy, k, IFT, Pcap)!radio de poro en metros
        Pvap = P
        Pliq = Pvap - Pcap

        print*, "T inicial", T
        print*, "Vx inicial", Vx
        print*, "Vy inicial", Vy
        print*, "k inicial", k
        print*, "P capilar inicial", Pcap   
        print*, "P liq inicial", Pliq   
        print*, "P vap inicial", Pvap  


        

        !call envelope2( &
        !   1, nc, z, T, P, k, &
        !   n_points, Tv, Pv, Dv, ncri, icri, Tcri, Pcri, Dcri, &
        !   pt_bub &
        !)    
        print*, "-----------------------Nano-------------------------------"
        !call Prueba_nano(1,nc,z,T,Pliq,Pvap,k)
        !call pepe

    end subroutine
    subroutine iniciador_bulk(ichoice, n, z, T, P, KFACT, Vx, Vy) ! This will probably always exist


        use envelopes, only: fix_delx
        use linalg, only: solve_system

         ! number of compounds in the system and starting point type
        integer, intent(in) :: n, ichoice

        ! estimated T and P for first point (then used for every point)
  
        ! Maximun pressure
        real(pr) :: maxP
  
        ! estimated K factors for first point (then used for every point)
        real(pr), intent(in out) :: KFACT(n), T, P
  
        ! composition of the system
        real(pr), intent(in) :: z(n)
  

         ! Intermediate variables during calculation process
        real(pr), dimension(n) :: y
        real(pr), dimension(n + 2) :: X, Xold, Xold2, delX, bd, F, dXdS
        real(pr), dimension(n + 2, n + 2) :: JAC, AJ
        real(pr), intent(out) :: Vx, Vy
        logical :: run, passingcri, minT, minmaxT
  
        character(len=:), allocatable :: incipient_phase    
        
        integer :: i
        integer :: iy, ix ! Vapor or liquid selectors
  
        ! Specification value, delta and index
        real(pr) :: S, delS
        integer :: ns
  
        real(pr) :: Told2, Told
        real(pr) :: frac
  
        ! Netwon method
        integer :: iter ! Iteration
        integer :: max_iter   



        minT = .false.
        minmaxT = .false.
        passingcri = .false.
        Told2 = 0.0
        Told = 10.0
        maxP = 0.d0
  
        !-----------------------------------------------------------
        ! Continuation method for tracing the envelope starts here
        run = .true.
        i = 0
        JAC(n + 1, :) = 0.d0
        X(:n) = log(KFACT)
        X(n + 1) = log(T)
        X(n + 2) = log(P)
        iy = 1
        ix = 1
  
        select case(ichoice)
        case (1)
           incipient_phase = "vapor"
           iy = -1
        case (2)
           incipient_phase = "liquid"
           ix = -1
        case (3)
           incipient_phase = "2ndliquid"
        end select

        if (ichoice <= 2) then
            ! low T bub (1) or dew (2)
            ! x will be vapor phase during the first part, 
            ! and liquid after a critical point is crossed
            if (ichoice == 1) iy = -1
            if (ichoice == 2) ix = -1
            ns = n + 1
            S = log(T)
            delS = 0.005
   
            ! Wilson estimate for vapor (or liquid) composition
            y = KFACT*z
        else
            ! (ichoice==3) high P L-L sat
            ! PmaxDewC = maxval(PdewC(1:ilastDewC))
            ns = n + 2
            S = log(P)
            delS = -0.05
            y = kfact * z
            ! y = 0.d0
            ! y(n) = 1.d0
        end if
   
        Xold = 0.d0

   
        i=0
        do while (run)
            i=1
            ! Newton starts here
            delX = 1.0
            iter = 0
            max_iter = 500
   
            do while (maxval(abs(delX)) > 1.d-9 .and. iter <= max_iter)
               ! Solve point with full Newton method
               call Fbulk(incipient_phase, z, y, X, S, ns, F, JAC, Vx, Vy)
   
               iter = iter + 1
   
               bd = -F
               AJ = JAC
               delX = solve_system(AJ, bd)
               call fix_delX(i, iter, 3, 10.0_pr, 0.08_pr, delX)
   
               X = X + delX
   
               if (.not. passingcri .and. i /= 1 &
                   .and. iter > 10 &
                   .and. maxval(abs(delX)) > 0.001) then 
                  ! Too many iterations --> Reduce step to new point
   
                  delS = delS*2.0/4.0
                  S = S - delS
                  X = Xold + dXdS*delS
               end if
   
               KFACT = exp(X(:n))
               y = z*KFACT
               T = exp(X(n + 1))
               P = exp(X(n + 2))
               

            end do
            ! Point converged (unless it jumped out because of high number of iterations)
            if (iter > max_iter) run = .false.
            if (P > maxP) maxP = P
   
            !print *, "log P y log T: ", p,t
            !print *, "y: ",y
            !print *, "log kfact", KFACT
            !print *, "Vx, Vy: ", Vx, Vy
            !print *, "iteraciones finales : ", iter

            !Tv(i) = T
            !Pv(i) = P
            run = .false.
   
        end do
        print*, "iteraciones bulk", iter
        print*, "convergencia bulk", F        
 
    end subroutine iniciador_bulk
    subroutine Fbulk(incipient, z, y, X, S, ns, F, dF, Vx, Vy)
        character(len=*), intent(in) :: incipient
        real(pr), intent(in) :: z(:)
        real(pr), intent(in) :: X(nc + 2)
        real(pr), intent(in) :: y(nc)
        real(pr), intent(in) :: S
        integer, intent(in) :: ns
  
        real(pr), intent(out) :: F(nc + 2)
        real(pr), intent(out) :: dF(nc + 2, nc + 2)
        real(pr), intent(out) :: Vx, Vy
  
        real(pr) :: lnfug_x(nc), lnfug_y(nc)
        real(pr) :: dlnphi_dt_x(nc), dlnphi_dt_y(nc)
        real(pr) :: dlnphi_dp_x(nc), dlnphi_dp_y(nc)
        real(pr) :: dlnphi_dn_x(nc, nc), dlnphi_dn_y(nc, nc)
  
        real(pr) :: T, P
  
        integer :: ix, iy, n, j
  
        n = size(z)
        F = 0
        dF = 0
  
        T = exp(X(n+1))
        P = exp(X(n+2))
  
        select case(incipient)
        case ("liquid")
           ix = -1
           iy = 1
        case ("vapor")
           ix = 1
           iy = -1
        case ("2ndliquid")
           ix = 1
           iy = 1
        case default
           ix = 0
           iy = 0
        end select
  
        call TERMO(n, iy, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y)
        call TERMO(n, ix, 2, T, P, z, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x)
  
        F(:n) = X(:n) + lnfug_y - lnfug_x  ! X(:n) are LOG_K
        F(n + 1) = sum(y - z)
        F(n + 2) = X(ns) - S
  
        ! Jacobian Matrix
        do j=1,n
           df(:n, j) = dlnphi_dn_y(:, j) * y(j)
           df(j, j) = dF(j, j) + 1
        end do
  
        df(:n, n + 1) = T * (dlnphi_dt_y - dlnphi_dt_x)
        df(:n, n + 2) = P * (dlnphi_dp_y - dlnphi_dp_x)
  
        df(n + 1, :n) = y
  
        df(n + 2, :) = 0
        df(n + 2, ns) = 1
    end subroutine Fbulk
    subroutine Laplace(r_poro,ang_cont,Vx,Vy,k,IFT,Pcap,Parachor)
        use legacy_ar_models, only: tc,pc,w,nc,z
        real(pr), optional, intent(out) :: Parachor(nc)
        real(pr), intent(out) :: IFT, Pcap
        real(pr), intent(in) :: Vx, Vy, k(nc)
        real(pr), intent(out) :: r_poro, ang_cont
        integer :: i
        real(pr) :: y(nc)
        real(pr) :: Par(nc)

        

        r_poro=0.0000001 !radio cualquiera de 100 nm
        ang_cont=0.0 !angulo cualquiera de completamente mojado
        ang_cont=ang_cont*3.14/180.0 !la variable esta en º y se necesita en radianes
        y = z*k

        !Par cm^3/mol * (mN/m)^1/4, P bar, T K
        do i=1,nc
            Par(i)= 40.1684*(0.151-0.0464*w(i))*(Tc(i)**(13.0/12.0))/(Pc(i)**(5.0/6.0))  ! https://doi.org/10.1002/cjce.5450750617 // eq(12)
        end do
        IFT=0.0
        !print*, Par

        !IFT (mN/m)^1/4
        do i=1,nc
            IFT=IFT+((Par(i)/1000.0)*(z(i)/Vx-y(i)/Vy))
        end do

        !Pcap bar
        Pcap = (0.00000001*2.0*(IFT**4)*cos(ang_cont))/r_poro !E=3.6
        
        !if(present(Parachor)) then
        !    do i=1,nc
        !        Parachor(i)=Par(i)
        !    end do
        !end if
        if(present(Parachor)) Parachor=Par

    end subroutine Laplace
    subroutine Prueba_nano(ichoice, n, z, T, Pliq, Pvap, KFACT) ! This will probably always exist
        

        use envelopes, only: fix_delx
        use linalg, only: solve_system

         ! number of compounds in the system and starting point type
        integer, intent(in) :: n, ichoice

        ! estimated T and P for first point (then used for every point)
        real(pr), intent(in out) :: T, Pliq, Pvap
  
        ! estimated K factors for first point (then used for every point)
        real(pr), intent(in out) :: KFACT(n)
  
        ! composition of the system
        real(pr), intent(in) :: z(n)

         ! Intermediate variables during calculation process
        real(pr), dimension(n) :: y
        real(pr), dimension(n + 3) :: X, Xold, Xold2, delX, bd, F, dXdS
        real(pr), dimension(n + 3, n + 3) :: JAC, AJ
        real(pr) :: Vx, Vy
        logical :: run, passingcri, minT, minmaxT
  
        character(len=:), allocatable :: incipient_phase    
        
        integer :: i
        integer :: iy, ix ! Vapor or liquid selectors
  
        ! Specification value, delta and index
        real(pr) :: S, delS
        integer :: ns
  
        real(pr) :: Told2, Told
        real(pr) :: frac
  
        ! Netwon method
        integer :: iter ! Iteration
        integer :: max_iter   

        ! Capilar Preasure variables
        real(pr) :: IFT, Pcap, r_poro, ang_cont, Par(nc) 


  
        minT = .false.
        minmaxT = .false.
        passingcri = .false.
        Told2 = 0.0
        Told = 10.0
  
        !-----------------------------------------------------------
        ! Continuation method for tracing the envelope starts here
        run = .true.
        i = 0
        JAC(n + 1, :) = 0.d0
        X(:n) = log(KFACT)
        X(n + 1) = log(T)
        X(n + 2) = Pliq
        X(n + 3) = Pvap
        iy = 1
        ix = 1
  
        select case(ichoice)
        case (1)
           incipient_phase = "vapor"
           iy = -1
        case (2)
           incipient_phase = "liquid"
           ix = -1
        case (3)
           incipient_phase = "2ndliquid"
        end select

        if (ichoice <= 2) then
            ! low T bub (1) or dew (2)
            ! x will be vapor phase during the first part, 
            ! and liquid after a critical point is crossed
            if (ichoice == 1) iy = -1
            if (ichoice == 2) ix = -1
            ns = n + 1
            S = log(T)
            delS = 0.005
   
            ! Wilson estimate for vapor (or liquid) composition
            y = KFACT*z
        else
            ! (ichoice==3) high P L-L sat
            ! PmaxDewC = maxval(PdewC(1:ilastDewC))
            !ns = n + 2
            !S = log(P)!hay que cambiar esto
            !delS = -0.05
            !y = kfact * z
            ! y = 0.d0
            ! y(n) = 1.d0
        end if
   
        Xold = 0.d0

   
        i=0
        do while (run)
            i=1
            ! Newton starts here
            delX = 1.0
            iter = 0
            max_iter = 500
   
            do while (maxval(abs(delX)) > 1.d-9 .and. iter <= max_iter)
               ! Solve point with full Newton method
               call Fnano(incipient_phase, z, y, X, S, ns, F, JAC, Pcapilar=Pcap, Vliq=Vx, Vvap=Vy)
   
               iter = iter + 1
   
               bd = -F
               AJ = JAC
               delX = solve_system(AJ, bd)
               call fix_delX(i, iter, 3, 10.0_pr, 0.08_pr, delX)
               !print*, "del x",delX
               X = X + delX
                !print*, "x", x
               if (.not. passingcri .and. i /= 1 &
                   .and. iter > 10 &
                   .and. maxval(abs(delX)) > 0.001) then 
                  ! Too many iterations --> Reduce step to new point
   
                  delS = delS*2.0/4.0
                  S = S - delS
                  X = Xold + dXdS*delS
               end if
   
               KFACT = exp(X(:n))
               y = z*KFACT
               T = exp(X(n + 1))
               Pliq = X(n + 2)
               Pvap = X(n + 3)
               
               !print*, "convergencia nano", F
               !print*, "Pcap", Pcap

            end do
            ! Point converged (unless it jumped out because of high number of iterations)
            if (iter > max_iter) run = .false.
   
            !print *, "log P y log T: ", p,t
            !print *, "y: ",y
            !print *, "log kfact", KFACT
            !print *, "Vx, Vy: ", Vx, Vy


            !Tv(i) = T
            !Pv(i) = P
            run = .false.
   
        end do
        print *, "iteraciones finales nano: ", iter
        print*, "convergencia nano", F
        !print*, "Vx final", Vx 
        !print*, "T final", Vy
        print*, "T final", T 
        !print*, "k final", KFACT  
        !print*, "P capilar final", Pcap    
        !print*, "P liq final", Pliq   
        !print*, "P vap final", Pvap 
    end subroutine Prueba_nano

    subroutine Fnano(incipient, z, y, X, S, ns, F, dF, Pcapilar, Vliq, Vvap)
        character(len=*), intent(in) :: incipient
        real(pr), intent(in) :: z(:)
        real(pr), intent(in) :: X(nc + 3)
        real(pr), intent(in) :: y(nc)
        real(pr), intent(in) :: S
        integer, intent(in) :: ns
  
        real(pr), intent(out) :: F(nc + 3)
        real(pr), intent(out) :: dF(nc + 3, nc + 3)
        real(pr) :: Vx, Vy
  
        real(pr) :: lnfug_x(nc), lnfug_y(nc)
        real(pr) :: dlnphi_dt_x(nc), dlnphi_dt_y(nc)
        real(pr) :: dlnphi_dp_x(nc), dlnphi_dp_y(nc)
        real(pr) :: dlnphi_dn_x(nc, nc), dlnphi_dn_y(nc, nc)
  
        real(pr) :: T, Pliq, Pvap
  
        integer :: ix, iy, n, j, i
        
        ! Capilar Preasure variables
        real(pr), optional, intent(out) ::  Pcapilar, Vliq, Vvap
        real(pr) :: IFT, r_poro, ang_cont, Pcap
        real(pr) :: Par(nc)
        real(pr) :: dPvap_dV_y, dPliq_dV_x, dPvap_dT_y, dPliq_dT_x, dV_x_dT, dV_y_dT
        real(pr) :: dPvap_dn_y(nc), dPliq_dn_x(nc), dV_y_dn_y(nc)
        
        !Capilar jacobian intermediate variables
        real(pr) :: var_dFn2, var_dFn2_dK, var_dFn2_dT, var_dFn2_dPliq, var_dFn2_dPvap

        n = size(z)
        F = 0
        dF = 0
  
        T = exp(X(n+1))
        Pliq = X(n+2)
        Pvap = X(n+3)
  
        select case(incipient)
        case ("liquid")
           ix = -1
           iy = 1
        case ("vapor")
           ix = 1
           iy = -1
        case ("2ndliquid")
           ix = 1
           iy = 1
        case default
           ix = 0
           iy = 0
        end select
  
        call TERMO(n, iy, 4, T, Pvap, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y,&
        dp_dv=dPvap_dV_y, dp_dt=dPvap_dT_y, dp_dn=dPvap_dn_y)! dPvap_dV_y, dPvap_dT_y, dPvap_dn_y
        call TERMO(n, ix, 2, T, Pliq, z, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x,&
        dp_dv=dPliq_dV_x, dp_dt=dPliq_dT_x, dp_dn=dPliq_dn_x)! dPliq_dV_x, dPliq_dT_x, dPliq_dn_x)
        call Laplace(r_poro,ang_cont,Vx,Vy,exp(X(:n)),IFT,Pcap,Parachor=Par)

        F(:n) = X(:n) + (lnfug_y+log(Pvap)) - (lnfug_x+log(Pliq))  ! X(:n) are LOG_K
        F(n + 1) = sum(y - z)
        F(n + 2) = Pliq - Pvap + Pcap
        F(n + 3) = X(ns) - S
  
        ! Jacobian Matrix
        do j=1,n
           df(:n, j) = dlnphi_dn_y(:, j) * y(j) + 1.0
           !df(j, j) = dF(j, j) + 1.0
        end do
  
        df(:n, n + 1) = T * (dlnphi_dt_y - dlnphi_dt_x)
        df(:n, n + 2) = -(1/Pliq) - (dlnphi_dp_x)
        df(:n, n + 3) = (1/Pvap) + (dlnphi_dp_y)
  
        df(n + 1, :n) = y
        
        var_dFn2 = (8.0*cos(ang_cont)/r_poro)*(IFT**3.0)

        do i=1,n
            dV_y_dn_y(i) = -(dPvap_dn_y(i)/dPvap_dV_y)
            var_dFn2_dK = var_dFn2_dK+(Par(i)/1000.0*((y(i)*dV_y_dn_y(i)/(Vy**2.0))-(1.0/Vy)))
        end do
        do i=1,n
            df(n + 2, i) = 0.00000001*y(i)*var_dFn2*var_dFn2_dK
        end do

        dV_x_dT = -(dPliq_dT_x/dPliq_dV_x)
        dV_y_dT = -(dPvap_dT_y/dPvap_dV_y)
        do i=1,n
            var_dFn2_dT = var_dFn2_dT+(Par(i)/1000.0*((y(i)*dV_y_dT/(Vy**2.0))-(z(i)*dV_x_dT/(Vx**2.0))))
        end do
        df(n + 2, n + 1) = 0.00000001*var_dFn2*T*var_dFn2_dT
        do i=1,n
            var_dFn2_dPliq = var_dFn2_dPliq + (-Par(i)/1000.0*z(i))
            var_dFn2_dPvap = var_dFn2_dPvap + (Par(i)/1000.0*y(i))
        end do
        df(n + 2, n + 2) = 1.0 + 0.00000001*var_dFn2*var_dFn2_dPliq*(1.0/(dPliq_dV_x*(Vx**2.0))) 
        df(n + 2, n + 3) = - 1.0 + 0.00000001*var_dFn2*var_dFn2_dPvap*(1.0/(dPvap_dV_y*(Vy**2.0)))

        df(n + 3, :) = 0.0
        df(n + 3, ns) = 1.0

        if (present(Pcapilar)) Pcapilar=Pcap
        if (present(Vliq)) Vliq=Vx
        If (present(Vvap)) Vvap=Vy

    end subroutine Fnano
    subroutine pepe
        integer :: j,n,p,i
        real :: df(4,4),asd(4,4)
        n=4
        ! Jacobian Matrix
        p=0
        do i=1,n
            do j=1,n
                p=p+1
                asd(i,j)=p
            end do
        end do
        do j=1,n
            df(:n, j) = asd(:, j)+ 0.1
            !df(j, j) = dF(j, j) + 0.1
        end do
        print*, asd(1,1),asd(1,2),asd(1,3),asd(1,4)
        print*, asd(2,1),asd(2,2),asd(2,3),asd(2,4)
        print*, asd(3,1),asd(3,2),asd(3,3),asd(3,4)
        print*, asd(4,1),asd(4,2),asd(4,3),asd(4,4)

        print*, "------------------------------------"

        print*, df(1,1),df(1,2),df(1,3),df(1,4)
        print*, df(2,1),df(2,2),df(2,3),df(2,4)
        print*, df(3,1),df(3,2),df(3,3),df(3,4)
        print*, df(4,1),df(4,2),df(4,3),df(4,4)
        !print*, df
    end subroutine pepe
  

end program main