module inj_envelopes
   !! Module to calculate Px phase envelopes
   use constants, only: pr, R
   use dtypes, only: envelope, critical_point
   use linalg, only: solve_system, interpol, full_newton
   !use progress_bar_module, only: progress_bar

   implicit none

   type, extends(envelope) :: injelope
      real(pr), allocatable :: alpha(:) !! Ammount of injected fluid
      real(pr), allocatable :: z_inj(:) !! Injected fluid composition
      real(pr), allocatable :: z_mix(:, :) !! Composition at each step
   end type

   ! ===========================================================================
   !  Parameters
   ! ---------------------------------------------------------------------------
   integer :: env_number = 0 !! Number of calculated envelope
   integer :: max_iters = 100 !! Maximum number of iterations for a newton step
   real(pr) :: solve_tol =1e-5 !! Newton solver tolerance
   integer, parameter :: max_points = 1000 !! Maximum number of points for each envelope
   real(pr), allocatable :: z_0(:) !! Original fluid composition
   real(pr), allocatable :: z_injection(:) !! Injection fluid composition
   real(pr) :: T !! Temperature of injection
   character(len=10) :: injection_case !! Kind of injection displace|dilute

   ! Two-phase settings
   real(pr) :: del_S_multiplier = 1.5_pr
   real(pr) :: max_dp=100.0_pr
   real(pr) :: max_dalpha=0.05_pr
   real(pr) :: critical_multiplier = 1.5_pr
   
   ! Three phase parameters
   real(pr) :: del_S_multiplier_three_phase = 1.7_pr
   real(pr) :: critical_fact = 3.0_pr
   ! ===========================================================================
contains
   subroutine from_nml(filepath)
      use legacy_ar_models, only: nc
      character(len=*), intent(in) :: filepath
      integer :: funit

      namelist /nml_px/ T, z_0, z_injection, injection_case, &
                        del_S_multiplier, max_dp, max_dalpha, critical_multiplier,&
                        del_S_multiplier_three_phase, critical_fact, &
                        max_iters, solve_tol

      allocate (z_0(nc), z_injection(nc))

      open (newunit=funit, file=filepath)
         read (funit, nml=nml_px)
      close (funit)

      z_injection = z_injection/sum(z_injection)
   end subroutine

   ! ===========================================================================
   ! Two-phases
   ! ---------------------------------------------------------------------------
   subroutine F_injection(X, ns, S, F, dF)
      !! Function to solve at each point of the phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnK_i, lnP, \alpha] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [lnK_i - ln \phi_i(y, P, T) + ln \phi_i(z, P, T),
      !!         \sum_{i=1}^N, X_{ns} - S] \)
      !!
      !! The injection can be considered as two kinds of injection:
      !! - Displacement: \( z = \alpha z_i + (1-\alpha) z_0 \)
      !! - Addition:  \( z = \frac{\alpha z_i + (1-\alpha) z_0}{\sum_{i=1}^N \alpha z_i + (1-\alpha) z_0} \)
      !!
      use iso_fortran_env, only: error_unit
      use legacy_thermo_properties, only: TERMO
      real(pr), intent(in)  :: X(:) !! Vector of variables
      integer, intent(in)   :: ns !! Number of specification
      real(pr), intent(in)  :: S !! Specification value
      real(pr), intent(out) :: F(size(X)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(x), size(X)) !! Jacobian matrix

      ! X variables
      real(pr) :: K(size(X) - 2)
      real(pr) :: alpha
      real(pr) :: P

      ! Main phase variables
      real(pr) :: Vz
      real(pr), dimension(size(X) - 2) :: z, lnfug_z, dlnphi_dt_z, dlnphi_dp_z
      real(pr), dimension(size(X) - 2, size(X) - 2) :: dlnphi_dn_z

      ! Incipient phase variables
      real(pr) :: Vy
      real(pr), dimension(size(X) - 2) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension(size(X) - 2, size(X) - 2) :: dlnphi_dn_y

      ! Derivative of z wrt alpha
      real(pr) :: dzda(size(X) - 2)

      integer :: i, j, n

      n = size(X) - 2
      K = exp(X(1:n))
      P = exp(X(n + 1))
      alpha = X(n + 2)

      call get_z(alpha, z, dzda)

      ! if (any(z < 0)) z = 0

      y = K*z

      call TERMO( &
         n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y &
         )
      call TERMO( &
         n, 0, 4, T, P, z, Vz, lnfug_z, dlnphi_dp_z, dlnphi_dt_z, dlnphi_dn_z &
         )

      F(1:n) = X(:n) + lnfug_y - lnfug_z
      F(n + 1) = sum(y - z)
      F(n + 2) = X(ns) - S

      df = 0

      do i = 1, n
         do j = 1, n
            df(i, j) = y(j)*dlnphi_dn_y(i, j)
         end do
         df(i, i) = df(i, i) + 1
         df(i, n + 2) = sum(K*dlnphi_dn_y(i, :)*dzda - dlnphi_dn_z(i, :)*dzda)
      end do

      df(:n, n + 1) = P*(dlnphi_dp_y - dlnphi_dp_z)
      df(n + 1, :n) = y
      df(n + 1, n + 2) = sum(dzda*(K - 1))

      df(n + 2, :) = 0
      df(n + 2, ns) = 1
   end subroutine

   subroutine injection_envelope(X0, spec_number, del_S0, envels)
      !! Subroutine to calculate Px phase envelopes via continuation method
      use constants, only: ouput_path
      use ftools__io
      real(pr), intent(in) :: X0(:) !! Vector of variables
      integer, intent(in) :: spec_number !! Number of specification
      real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
      type(injelope), intent(out) :: envels !! Calculated envelopes

      type(critical_point), allocatable :: cps(:)

      real(pr) :: X(size(X0))
      integer :: ns
      real(pr) :: S
      real(pr) :: XS(max_points, size(X0))
      real(pr) :: del_S

      real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

      real(pr) :: z(size(X0) - 2)

      integer :: point, iters, n
      integer :: i
      integer :: funit_output
      character(len=254) :: fname_env

      allocate (cps(0))
      X = X0
      n = size(X0) - 2
      ns = spec_number
      S = X(ns)
      del_S = del_S0

      call get_z(X(n+2), z)

      ! ========================================================================
      !  Output file
      ! ------------------------------------------------------------------------
      env_number = env_number + 1

      write (fname_env, *) env_number
      fname_env = "env-2ph-PX"//"_"//trim(adjustl(fname_env))
      fname_env = trim(adjustl(ouput_path))//trim(fname_env)//".dat"

      open (newunit=funit_output, file=fname_env)
      write (funit_output, *) "#", T
      write (funit_output, *) "STAT", " iters", " ns", " alpha", " P", &
         (" lnK"//str(i), i=1,n), (" z"//str(i), i=1,n)
      write (funit_output, *) "X0", iters, ns, X(n + 2), exp(X(n + 1)), X(:n), z
      ! ========================================================================

      enveloop: do point = 1, max_points
         ! call progress_bar(point, max_points, advance=.false.)
         call full_newton(&
            f_injection, iters, X, ns, S, max_iters, F, dF, solvetol=solve_tol &
         )

         if (iters >= max_iters) then
            ! call progress_bar(point, max_points, advance=.true.)
            print *, "Breaking: Above max iterations"
            exit enveloop
         end if

         XS(point, :) = X

         call update_spec(X, ns, del_S, dF, dXdS)
         call fix_step_two_phases(X, ns, S, iters, del_S, dXdS)

         detect_critical: block
            real(pr) :: K(size(X0) - 2), Knew(size(X0) - 2), &
                        Xnew(size(X0)), fact
            real(pr) :: pc, alpha_c, dS_c
            integer :: max_changing
            fact = critical_fact ! 2.5

            Xnew = X + fact*dXdS*del_S

            K = X(:n)
            Knew = Xnew(:n)

            if (all(K*Knew < 0)) then
               max_changing = maxloc(abs(K - Knew), dim=1)

               dS_c = ( &
                      -K(max_changing)*(Xnew(max_changing) - X(max_changing)) &
                      /(Knew(max_changing) - K(max_changing)) &
                      )

               Xnew = X + dXdS * dS_c
               alpha_c = Xnew(n + 2)
               pc = exp(Xnew(n + 1))
               cps = [cps, critical_point(t, pc, alpha_c)]
               
               del_S = dS_c + critical_multiplier * dS_c

               write (funit_output, *) ""
               write (funit_output, *) ""
            end if
         end block detect_critical

         X = X + dXdS*del_S
         S = X(ns)

         if (any(break_conditions(X, ns, S, del_S))) then
            print *, "Breaking: ", break_conditions(X, ns, S, del_S)
            exit enveloop
         end if

         call get_z(X(n+2), z)
         write (funit_output, *) "SOL", iters, ns, X(n + 2), exp(X(n + 1)), &
            X(:n), z
      end do enveloop

      write (funit_output, *) ""
      write (funit_output, *) ""
      write (funit_output, *) "#critical"
      if (size(cps) > 0) then
         do i = 1, size(cps)
            write (funit_output, *) cps(i)%alpha, cps(i)%p
         end do
      else
         write (funit_output, *) "NaN NaN"
      end if

      close (funit_output)

      point = point - 1
      envels%z = z_0
      envels%z_inj = z_injection
      envels%logk = XS(:point, :n)
      envels%alpha = XS(:point, n + 2)
      envels%p = exp(XS(:point, n + 1))
      envels%critical_points = cps
   end subroutine

   subroutine update_spec(X, ns, del_S, dF, dXdS)
      real(pr), intent(in) :: X(:)
      integer, intent(in out) :: ns
      real(pr), intent(in out) :: del_S
      real(pr), intent(in) :: dF(size(X), size(X))
      real(pr), intent(in out) :: dXdS(size(X))

      real(pr) :: dFdS(size(X))
      integer  :: ns_new
      integer :: n

      dFdS = 0
      dFdS(size(dFdS)) = 1

      dXdS = solve_system(dF, dFdS)

      ns_new = maxloc(abs(dXdS), dim=1)

      if (ns_new /= ns) then
         ! translation of delS and dXdS to the  new specification variable
         del_S = dXdS(ns_new)*del_S
         dXdS = dXdS/dXdS(ns_new)
         ns = ns_new
      end if
   end subroutine

   subroutine fix_step_two_phases(X, ns, S, solve_its, del_S, dXdS)
      real(pr), intent(in) :: X(:)
      integer, intent(in) :: ns
      real(pr), intent(in) :: S
      integer, intent(in) :: solve_its
      real(pr), intent(in out) :: del_S
      real(pr), intent(in out) :: dXdS(size(X))


      real(pr) :: Xnew(size(X-1))
      real(pr) :: dP, dalpha
      real(pr) :: dP_tol, dalpha_tol
      integer :: n

      n = size(X) - 2

      del_S = sign(del_S_multiplier, del_S)*minval([ &
                                         max(sqrt(abs(X(ns)))/10, 0.1), &
                                         abs(del_S)*3/solve_its &
                                         ] &
                                 )

      if (injection_case == "dilution") del_S = 50*del_S

      Xnew = X + dXdS*del_S
      dP = exp(Xnew(n + 1)) - exp(X(n + 1))
      dalpha = Xnew(n + 2) - X(n + 2)

      if (X(n+2) > 1.9) then 
         dP_tol = max_dP/2
         dalpha_tol = max_dalpha/2
      else
         dP_tol = max_dp
         dalpha_tol = max_dalpha
      end if

      do while (abs(dP) > dp_tol .or. abs(dalpha) > dalpha_tol)
         dXdS = dXdS/2.0_pr

         Xnew = X + dXdS*del_S
         dP = exp(Xnew(n + 1)) - exp(X(n + 1))
         dalpha = Xnew(n + 2) - X(n + 2)
      end do
   end subroutine
   
   function break_conditions(X, ns, S, del_S)
      !! Set of conditions to break the tracing of a two phase line.
      real(pr) :: X(:) !! Vector of variables
      integer :: ns !! Number of specification
      real(pr) :: S !! Specification value
      real(pr) :: del_S !! \(\Delta S\)

      integer :: n
      real(pr) :: p, alpha
      logical, allocatable :: break_conditions(:)

      n = size(X) - 2
      p = exp(X(n + 1))
      alpha = X(n + 2)

      break_conditions = [ .false. &
                         ! p < 1e-15 .or. p > 5000, &
                         ! abs(del_S) < 1e-3 &
                         ]
   end function
   ! ===========================================================================

   ! ===========================================================================
   ! Three-phases
   ! ---------------------------------------------------------------------------
   subroutine F_injection_three_phases(Xvars, ns, S, F, dF)
      !! Function to solve at each point of a three phase envelope.
      !!
      !! The vector of variables X corresponds to:
      !! \( X = [lnKx_i, lnKy_i lnP, \alpha, \beta] \)
      !!
      !! While the equations are:
      !!
      !! \( F = [
      !!        lnKx_i - ln \phi_i(x, P, T) + ln \phi_i(w, P, T),
      !!        lnKy_i - ln \phi_i(y, P, T) + ln \phi_i(w, P, T),
      !!        \sum_{i=1}^N (w_i) - 1,
      !!        \sum_{i=1}^N (x_i - y_i),
      !!        X_{ns} - S
      !! ] \)
      use legacy_thermo_properties, only: TERMO
      use iso_fortran_env, only: error_unit
      real(pr), intent(in)  :: Xvars(:) !! Vector of variables
      integer, intent(in)  :: ns   !! Number of specification
      real(pr), intent(in)  :: S    !! Specification value
      real(pr), intent(out) :: F(size(Xvars)) !! Vector of functions valuated
      real(pr), intent(out) :: df(size(Xvars), size(Xvars)) !! Jacobian matrix

      ! Xvars variables
      real(pr) :: z((Size(Xvars)-3)/2)
      real(pr) :: Kx((Size(Xvars)-3)/2)
      real(pr) :: Ky((Size(Xvars)-3)/2)
      real(pr) :: P
      real(pr) :: beta
      real(pr) :: alpha

      ! Main phase 1 variables
      real(pr) :: Vx
      real(pr), dimension((Size(Xvars)-3)/2) :: x, lnfug_x, dlnphi_dt_x, dlnphi_dp_x
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_x

      ! Main phase 2 variables
      real(pr) :: Vy
      real(pr), dimension((Size(Xvars)-3)/2) :: y, lnfug_y, dlnphi_dt_y, dlnphi_dp_y
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_y

      ! Incipient phase variables
      real(pr) :: Vw
      real(pr), dimension((Size(Xvars)-3)/2) :: w, lnfug_w, dlnphi_dt_w, dlnphi_dp_w
      real(pr), dimension((Size(Xvars)-3)/2, (Size(Xvars)-3)/2) :: dlnphi_dn_w

      ! Derivative of z wrt alpha
      real(pr) :: dzda((Size(Xvars)-3)/2), dwda((Size(Xvars)-3)/2)

      ! Derivative of w wrt beta
      real(pr) :: dwdb((Size(Xvars)-3)/2)

      real(pr) :: dwdKx((Size(Xvars)-3)/2), dxdKx((Size(Xvars)-3)/2), dydKx((Size(Xvars)-3)/2)
      real(pr) :: dwdKy((Size(Xvars)-3)/2), dxdKy((Size(Xvars)-3)/2), dydKy((Size(Xvars)-3)/2)

      integer :: i, j, n

      n = (Size(Xvars)-3)/2

      Kx = exp(Xvars(1:n))
      Ky = exp(Xvars(n + 1:2*n))
      P = exp(Xvars(2*n + 1))
      alpha = Xvars(2*n + 2)
      beta = Xvars(2*n + 3)

      call get_z(alpha, z, dzda)
      ! if (any(z < 0)) z = 0

      w = z/(beta*Ky + (1 - beta)*Kx)
      x = w*Kx
      y = w*Ky

      call TERMO( &
         n, 0, 4, T, P, x, Vx, lnfug_x, dlnphi_dp_x, dlnphi_dt_x, dlnphi_dn_x &
         )
      call TERMO( &
         n, 0, 4, T, P, y, Vy, lnfug_y, dlnphi_dp_y, dlnphi_dt_y, dlnphi_dn_y &
         )
      call TERMO( &
         n, 0, 4, T, P, w, Vw, lnfug_w, dlnphi_dp_w, dlnphi_dt_w, dlnphi_dn_w &
         )

      F(1:n) = Xvars(1:n) + lnfug_x - lnfug_w
      F(n + 1:2*n) = Xvars(n + 1:2*n) + lnfug_y - lnfug_w

      F(2*n + 1) = sum(w) - 1
      F(2*n + 2) = sum(x - y)
      F(2*n + 3) = Xvars(ns) - S

      df = 0
      dwda = 1.0_pr/(beta*Ky + (1 - beta)*Kx)*dzda
      dwdb = z*(Kx - Ky)/((1 - beta)*Kx + beta*Ky)**2

      dwdKx = -z*(1 - beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKx = Kx*dwdKx + w
      dydKx = Ky*dwdKx

      dwdKy = -z*(beta)/(Ky*beta + (1 - beta)*Kx)**2
      dxdKy = Kx*dwdKy
      dydKy = Ky*dwdKy + w

      do i = 1, n
         do j = 1, n
            df(i, j) = Kx(j)*(dlnphi_dn_x(i, j)*dxdKx(j) &
                              - dlnphi_dn_w(i, j)*dwdKx(j))
            df(i + n, j) = Kx(j)*(dlnphi_dn_y(i, j)*dydKx(j) &
                                  - dlnphi_dn_w(i, j)*dwdKx(j))

            df(i, j + n) = Ky(j)*(dlnphi_dn_x(i, j)*dxdKy(j) &
                                  - dlnphi_dn_w(i, j)*dwdKy(j))
            df(i + n, j + n) = Ky(j)*(dlnphi_dn_y(i, j)*dydKy(j) &
                                      - dlnphi_dn_w(i, j)*dwdKy(j))
         end do

         df(i, i) = df(i, i) + 1
         df(i + n, i + n) = df(i + n, i + n) + 1
         df(i, 2*n + 2) = sum( &
                          Kx*dlnphi_dn_x(i, :)*dwda - dlnphi_dn_w(i, :)*dwda &
                          )
         df(i + n, 2*n + 2) = sum(Ky*dlnphi_dn_y(i, :)*dwda &
                                  - dlnphi_dn_w(i, :)*dwda)
         ! Derivatives wrt beta
         df(i, 2*n + 3) = sum(Kx*dlnphi_dn_x(i, :)*dwdb &
                              - dlnphi_dn_w(i, :)*dwdb)
         df(i + n, 2*n + 3) = sum(Ky*dlnphi_dn_y(i, :)*dwdb &
                                  - dlnphi_dn_w(i, :)*dwdb)
         
         df(2*n + 1, i) = Kx(i)*dwdKx(i)
         df(2*n + 1, i + n) = Ky(i)*dwdKy(i)

         df(2*n + 2, i) = Kx(i)*dxdKx(i) - Kx(i)*dydKx(i)
         df(2*n + 2, i + n) = Ky(i)*dxdKy(i) - Ky(i)*dydKy(i)
      end do

      ! Derivatives wrt P
      df(:n, 2*n + 1) = P*(dlnphi_dp_x - dlnphi_dp_w)
      df(n + 1:2*n, 2*n + 1) = P*(dlnphi_dp_y - dlnphi_dp_w)

      ! Derivatives wrt alpha
      df(2*n + 1, 2*n + 2) = sum(dwda)
      df(2*n + 2, 2*n + 2) = sum(Kx*dwda - Ky*dwda)

      ! Derivatives wrt beta
      df(2*n + 1, 2*n + 3) = sum(dwdb)
      df(2*n + 2, 2*n + 3) = sum(Kx*dwdb - Ky*dwdb)

      ! Derivatives wrt Xs
      df(2*n + 3, :) = 0
      df(2*n + 3, ns) = 1
   end subroutine

   subroutine injection_envelope_three_phase(X0, spec_number, del_S0, envels)
      use constants, only: ouput_path
      use ftools__io
      !! Subroutine to calculate Px phase envelopes via continuation method.
      !! Three phases version.
      real(pr), intent(in) :: X0(:) !! Vector of variables
      integer, intent(in) :: spec_number !! Number of specification
      real(pr), intent(in) :: del_S0 !! \(\Delta S_0\)
      type(injelope), intent(out) :: envels !! Calculated envelopes

      type(critical_point), allocatable :: cps(:)

      real(pr), target :: X(size(X0))
      integer :: ns
      real(pr) :: S
      real(pr) :: XS(max_points, size(X0))
      real(pr) :: del_S

      real(pr) :: F(size(X0)), dF(size(X0), size(X0)), dXdS(size(X0))

      real(pr), pointer :: lnP
      real(pr), pointer :: alpha
      real(pr), pointer :: beta
      real(pr), pointer :: lnKx(:)
      real(pr), pointer :: lnKy(:)
      real(pr) :: z(size(z_0))

      integer :: point, iters, n
      integer :: i
      integer :: funit_output
      character(len=254) :: fname_env

      allocate (cps(0))
      X = X0
      n = (size(X0) - 3)/2
      ns = spec_number
      S = X(ns)
      del_S = del_S0

      lnKx => X(:n)
      lnKy => X(n+1:2*n)
      lnP => X(2*n+1)
      alpha => X(2*n+2)
      beta => X(2*n+3)
      call get_z(alpha, z)

      ! ======================================================================
      !  Output file
      ! ----------------------------------------------------------------------
      env_number = env_number + 1

      write (fname_env, *) env_number
      fname_env = "env-3ph-PX"//"_"//trim(adjustl(fname_env))
      fname_env = trim(adjustl(ouput_path))//trim(fname_env)//".dat"

      open (newunit=funit_output, file=fname_env)
      write (funit_output, *) "#", T
      write (funit_output, *) "STAT", " iters", " ns", " alpha", " P", &
         " beta", (" lnKx"//str(i), i=1,n), (" lnKy"//str(i), i=1,n), (" z" //str(i), i=1,n)
      write (funit_output, *) "X0", iters, ns, X(2*n + 2), exp(X(2*n + 1)), &
         X(2*n + 3), X(:2*n)
      ! ======================================================================

      enveloop: do point = 1, max_points
         ! call progress_bar(point, max_points, advance=.false.)
         call full_newton(&
            F_injection_three_phases, iters, X, ns, S, max_iters, F, dF, solvetol=solve_tol &
         )
         if (iters >= max_iters) then
            ! call progress_bar(point, max_points, advance=.true.)
            print *, "Breaking: Above max iterations"
            exit enveloop
         end if

         call get_z(alpha, z)
         write (funit_output, *) "SOL", iters, ns, alpha, &
            exp(X(2*n + 1)), X(2*n + 3), X(:2*n), z
         XS(point, :) = X

         call update_spec(X, ns, del_S, dF, dXdS)

         fix_step: block
            real(pr) :: Xnew(size(X0))
            real(pr) :: dP, dalpha

            del_S = sign(del_S_multiplier_three_phase, del_S)*minval([ &
                                               max(abs(sqrt(X(ns))/10), 0.1_pr), &
                                               abs(del_S)*3/iters &
                                               ] &
                                               )

            if (injection_case == "dilution") del_S = 50*del_S

            Xnew = X + dXdS*del_S
            dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
            dalpha = Xnew(2*n + 2) - X(2*n + 2)

            do while (abs(dP) > max_dp .or. abs(dalpha) > max_dalpha)
               dXdS = dXdS/2.0_pr

               Xnew = X + dXdS*del_S
               dP = exp(Xnew(2*n + 1)) - exp(X(2*n + 1))
               dalpha = Xnew(2*n + 2) - X(2*n + 2)
            end do
         end block fix_step

         detect_critical: block
            real(pr) :: K((size(X0) - 3)/2), Knew((size(X0) - 3)/2), &
                        Xnew(size(X0)), fact
            real(pr) :: pc, alpha_c, dS_c, dXdS_in(size(X0))
            real(pr) :: bu_solve_tol
            integer :: max_changing, i
            logical :: found_critical

            fact = critical_fact
            found_critical = .false.

            loop: do i = 0, 1
               Xnew = X + fact*dXdS*del_S

               K = X(i*n + 1:(i + 1)*n)
               Knew = Xnew(i*n + 1:(i + 1)*n)

               max_changing = maxloc(abs(Knew - K), dim=1)

               if (all(K*Knew < 0)) then

                  print *, ""
                  print *, "CP ENV", env_number
                  print *, ""

                  dS_c = ( &
                         -k(max_changing) * (Xnew(ns) - X(ns)) &
                          /(Knew(max_changing) - K(max_changing)) &
                         )

                  Xnew = X + dXdS*dS_c
                  alpha_c = Xnew(2*n + 2)
                  pc = exp(Xnew(2*n + 1))
                  cps = [cps, critical_point(t, pc, alpha_c)]

                  del_S = dS_c + 1.7_pr * dS_c

                  write (funit_output, *) ""
                  write (funit_output, *) ""
                  exit loop
               end if
            end do loop
         end block detect_critical

         if (x(2*n + 3) > 1 .or. (x(2*n+3) < 0)) then
            ! call progress_bar(point, max_points, .true.)
            print *, "Breaking: positive 𝜷"
            exit enveloop
         end if

         X = X + dXdS*del_S
         S = X(ns)
         if (any(break_conditions_three_phases(X, ns, S, del_S)) .and. point > 10) then
            ! call progress_bar(point, max_points, .true.)
            print *, "Breaking: ", break_conditions_three_phases(X, ns, S, del_S)
            exit enveloop
         end if
      end do enveloop

      write (funit_output, *) ""
      write (funit_output, *) ""
      write (funit_output, *) "#critical"
      if (size(cps) > 0) then
         do i = 1, size(cps)
            write (funit_output, *) cps(i)%alpha, cps(i)%p
         end do
      else
         write (funit_output, *) "NaN NaN"
      end if

      close (funit_output)
      envels%z = z_0
      envels%z_inj = z_injection
      envels%logk = XS(:point, :n)
      envels%alpha = XS(:point, n + 2)
      envels%p = exp(XS(:point, n + 1))
      envels%critical_points = cps
   end subroutine
   
   function break_conditions_three_phases(X, ns, S, del_S)
      !! Set of conditions to break the tracing.
      real(pr) :: X(:) !! Variables vector
      integer :: ns !! Number of specification
      real(pr) :: S !! Value of specification
      real(pr) :: del_S

      integer :: n
      real(pr) :: p, alpha
      logical, allocatable :: break_conditions_three_phases(:)

      n = (size(X) - 3)/2
      p = exp(X(2*n + 1))
      alpha = X(2*n + 2)

      break_conditions_three_phases = [ &
                                       ! p < 0.001_pr .or. p > 5000, &
                                       abs(del_S) < 1e-5 &
                                      ]
   end function
   ! ===========================================================================
   
   subroutine get_z(alpha, z, dzda)
      !! Calculate the fluid composition based on an amount of addition
      !! of second fluid.
      !!
      !! The injection can be considered as two kinds of injection:
      !! - Displacement: \( z = \alpha z_i + (1-\alpha) z_0 \)
      !! - Addition:  \( z = \frac{\alpha z_i + (1-\alpha) z_0}{\sum_{i=1}^N \alpha z_i + (1-\alpha) z_0} \)
      real(pr), intent(in)  :: alpha !! Addition percentaje \( \alpha \)
      real(pr), intent(out) :: z(size(z_0)) !! New composition
      real(pr), optional, intent(out) :: dzda(size(z_0)) !! Derivative wrt \(\alpha\)

      select case (injection_case)
      case ("displace")
         z = (z_injection*alpha + (1.0_pr - alpha)*z_0)
         if (present(dzda)) dzda = z_injection - z_0
      case ("dilute")
         z = (z_injection*alpha + z_0)/sum(z_injection*alpha + z_0)
         if (present(dzda)) dzda = -(alpha*z_injection + z_0) &
                *sum(z_injection)/sum(alpha*z_injection + z_0)**2 &
                + z_injection/sum(alpha*z_injection + z_0)
      case default
         z = (z_injection*alpha + (1.0_pr - alpha)*z_0)
         if (present(dzda)) dzda = z_injection - z_0
      end select
   end subroutine

   function get_case(dew_envel, bub_envel) result(n_case)
      type(injelope), intent(in) :: dew_envel
      type(injelope), intent(in) :: bub_envel
      integer :: n_case
   end function

   function remove_duplicates(envels) result(clean_envels)
      !! From a set of envelopes check if they are the same line
      class(injelope) :: envels(:)
      type(injelope), allocatable :: clean_envels(:)

      if (size(envels) /= 1) then
      else
         clean_envels = envels
      end if
   end function

   function same_line(env1, env2)
      !! 
      class(injelope), intent(in) :: env1, env2
      logical :: same_line
   end function

   ! ===========================================================================
   ! Initialization procedures
   ! ---------------------------------------------------------------------------
   function px_two_phase_from_pt(&
      t_inj, pt_env_2, t_tol, del_S0) result(px_envels)
      !! Calculate two phase Px envelopes at a given injection temperature.
      !!
      !! Given an injection temperature `t_inj` and a base PT envelope 
      !! `pt_env_2`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.
      
      use linalg, only: interpol
      use stdlib_optval, only: optval
      
      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(envelope), intent(in) :: pt_env_2 !! Base PT envelope
      real(pr), intent(in) :: t_tol !! Absolute temperature tolerance
      real(pr), optional, intent(in) :: del_S0 !! First point \(\Delta S\)
      type(injelope), allocatable :: px_envels(:) !! Output Px envelope
      
      real(pr), allocatable :: ts_envel(:) !! Temperatures under tolerance 
      real(pr), allocatable :: k(:) !! K values
      real(pr), allocatable :: X(:) !! Vector of variables
      real(pr) :: alpha !! Amount of injection
      real(pr) :: p !! Pressure of ocurrence
      real(pr) :: pold !! Old pressure, used to assure no repeats

      integer :: i, idx, ns
      real(pr) :: del_S

      type(injelope) :: px_envel

      allocate(px_envels(0))
      del_S = optval(del_S0, 0.1_pr)
      pold = 1e9

      ts_envel = pack(pt_env_2%t, mask=abs(pt_env_2%t - t_inj) < t_tol)
      do i = 1, size(ts_envel)
         idx = findloc(pt_env_2%t, value=ts_envel(i), dim=1)
         p = interpol( &
               pt_env_2%t(idx), pt_env_2%t(idx + 1), &
               pt_env_2%p(idx), pt_env_2%p(idx + 1), &
               t_inj)

         if (abs(p - pold) < 2) cycle
         pold = p

         k = exp(interpol( &
                  pt_env_2%t(idx), pt_env_2%t(idx + 1), &
                  pt_env_2%logk(idx, :), pt_env_2%logk(idx + 1, :), &
                  t_inj))

         alpha = 0
         X = [log(K), log(P), alpha]
         ns = size(X)

         call injection_envelope(X, ns, del_S, px_envel)
         px_envels = [px_envels, px_envel]
      end do
   end function

   function px_three_phase_from_pt(t_inj, pt_env_3, t_tol, del_S0) result(envel)
      !! Calculate three phase Px envelopes at a given injection temperature.
      !!
      !! Given an injection temperature `t_inj` and a base PT envelope 
      !! `pt_env_3`, finds all the points on the PT envelope near `t_inj`, based
      !! on an absolute tolerance `t_tol`. These points are used as 
      !! initialization for calculation of Px envelopes.

      use linalg, only: interpol
      use stdlib_optval, only: optval
      use envelopes, only: PTEnvel3

      real(pr), intent(in) :: t_inj !! Injection temperature [K]
      type(PTEnvel3), intent(in) :: pt_env_3(:) !! Base PT envelopes
      real(pr), intent(in) :: t_tol !! Absolute temperature tolerance
      real(pr), optional, intent(in) :: del_S0 !! First point \(\Delta S\)
      type(injelope) :: envel !! Output Px envelope

      real(pr), allocatable :: ts_envel(:) !! Temperatures under tolerance 
      real(pr), allocatable :: kx(:), ky(:) !! K values
      real(pr), allocatable :: X(:) !! Vector of variables
      real(pr) :: alpha !! Amount of injection
      real(pr) :: beta  !! Main phases molar fraction
      real(pr) :: p !! Pressure of ocurrence
      real(pr) :: pold !! Old pressure, used to assure no repeats

      integer :: i, idx, ns, i_envel
      real(pr) :: del_S

      del_S = optval(del_S0, 0.05_pr)
      pold = 0

      do i_envel = 1, size(pt_env_3)
         print *, i_envel
         associate(pt => pt_env_3(i_envel))
         ts_envel = pack(pt%t, mask=abs(pt%t - t_inj) < t_tol)
         do i = 1, size(ts_envel)
            idx = findloc(pt%t, value=ts_envel(i), dim=1)
            p = interpol( &
                  pt%t(idx), pt%t(idx + 1), &
                  pt%p(idx), pt%p(idx + 1), &
                  t_inj)

            if (abs(p - pold) < 5) cycle
            pold = p

            kx = exp(interpol( &
                     pt%t(idx), pt%t(idx + 1), &
                     pt%lnkx(idx, :), pt%lnkx(idx + 1, :), &
                     t_inj))
            ky = exp(interpol( &
                     pt%t(idx), pt%t(idx + 1), &
                     pt%lnky(idx, :), pt%lnky(idx + 1, :), &
                     t_inj))
            beta = interpol( &
                     pt%t(idx), pt%t(idx + 1), &
                     pt%beta(idx), pt%beta(idx + 1), &
                     t_inj)
            
            alpha = 0
            X = [log(Kx), log(Ky), log(P), alpha, beta]
            ns = size(X) - 1

            print *, "Running isolated PX"
            call injection_envelope_three_phase(X, ns, del_S, envel)
         end do
         end associate
      end do
   end function

   function px_three_phase_from_inter(&
         inter, px_1, px_2, del_S0, beta0 &
         ) result(envels)
      use legacy_ar_models, only: nc
      use stdlib_optval, only: optval
      use linalg, only: point, interpol
      type(point), intent(in) :: inter
      type(injelope), intent(in) :: px_1, px_2
      type(injelope) :: envels(2)
      real(pr), optional :: del_S0
      real(pr), optional :: beta0

      integer :: i, j

      real(pr) :: lnKx(nc), lnKy(nc), alpha, p, beta, X(2*nc+3)
      real(pr) :: phase_y(nc), phase_x(nc)
      real(pr) :: del_S
      real(pr) :: z(nc), dzda(nc)
      integer :: ns

      del_S = optval(del_S0, -0.01_pr)
      beta = optval(beta0, 1.0_pr)

      i = inter%i
      j = inter%j

      alpha = inter%x
      p = inter%y

      lnKx = interpol( &
               px_1%alpha(i), px_1%alpha(i + 1), &
               px_1%logk(i, :), px_1%logk(i + 1, :), &
               alpha &
               )

      lnKy = interpol( &
               px_2%alpha(j), px_2%alpha(j + 1), &
               px_2%logk(j, :), px_2%logk(j + 1, :), &
               alpha &
               )

      call get_z(alpha, z, dzda)

      ! Bubble line composition
      phase_y = exp(lnKy)*z
      ! Dew line composition
      phase_x = exp(lnKx)*z

      ns = 2*nc + 3

      ! ==================================================================
      !  Line with incipient phase gas
      ! ------------------------------------------------------------------
      print *, "Three Phase: Gas"
      lnKx = log(phase_x/phase_y)
      lnKy = log(z/phase_y)
      X = [lnKx, lnKy, log(p), alpha, beta]
      call injection_envelope_three_phase(X, ns, del_S, envels(1))
      ! ==================================================================

      ! ==================================================================
      !  Line with incipient phase liquid
      ! ------------------------------------------------------------------
      print *, "Three Phase: Liquid"
      lnKx = log(phase_y/phase_x)
      lnKy = log(z/phase_x)
      X = [lnKx, lnKy, log(p), alpha, beta]
      call injection_envelope_three_phase(X, ns, del_S, envels(2))
   end function
   
   function px_hpl_line(alpha_0, p)
      ! Find a HPLL PX line at a given pressure, starting from a given alpha
      use legacy_ar_models, only: nc
      use legacy_thermo_properties, only: termo
      use linalg, only: solve_system
      real(pr), intent(in) :: alpha_0 !! Staring \(\alpha\) to search
      real(pr), intent(in) :: p !! Pressure of HPLL
      type(injelope) :: px_hpl_line !! Resulting HPLL line

      real(pr) :: diff
      real(pr) :: lnfug_z(nc), lnfug_y(nc), &
                  dlnphi_dz(nc, nc), dlnphi_dy(nc, nc), &
                  dlnphi_dp_z(nc), dlnphi_dp_y(nc)
      real(pr) :: z(nc), y(nc), v, k(nc)
      real(pr) :: p_in, alpha_in

      real(pr), allocatable :: x(:)
      real(pr) :: del_S0
      integer :: ns, ncomp

      find_hpl: do ncomp = nc, 1, -1
         alpha_in = alpha_0
         p_in = p

         y = 0
         y(ncomp) = 1.0

         diff = foo([alpha_in, p_in, y])
         return 
         do while (abs(diff) > 0.001 .and. p_in < 5000)
            p_in = p_in + 10.0_pr
            diff = foo([alpha_in, p_in, y])
         end do

         if (p_in >= 5000) then
            do while (abs(diff) > 0.001 .and. alpha_in > 0)
               alpha_in = alpha_in - 0.01_pr
               diff = foo([alpha_in, p_in, y])
            end do
         end if

         if (.true.) then!(alpha_in > 0 .or. p_in < 5000) then
            call optim
            print *, alpha_in, p_in
            k = 1/exp(lnfug_y - lnfug_z)

            X = [log(K), log(P_in), alpha_in]

            del_S0 = 0.1_pr
            ns = size(X) - 1
            call injection_envelope(X, ns, del_S0, px_hpl_line)
            exit find_hpl
         end if
      end do find_hpl

      contains 
         function foo(x) result(f)
            real(pr) :: x(:)
            real(pr) :: f

            if (x(1) > 1) x(1) = 0.97_pr
            ! alpha_in = x(1)
            p_in = x(2)
            y = abs(x(3:))

            call get_z(alpha_in, z)
            call termo(nc, 4, 1, t, p_in, z, v, philog=lnfug_z, dlphip=dlnphi_dp_z, fugn=dlnphi_dz)
            call termo(nc, 4, 1, t, p_in, y, v, philog=lnfug_y, dlphip=dlnphi_dp_y, fugn=dlnphi_dy)
            f = abs((log(z(ncomp)) + lnfug_z(ncomp)) - (log(y(ncomp)) + lnfug_y(ncomp)))
         end function

         subroutine optim
            use optimization, only: nm_opt
            real(pr) :: x0(size(z) + 2)
            integer :: stat
            real(pr) :: step(size(x0))

            x0 = [alpha_in, p_in, y]
            step = 0.1
            step(1) = 0.1
            step(2) = -100
            call nm_opt(foo, x0, stat, step)
         end subroutine
   end function
   ! ===========================================================================
end module