! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module kde_distribution_mod

use types_mod,               only : r8, missing_r8

use utilities_mod,           only : E_ERR, E_MSG, error_handler

use normal_distribution_mod, only : inv_cdf

use sort_mod,                only : sort, index_sort

use distribution_params_mod, only : distribution_params_type, deallocate_distribution_params

implicit none
private

public :: kde_cdf, kde_cdf_params, inv_kde_cdf, inv_kde_cdf_params,     &
          test_kde, obs_dist_types

type available_obs_dist_types
   integer  :: uninformative, normal, binomial, gamma, &
               inv_gamma, lognormal
end type

type(available_obs_dist_types), parameter :: obs_dist_types = &
available_obs_dist_types(uninformative=999, normal=998, binomial=997, &
   gamma=996, inv_gamma=995, lognormal=994)
character(len=512)          :: errstring
character(len=*), parameter :: source = 'kde_distribution_mod.f90'

contains

!---------------------------------------------------------------------------

function likelihood(x, y, obs_param, obs_dist_type) result(l)
   real(r8)             :: l  ! likelihood value
   real(r8), intent(in) :: x  ! state value
   real(r8), intent(in) :: y  ! obs value
   real(r8), intent(in) :: obs_param  ! meaning depends on obs_dist_type
   integer,  intent(in) :: obs_dist_type ! obs distribution type

   ! Evaluates the likelihood pdf(y | x) for various kinds of observation
   ! distributions. The value returned is not equal to the observation
   ! pdf evaluated at y, because normalization constants that don't depend
   ! on x are omitted.

   real(r8) :: gamma_shape, gamma_scale
   real(r8) :: inv_gamma_shape, inv_gamma_scale

   select case (obs_dist_type)
      case (obs_dist_types%uninformative)
         ! Uninformative observations have a likelihood equal to one.
         l = 1._r8
      case (obs_dist_types%normal)
         ! For a normal obs distribution, like_param is the obs error variance
         if (obs_param .le. 0._r8) then
            write(errstring, *) 'obs error sd .le. 0', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = exp( -0.5_r8 * (x - y)**2 / obs_param )
         end if
      case (obs_dist_types%binomial)
         ! For a binomial obs distribution 0<=x<=1 is the probability of
         ! observing y successes of obs_param total trials
         if (y .lt. 0._r8) then
            write(errstring, *) 'y value is negative with a binomial obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y .gt. obs_param) then
            write(errstring, *) 'successes greater than total trials with a binomial obs model ', y, obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif ((x .lt. 0._r8) .or. (x .gt. 1._r8)) then
            write(errstring, *) 'x outside [0,1] with a binomial obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = (x**y) * (1._r8 - x)**(obs_param - y)
         end if
      case (obs_dist_types%gamma)
         ! For a gamma obs distribution, the mean is x and the variance is obs_param * x^2, i.e.
         ! the obs error sd is sqrt(obs_param) times  the true value. If the error sd is p% of x,
         ! set obs_param = (p/100._r8)**2.
         ! For a gamma obs distribution, the likelihood is inverse gamma.
         if (x .lt. 0._r8) then
            write(errstring, *) 'x value is negative with a gamma obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y .le. 0._r8) then
            write(errstring, *) 'y value is non-positive with a gamma obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param .le. 0._r8) then
            write(errstring, *) 'obs variance is non-positive with a gamma obs model ', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            gamma_shape = 1._r8 / obs_param ! k
            gamma_scale = x * obs_param     ! theta
            if (x .eq. 0._r8) then
               l = 0._r8 ! Technically x must be > 0, but just setting l=0 is convenient.
            else
               l = (y / gamma_scale)**gamma_shape * exp(-y / gamma_scale)
            end if
         end if
      case (obs_dist_types%inv_gamma)
         ! For an inverse gamma obs distribution, the mean is x and the variance is obs_param * x^2,
         ! i.e. the obs error sd is sqrt(obs_param) times  the true value. If the error sd is p% of x,
         ! set obs_param = (p/100._r8)**2.
         ! For an inverse gamma obs distribution, the likelihood is gamma.
         if (x .lt. 0._r8) then
            write(errstring, *) 'x value is negative with an inverse gamma obs model ', x
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (y .le. 0._r8) then
            write(errstring, *) 'y value is non-positive with an inverse gamma obs model ', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param .le. 0._r8) then
            write(errstring, *) 'obs variance is non-positive with an inverse gamma obs model ', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            inv_gamma_shape = (1._r8 / obs_param) + 2._r8   ! alpha
            inv_gamma_scale = x * (inv_gamma_shape - 1._r8) ! beta
            if (x .eq. 0._r8) then
               l = 0._r8 ! Technically x must be > 0, but just setting l=0 is convenient.
            else
               l = (inv_gamma_scale**inv_gamma_shape) * exp( - inv_gamma_scale / y )
            end if
         end if
      case (obs_dist_types%lognormal)
         ! For a lognormal obs distribution, ln(y) is normal with mean x and variance obs_param.
         ! The likelihood is normal.
         if (y .le. 0._r8) then
            write(errstring, *) 'y value is non-positive with a lognormal obs model', y
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         elseif (obs_param .le. 0._r8) then
            write(errstring, *) 'obs error sd .le. 0', obs_param
            call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
         else
            l = exp( -0.5_r8 * (x - log(y))**2 / obs_param )
         end if
      case DEFAULT
         write(errstring, *) 'likelihood called with unrecognized obs_dist_type ', obs_dist_type
         call error_handler(E_MSG, 'kde_distribution_mod:likelihood', errstring, source)
   end select

end function likelihood

!---------------------------------------------------------------------------

elemental function biweight_kernel(x) ! The fact that this is elemental is not (yet) used
   real(r8) :: biweight_kernel
   real(r8), intent(in) :: x
   real(r8), parameter  :: norm_const = 15._r8 / 16._r8
   biweight_kernel = norm_const * max(0._r8, 1._r8 - x**2)**2

end function biweight_kernel

!---------------------------------------------------------------------------

elemental function biweight_cdf(x) ! The fact that this is elemental is not (yet) used
   real(r8) :: biweight_cdf
   real(r8), intent(in) :: x
   real(r8), parameter  :: norm_const = 1._r8 / 16._r8
   biweight_cdf = min(1._r8, max(0._r8,                                    &
         norm_const * ((1._r8 + x)**3) * (8._r8 + 3._r8 * x * (x - 3._r8)) &
                                                                          ))

end function biweight_cdf

!---------------------------------------------------------------------------

elemental subroutine boundary_correction(x, lx, mx)
   real(r8), intent(in)  :: x
   real(r8), intent(out) :: lx, mx

   ! Boundary correction for kde, from Jones (1993) and Jones & Foster (1996)
   ! The fact that this is elemental is not (yet) used

   real(r8) :: denom

   denom = (5._r8 * x**4 - 40._r8 * x**3 + 126._r8 * x**2 - 168._r8 * x + 81._r8) * (x+1._r8)**5
   lx = 64._r8 * (15._r8 * x**4 - 45._r8 * x**3 + 48._r8 * x**2 - 24._r8 * x + 8._r8) / denom
   mx =-1120._r8 * ((x - 1._r8)**3) / denom

end

!---------------------------------------------------------------------------

function kde_pdf(x, p)
   real(r8)                                   :: kde_pdf
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p

   ! Returns the kernel estimate of the pdf evaluated at x.

   integer  :: i
   real(r8) :: u_lower, lx_lower, mx_lower
   real(r8) :: u_upper, lx_upper, mx_upper

   kde_pdf = 0._r8 ! Initialize
   do i=1, p%ens_size ! Reduction loop - parallelizable
      u_lower = 0._r8; lx_lower = 1._r8 ; mx_lower = 0._r8
      u_upper = 0._r8; lx_upper = 1._r8 ; mx_upper = 0._r8
      if (p%bounded_below) then ! Bounded below
         u_lower = min( 1._r8, max( 0._r8, (x - p%lower_bound) / p%more_params(i) ) ) ! p%more_params(i) holds kernel width for ensemble member i
         call boundary_correction(u_lower, lx_lower, mx_lower)
      end if
      if (p%bounded_above) then ! Bounded above
         u_upper = min( 1._r8, max( 0._r8, (p%upper_bound - x) / p%more_params(i) ) )
         call boundary_correction(u_upper, lx_upper, mx_upper)
      end if
      kde_pdf = kde_pdf + (1._r8 / p%more_params(i)) * &
                          (lx_lower + mx_lower * u_lower) * &
                          (lx_upper + mx_upper * u_upper) * &
                          biweight_kernel( (x - p%ens(i)) / p%more_params(i) )
   end do
   kde_pdf = kde_pdf / (p%ens_size * p%more_params(p%ens_size + 1)) ! p%more_params(end) normalizes the pdf

end function kde_pdf

!-----------------------------------------------------------------------

subroutine get_kde_bandwidths(ens_size, ens, bandwidths)
   integer,         intent(in)   :: ens_size
   real(r8),        intent(in)   :: ens(ens_size)
   real(r8),        intent(out)  :: bandwidths(ens_size)

   real(r8)                      :: ens_mean, ens_sd, h0, g
   real(r8), dimension(ens_size) :: f_tilde, d, lambda
   integer                       :: i, k

   ens_mean = sum(ens) / ens_size
   ens_sd   = sqrt( sum( (ens - ens_mean)**2 ) / (ens_size - 1._r8) )
   h0 = 2.36_r8 * ens_sd / (ens_size**0.2_r8) ! This would be the kernel width if the widths were not adaptive.
                                              ! It would be better to use min(sd, iqr/1.34) but don't want to compute iqr
   k = floor( sqrt( real(ens_size) ) ) ! distance to kth nearest neighbor used to set bandwidth
   do i=1,ens_size
      d(:) = sort( abs( ens(:) - ens(i) ) ) ! Sorted neighbor distances
      f_tilde(i) = 0.5_r8 * real(k, r8) / (real(ens_size, r8) * d(k+1)) ! Initial density estimate
   end do
   g = product( f_tilde )**(1._r8 / real(ens_size, r8) )
   lambda(:) = sqrt( g / f_tilde(:) )
   bandwidths(:) = h0 * lambda(:)

end subroutine get_kde_bandwidths

!---------------------------------------------------------------------------

function integrate_pdf(x, p) result(q)
   real(r8)                                   :: q
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p

   ! Uses quadrature to approximate \int_{-\infty}^x l(x; y) p(s) ds where
   ! p(s) is the prior pdf and l(x; y) is the likelihood. The interval is
   ! broken up into sub-intervals whose boundaries are either one of the bounds,
   ! or the edge of support of one of the kernels, or the value of x. On each
   ! sub-interval the integral is approximated using Gauss-Legendre quadrature
   ! with 5 points. When the likelihood is flat and the boundaries are far
   ! from the ensemble, the result is exact up to roundoff error.

   real(r8) :: y
   real(r8) :: obs_param ! See likelihood function for interpretation
   integer  :: obs_dist_type  ! See likelihood function for interpretation
   real(r8) :: edges(2*p%ens_size)
   real(r8) :: left, right, xi ! edges of current sub-interval, quadrature point
   ! real(r8), parameter :: chi(3) = [-sqrt(3._r8/5._r8), 0._r8, sqrt(3._r8/5._r8)] ! Gauss quadrature points
   ! real(r8), parameter :: w(3)   = [     5._r8/9._r8, 8._r8/9._r8, 5._r8/9._r8 ] ! GQ weights
   real(r8), parameter :: chi(5) = [-sqrt(5._r8 + 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                                    -sqrt(5._r8 - 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                                    0._r8, &
                                     sqrt(5._r8 - 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8, &
                                     sqrt(5._r8 + 2._r8 * sqrt(10._r8 / 7._r8)) / 3._r8] ! Gauss quadrature points
   real(r8), parameter :: w(5)   = [(322._r8 - 13._r8 * sqrt(70._r8)) / 900._r8, &
                                    (322._r8 + 13._r8 * sqrt(70._r8)) / 900._r8, &
                                    128._r8/225._r8, &
                                    (322._r8 + 13._r8 * sqrt(70._r8)) / 900._r8, &
                                    (322._r8 - 13._r8 * sqrt(70._r8)) / 900._r8] ! GQ weights
   integer  :: i, k

   ! Unpack obs info from param struct
   y         = p%more_params(p%ens_size + 2)
   obs_param = p%more_params(p%ens_size + 3)
   obs_dist_type  = p%more_params(p%ens_size + 4)

   edges(1:p%ens_size)                = p%ens(:) - p%more_params(1:p%ens_size)
   edges(p%ens_size + 1:2*p%ens_size) = p%ens(:) + p%more_params(1:p%ens_size)
   edges(:) = sort(edges(:)) ! If bandwidths were constant we would not need to sort

   ! If x is outside the support of the pdf then we can skip the quadrature.
   ! Note that when bounded_below=false I still set lower_bound=edges(1).
   left = max(edges(1), p%lower_bound)
   if (x .le. left) then
      q = 0._r8
      return
   end if

   ! Important to use x .gt. upper_bound here because I use
   ! x = upper_bound to compute the normalization constant.
   ! When bounded_above=false I still set upper_bound=maxval(edges).
   if (x .gt. p%upper_bound) then
      q = 1._r8
      return
   end if

   ! If we haven't returned yet, then there is at least one subinterval.
   q = 0._r8
   i = 1
   right = min(x, edges(2)) ! left was computed above
   do k=1,5
      xi = 0.5_r8 * ((right - left) * chi(k) + left + right)
      q  = q + 0.5_r8 * (right - left) * w(k) * kde_pdf(xi, p) * &
           likelihood(xi, y, obs_param, obs_dist_type)
   end do
   do while ((x .gt. right) .and. (i+1 .lt. 2*p%ens_size))
      i     = i + 1
      left  = right
      right = min(x, edges(i+1))
      do k=1,5
         xi = 0.5_r8 * ((right - left) * chi(k) + left + right)
         q  = q + 0.5_r8 * (right - left) * w(k) * kde_pdf(xi, p) * &
              likelihood(xi, y, obs_param, obs_dist_type)
      end do
   end do
   ! Note that it is possible to have maxval(edges) < x < upper_bound,
   ! but that last sub-interval from maxval(edges) to x has zero integral,
   ! so it can be safely skipped.

end function integrate_pdf

!-----------------------------------------------------------------------

subroutine pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
   ens, y, obs_param, obs_dist_type, p)
   integer,                        intent(in)  :: ens_size
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: ens(ens_size)
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   type(distribution_params_type), intent(out) :: p

   real(r8) :: bandwidths(ens_size)
   real(r8) :: edge
   logical  :: needs_normalization

   ! Set the fixed storage parameters in the distribution_params_type
   p%ens_size = ens_size

   ! Allocate space needed for the parameters
   allocate(p%ens(ens_size))
   allocate(p%more_params(ens_size+4))

   ! p%more_params(1:ens_size) are the kernel bandwidths
   ! p%more_params(ens_size + 1) is the normalization constant for the pdf
   ! p%more_params(ens_size + 2) is the observation value y; not used for prior
   ! p%more_params(ens_size + 3) is the observation parameter (see likelihood for details); not used for prior
   ! p%more_params(ens_size + 4) is the observation error distribution type. For prior set to uninformative

   ! Save the ensemble values, sorted now so that they don't have to be re-sorted for the first guess
   p%ens(:) = sort(ens(:))

   ! Store the kernel bandwidths in the more_params array
   ! Important to make sure that p%ens and p%more_params are sorted in same order by passing p%ens rather than ens
   call get_kde_bandwidths(ens_size, p%ens, bandwidths)
   p%more_params(1:ens_size) = bandwidths(:)

   ! If the ensemble is sufficiently far from the boundary, then
   ! we can use the unbounded code, which is cheaper.
   edge = minval(p%ens(:) - 2*bandwidths(:))
   if (bounded_below .and. (edge .gt. lower_bound)) then
      p%bounded_below = .false.
      p%lower_bound = minval(p%ens(:) - bandwidths(:))
   else
      p%bounded_below = bounded_below
      p%lower_bound = lower_bound
   end if

   edge = maxval(p%ens(:) + 2*bandwidths(:))
   if (bounded_above .and. (edge .lt. upper_bound)) then
      p%bounded_above = .false.
      p%upper_bound = maxval(p%ens(:) + bandwidths(:))
   else
      p%bounded_above = bounded_above
      p%upper_bound = upper_bound
   end if

   ! Pack obs information
   p%more_params(ens_size + 2) = y
   p%more_params(ens_size + 3) = obs_param
   p%more_params(ens_size + 4) = obs_dist_type

   ! Get the normalization constant
   p%more_params(ens_size + 1) = 1._r8 ! This value is used below
   needs_normalization = .false. .or. p%bounded_below
   needs_normalization = needs_normalization .or. p%bounded_above
   needs_normalization = needs_normalization .or. &
                         .not.(obs_dist_type .eq. obs_dist_types%uninformative)
   if (needs_normalization) then
      ! quadrature-based, so need to normalize. The call below uses
      ! quadrature to integrate the pdf from -infty to the upper end
      ! of the support using a normalization constant of 1, then
      ! packs the result into p%more_params(ens_size+1)
      p%more_params(ens_size + 1) = integrate_pdf(p%upper_bound, p)
   end if

end subroutine pack_kde_params

!-----------------------------------------------------------------------

function kde_cdf_params(x, p) result(quantile)
   real(r8),                       intent(in) :: x
   type(distribution_params_type), intent(in) :: p
   real(r8)                                   :: quantile

   ! Returns the cumulative distribution of the kernel density estimate
   ! at the value x.

   integer  :: i
   logical  :: use_analytical_cdf
   ! It is a waste of memory to unpack p%more_params, but it aids readability
   real(r8) :: bandwidths(p%ens_size)
   real(r8) :: y, obs_param
   integer  :: obs_dist_type

   bandwidths(:) = p%more_params(1:p%ens_size)
   y             = p%more_params(p%ens_size + 2)
   obs_param     = p%more_params(p%ens_size + 3)
   obs_dist_type = p%more_params(p%ens_size + 4)

   quantile = 0._r8
   ! If the likelihood is uninformative and the distribution is unbounded, we can evaluate
   ! the cdf analytically (instead of using quadrature) to save computation.
   use_analytical_cdf = ( (.not.p%bounded_below) .and. (.not.p%bounded_above) .and. &
                          (obs_dist_type .eq. obs_dist_types%uninformative) )
   if (use_analytical_cdf) then
      do i=1,p%ens_size
         quantile = quantile + biweight_cdf( (x - p%ens(i)) / bandwidths(i) ) / bandwidths(i)
      end do
   else ! Compute cdf using quadrature
      quantile  = integrate_pdf(x, p)
   end if

end function kde_cdf_params

!---------------------------------------------------------------------------

function kde_cdf(x, ens, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, y, obs_param, obs_dist_type) result(quantile)
   real(r8),                       intent(in)  :: x
   integer,                        intent(in)  :: ens_size
   real(r8),                       intent(in)  :: ens(ens_size)
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   real(r8)                                    :: quantile

   ! Returns the cumulative distribution of the kernel density estimate
   ! at the value x.

   type(distribution_params_type) :: p

   call pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
      ens, y, obs_param, obs_dist_type, p)
   quantile = kde_cdf_params(x,p)

end function kde_cdf

!---------------------------------------------------------------------------

function inv_kde_cdf(quantile, ens, ens_size, bounded_below, bounded_above, &
   lower_bound, upper_bound, y, obs_param, obs_dist_type) result(x)
   real(r8),                       intent(in)  :: quantile
   integer,                        intent(in)  :: ens_size
   real(r8),                       intent(in)  :: ens(ens_size)
   logical,                        intent(in)  :: bounded_below,        bounded_above
   real(r8),                       intent(in)  :: lower_bound,          upper_bound
   real(r8),                       intent(in)  :: y
   real(r8),                       intent(in)  :: obs_param
   integer,                        intent(in)  :: obs_dist_type
   real(r8)                                    :: x

   ! Returns the value x such that cdf(x) = quantile.

   type(distribution_params_type) :: p

   call pack_kde_params(ens_size, bounded_below, bounded_above, lower_bound, upper_bound, &
      ens, y, obs_param, obs_dist_type, p)

   x = inv_cdf(quantile, kde_cdf_params, inv_kde_first_guess_params, p)

end function inv_kde_cdf

!---------------------------------------------------------------------------

function inv_kde_cdf_params(quantile, p) result(x)
   real(r8)                                   :: x
   real(r8),                       intent(in) :: quantile
   type(distribution_params_type), intent(in) :: p

   ! Returns the value x such that cdf(x) = quantile.

   x = inv_cdf(quantile, kde_cdf_params, inv_kde_first_guess_params, p)

end function inv_kde_cdf_params

!---------------------------------------------------------------------------

function inv_kde_first_guess_params(quantile, p) result(x)
   real(r8)                                   :: x
   real(r8),                       intent(in) :: quantile
   type(distribution_params_type), intent(in) :: p

   ! This first-guess subroutine evaluates the cdf at the ensemble members,
   ! then finds a pair of ensemble members whose quantiles bracket the
   ! target quantile, then sets the first guess to a convex combination of
   ! these two ensemble members.

   real(r8) :: q0, q1
   integer  :: i

   q0 = kde_cdf_params(p%ens(1), p)
   if (q0 .ge. quantile) then
      x = p%ens(1)
      return
   end if
   do i=1,p%ens_size-1
      q1 = kde_cdf_params(p%ens(i+1), p)
      if (q0 .lt. quantile) then ; if (q1 .ge. quantile) then
         x = ((q1 - quantile) * p%ens(i) + (quantile - q0) * p%ens(i+1)) / &
             (q1 - q0)
         return
      end if ; end if
      q0 = q1
   end do
   x = p%ens(p%ens_size)

end function inv_kde_first_guess_params

!-----------------------------------------------------------------------

subroutine test_kde
   ! This routine provides limited tests of the numerics in this module. It tests
   ! the boundary correction function, the likelihood, the bandwidth selection,
   ! and the cdf inversion. It uses an ensemble [-1, 1]. It tests the bandwidth selection and the cdf inversion
   ! with zero, one, or two bounds at [-2, 2]. Failing these tests suggests a
   ! serious problem. Passing them does not indicate that there are acceptable
   ! results for all possible inputs.

   ! TODO: Add test for likelihood(s)

   integer,             parameter :: ens_size = 2
   real(r8),  dimension(ens_size) :: ensemble, bandwidths, target_bandwidths
   type(distribution_params_type) :: p
   real(r8)                       :: x, y, inv, max_diff, lx, mx, like, obs_param
   integer                        :: i, obs_dist_type

   ! Test the boundary correction code against a Matlab implementation (R2020a)
   x = 0.5_r8
   call boundary_correction(x, lx, mx)
   write(*, *) '----------------------------'
   write(*, *) 'test boundary correction'
   write(*, *) 'abs difference in lx is ', abs(1.172396660294006_r8 - lx)
   write(*, *) 'abs difference in mx is ', abs(0.7742242096281174_r8 - mx)
   write(*, *) 'abs differences should be less than 1e-15'

   ! Test uninformative likelihood
   max_diff = -1.0_r8
   y = 1._r8
   obs_param = 1._r8
   obs_dist_type = obs_dist_types%uninformative
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      like = likelihood(x, y, obs_param, obs_dist_type)
      max_diff = max(abs(1._r8 - like), max_diff)
   end do
   write(*, *) '----------------------------'
   write(*, *) 'Uninformative likelihood test'
   write(*, *) 'max difference in likelihood is ', max_diff
   write(*, *) 'max difference should be less than 1e-15'

   ! Test normal likelihood
   x = 0.5_r8
   y = 0._r8
   obs_param = 1._r8
   obs_dist_type = obs_dist_types%normal
   like = likelihood(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Normal likelihood test'
   write(*, *) 'abs difference in likelihood is ', abs(0.8824969025845955_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test binomial obs distribution (beta likelihood)
   x = 0.25_r8
   y = 3._r8
   obs_param = 5._r8
   obs_dist_type = obs_dist_types%binomial
   like = likelihood(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Binomial obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.0087890625_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test gamma obs distribution (inverse gamma likelihood)
   y = 1._r8
   x = 2._r8
   obs_param = 3._r8
   obs_dist_type = obs_dist_types%gamma
   like = likelihood(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Gamma obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.4658368455179406_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test inverse gamma obs distribution (gamma likelihood)
   y = 3._r8
   x = 2._r8
   obs_param = 4._r8
   obs_dist_type = obs_dist_types%inv_gamma
   like = likelihood(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Inverse gamma obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(3.415489474106968_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test lognormal obs distribution (normal likelihood)
   y = 1._r8
   x = 2._r8
   obs_param = 4._r8
   obs_dist_type = obs_dist_types%lognormal
   like = likelihood(x, y, obs_param, obs_dist_type)
   write(*, *) '----------------------------'
   write(*, *) 'Lognormal obs distribution test'
   write(*, *) 'abs difference in likelihood is ', abs(0.6065306597126334_r8 - like)
   write(*, *) 'abs difference should be less than 1e-15'

   ! Test bandwidth selection
   target_bandwidths(:) = 8.217997317515410_r8 * [1._r8, 1._r8]

   ensemble(:) = [-1._r8, 1._r8]

   call get_kde_bandwidths(ens_size, ensemble, bandwidths)

   ! Compare computed bandwidths to exact
   write(*, *) 'kde bandwidths test: Absolute value of differences should be less than 1e-15'
   do i = 1, ens_size
      write(*, *) i, abs(bandwidths(i) - target_bandwidths(i))
   end do

   ! Test the inversion of the cdf over the entire support of the pdf, unbounded
   call pack_kde_params(ens_size, .false., .false., 0._r8, 0._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      y = kde_cdf_params(x, p)
      inv = inv_kde_cdf_params(y, p)
      max_diff = max(abs(x-inv), max_diff)
   end do

   write(*, *) '----------------------------'
   write(*, *) 'Unbounded cdf/icdf test'
   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'max difference should be less than 1e-11'

   call deallocate_distribution_params(p)

   ! Test the inversion of the cdf over the entire support of the pdf, bounded below
   call pack_kde_params(ens_size, .true., .false., -2._r8, 0._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if (x .le. -2._r8) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test bounded below'
   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'max difference should be less than 1e-11'

   call deallocate_distribution_params(p)


   ! Test the inversion of the cdf over the entire support of the pdf, bounded above
   call pack_kde_params(ens_size, .false., .true., 0._r8, 2._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if (x .ge. 2._r8) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test bounded below'
   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'max difference should be less than 1e-11'

   call deallocate_distribution_params(p)

   ! Test the inversion of the cdf over the entire support of the pdf, doubly bounded
   call pack_kde_params(ens_size, .true., .true., -2._r8, 2._r8, &
      ensemble, 0._r8, 1._r8, obs_dist_types%uninformative, p)

   max_diff = -1.0_r8
   do i = 0, 1000
      x = ((real(i, r8) - 500.0_r8) / 500.0_r8) * (1._r8 + target_bandwidths(1))
      if ((x .le. -2._r8) .or. (x .ge. 2._r8)) then
         cycle
      else
         y = kde_cdf_params(x, p)
         inv = inv_kde_cdf_params(y, p)
         max_diff = max(abs(x-inv), max_diff)
      end if
   end do

   write(*, *) '----------------------------'
   write(*, *) 'cdf/icdf test bounded below'
   write(*, *) 'max difference in inversion is ', max_diff
   write(*, *) 'max difference should be less than 1e-11'

   ! Test the quadrature: Construct a case with bounds, but where the bounds are
   ! far enough from the data that they are not used. In this case the kernel
   ! density estimate should integrate exactly to one.
   p%lower_bound = -20._r8
   p%upper_bound =  20._r8

   y = integrate_pdf(20._r8, p)
   write(*, *) '----------------------------'
   write(*, *) 'test quadrature'
   write(*, *) 'abs difference is ', abs(1._r8 - y)
   write(*, *) 'abs difference should be less than 1e-15'

   call deallocate_distribution_params(p)

end subroutine test_kde

!-----------------------------------------------------------------------

end module kde_distribution_mod