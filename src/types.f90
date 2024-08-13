module dtypes
   use constants, only: pr
   use ftools__io

   implicit none

   private
   public :: AbsEnvel
   public :: envelope
   public :: critical_point

   type, abstract :: AbsEnvel
   end type
   
   type :: critical_point
      real(pr) :: t
      real(pr) :: p
      real(pr) :: alpha
   end type critical_point

   type :: envelope
      real(pr), allocatable :: vars(:, :)  !! Value of the set of variables at each point
      real(pr), allocatable :: z(:) !! Global composition
      real(pr), allocatable :: t(:) !! Temperature points
      real(pr), allocatable :: p(:) !! Pressure points
      real(pr), allocatable :: logk(:, :) !! ln(K) for each point
      real(pr), allocatable :: logphi(:, :) !! lnphi for each point
      type(critical_point), allocatable :: critical_points(:) !! Critical points
   end type envelope
   type :: envelope_nano
      real(pr), allocatable :: vars(:, :)  !! Value of the set of variables at each point
      real(pr), allocatable :: z(:) !! Global composition
      real(pr), allocatable :: t(:) !! Temperature points
      real(pr), allocatable :: pvap(:) !! Pressure of vapor phase (probably incipient) points
      real(pr), allocatable :: pliq(:) !! Pressure of liquid phase (probably constant) points
      real(pr), allocatable :: pcap(:) !! Capilar Pressure points
      real(pr), allocatable :: vy(:) !! Volume of vapor phase (probably incipient) points
      real(pr), allocatable :: vx(:) !! Volume of liquid phase (probably constant) points
      real(pr), allocatable :: logk(:, :) !! ln(K) for each point
      real(pr), allocatable :: logphi(:, :) !! lnphi for each point
      type(critical_point), allocatable :: critical_points(:) !! Critical points
   end type envelope_nano
end module