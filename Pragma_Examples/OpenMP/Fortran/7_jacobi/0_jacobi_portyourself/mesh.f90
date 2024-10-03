module mesh_mod
  use, intrinsic :: ISO_Fortran_env, only: stdout=>output_unit
  use prec_mod
  implicit none

  private
  public :: mesh_t, init_mesh

  real(kind=RK), parameter :: x_min=-0.5_RK, x_max=0.5_RK
  real(kind=RK), parameter :: y_min=-0.5_RK, y_max=0.5_RK
  real(kind=RK), parameter :: l_x=x_max-x_min, l_y=y_max-y_min

  type :: mesh_t
    integer(kind=IK) :: n_x, n_y
    real(RK) ::dx, dy, factor, invdx2, invdy2
    real(RK), allocatable :: x(:), y(:)
  end type mesh_t

contains

  subroutine init_mesh(this,n_x,n_y)
    type(mesh_t), intent(inout) :: this
    integer(kind=IK) :: n_x, n_y, i

    this%N_x = n_x
    this%N_y = n_y
    this%dx = l_x/(n_x+1._RK)
    this%dy = l_y/(n_y+1._RK)
    this%invdx2 = this%dx**-2
    this%invdy2 = this%dy**-2
    this%factor = (2._RK/this%dx**2+2._RK/this%dy**2)**-1

    allocate(this%x(0:n_x+1),this%y(0:n_y+1))
    this%x(0) = x_min
    do i = 1,n_x
      this%x(i) = x_min + i*this%dx
    end do
    this%x(n_x+1) = x_max
    this%y(0) = y_min
    do i = 1,n_y
      this%y(i) = y_min + i*this%dy
    end do
    this%y(n_y+1) = y_max

    write(stdout,'(A,I5,A,I5)') 'Domain size: ',n_x,' x ',n_y

  end subroutine init_mesh

end module mesh_mod
