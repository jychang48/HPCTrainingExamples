module norm_mod
  use prec_mod
  use mesh_mod, only: mesh_t
  implicit none

  private
  public :: norm

contains

  function norm(mesh, u) result(norm_val)
    type(mesh_t), intent(inout) :: mesh
    real(kind=RK), intent(inout) :: u(:,:)
    real(kind=RK) :: norm_val
    integer(kind=IK) :: i,j

    norm_val = 0._RK

    do j = 1,mesh%n_y
      do i = 1,mesh%n_x
        norm_val = norm_val + u(i,j)**2*mesh%dx*mesh%dy
      end do
    end do

    norm_val = sqrt(norm_val)/(mesh%n_x*mesh%n_y)
  end function norm

end module norm_mod
