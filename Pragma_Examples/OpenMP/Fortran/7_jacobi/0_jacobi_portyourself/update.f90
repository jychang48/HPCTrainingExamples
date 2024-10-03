module update_mod
  use prec_mod
  use mesh_mod, only: mesh_t
  implicit none

  private
  public :: update

contains

  subroutine update(mesh,rhs,au,u,res)
    type(mesh_t), intent(inout) :: mesh
    real(kind=RK), intent(inout) :: rhs(:,:)
    real(kind=RK), intent(inout) :: au(:,:)
    real(kind=RK), intent(inout) :: u(:,:)
    real(kind=RK), intent(inout) :: res(:,:)
    integer(kind=IK) :: i,j
    real(kind=RK) :: temp

    do j = 1,mesh%n_y
      do i = 1,mesh%n_x
        temp = rhs(i,j) - au(i,j)
        res(i,j) = temp
        u(i,j) = u(i,j) + temp*mesh%factor
      end do
    end do
  end subroutine

end module update_mod
