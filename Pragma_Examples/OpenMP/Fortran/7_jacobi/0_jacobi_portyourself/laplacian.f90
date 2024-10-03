module laplacian_mod
  use prec_mod
  use mesh_mod, only: mesh_t
  implicit none

  private
  public :: laplacian

contains

  subroutine laplacian(mesh,u,au)
    type(mesh_t), intent(inout) :: mesh
    real(kind=RK), intent(inout) :: u(:,:)
    real(kind=RK), intent(inout) :: au(:,:)
    integer(kind=IK) :: i,j

    do j = 2,mesh%n_y-1
      do i = 2,mesh%n_x-1
        au(i,j) = (-u(i-1,j)+2._RK*u(i,j)-u(i+1,j))*mesh%invdx2 &
                + (-u(i,j-1)+2._RK*u(i,j)-u(i,j+1))*mesh%invdy2
      end do
    end do

  end subroutine laplacian

end module laplacian_mod
