module prec_mod
  use, intrinsic :: ISO_Fortran_env, only: int32,real32,real64
  implicit none
#ifdef SINGLE_PRECISION
  integer,parameter,public :: RK = real32
#else
  integer,parameter,public :: RK = real64
#endif
  integer,parameter,public :: IK = int32
end module prec_mod
