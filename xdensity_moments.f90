program xdensity_moments
! 06/15/2021 10:02 PM computes moments for probability distributions
! 06/15/2021 06:09 PM branched from xxstudent_t_moments.f90 to xdensity_moments.f90
! 06/12/2021 09:14 PM specifieds standard deviation if set_sd is .true.
! 06/11/2021 09:15 PM branched from xstudent_t_moments.f90 to xxstudent_t_moments.f90
! 05/25/2014 09:17 AM computes moments of raw and standardized Student t distributions
! 05/25/2014 08:00 AM branched from xdensity.f90
! compute the moments of some distributions
use student_t_mod, only: student_t_density,tstd
use kind_mod     , only: dp
use pdf_mod      , only: std_density
implicit none
integer          , parameter :: n = 4000001,nmom=4,ndist = 5
real(kind=dp)                :: xx,xh,xmom(0:nmom),xmin,xmax,dens
integer                      :: i,imom,idist
character (len=*), parameter :: fmt_cr = "(a10,':',100f12.6)",fmt_ci = "(a10,':',100i12)"
character (len=20)           :: dist(ndist)
dist = [character (len=7) :: "normal","Laplace","sech","t5","t10"]
xmax = 100.0_dp
xmin = -xmax
xh   = (xmax-xmin)/(n-1)
write (*,"(/,12x,a12,/,a12,100i12)") "moment","density",(imom,imom=0,nmom)
do idist=1,ndist
   xmom = 0.0_dp
   do i=1,n
      xx          = xmin + (i-1)*xh
      dens        = std_density(xx,dist(idist))
      do imom=0,nmom
         if (imom == 0) then
            xmom(0)    = xmom(0)   + xh*dens
         else
            xmom(imom) = xmom(imom) + xh*(xx**imom)*dens
         end if
      end do
   end do
   write (*,"(a12,100f12.6)") trim(dist(idist)),xmom
end do
write (*,*)
write (*,fmt_ci) "n",n
write (*,fmt_cr) "xmin",xmin
write (*,fmt_cr) "xmax",xmax
write (*,fmt_cr) "xh",xh
write (*,"(/,a)") "(2) finished xdensity_moments.f90"
end program xdensity_moments
