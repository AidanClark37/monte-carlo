module interpolation
  implicit none
  contains
subroutine interpolate(ydata,xval,yval)
  use struct_funcs
  use mpi_modules
  implicit none
  real*8,intent(in)::ydata(:)
  real*8,intent(in)::xval
  real*8,intent(out)::yval
  integer::i,ii
  real*8::xmin,xmax,xrange,slope,L0,L1,L2,h
  
  ii=1
  xmin=xarray(1)
  xmax=xarray(msg%npoint)
  xrange=xmax-xmin
  h=xrange/(msg%npoint-1)
!  write(*,*)xval," > ",xarray(msg%npoint)
  if (xval.gt.xarray(msg%npoint))then
     yval=0
  else
     
  ii = int((xval-xmin)/h)+1
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !linear interpolation
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 ! ii=min(ii,msg%npoint-1)
 ! slope=(ydata(ii+1)-ydata(ii))/(xarray(ii+1)-xarray(ii))

 ! yval = slope*(xval-xarray(ii))+ydata(ii)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !quadratic interpolation
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ii=min(ii,msg%npoint-2)
  L0=(xval-xarray(ii+1))*(xval-xarray(ii+2))/((xarray(ii)-xarray(ii+1))*(xarray(ii)-xarray(ii+2)))
  L1=(xval-xarray(ii))*(xval-xarray(ii+2))/((xarray(ii+1)-xarray(ii))*(xarray(ii+1)-xarray(ii+2)))
  L2=(xval-xarray(ii))*(xval-xarray(ii+1))/((xarray(ii+2)-xarray(ii))*(xarray(ii+2)-xarray(ii+1)))
    yval=ydata(ii)*L0 + ydata(ii+1)*L1+ydata(ii+2)*L2
end if
return
end subroutine interpolate
end module interpolation
