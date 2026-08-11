function red_get_psf,nx,ny,cx,cy
  psf=dblarr(nx,ny)
  nx2=double(nx/2)
  ny2=double(ny/2)
  cyy=cy/(2.d0*sqrt(2.d0*alog(2.d0)))
  cxx=cx/(2.d0*sqrt(2.d0*alog(2.d0)))
  for yy=0.,ny-1. do for xx=0.,nx-1 do begin
     u=((xx-nx2)/cxx)^2.+((yy-ny2)/cyy)^2.
     psf[xx,yy]=exp(-u*0.5d0)
  endfor
  return,psf
end
