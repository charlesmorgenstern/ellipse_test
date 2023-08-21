using Printf
using Plots
using QuadGK: quadgk
using Roots
using LinearAlgebra: norm, eigen, dot, cross
##########################################################
##########################################################
function getrhos(n)
temp=collect(range(0,4,n+1))
rhos=temp[2:end]
rhos=rhos.^2
return rhos
end
##########################################################
##########################################################
function getc(n)

f(x)=sqrt(16*sin(x)^2+16*(1/2)*cos(x)^2)
length, error = quadgk(f, 0, pi/2)
d=length/n

b=4*sqrt(.5)
a=4

tstep=pi/10000

c=Matrix{Float64}(undef,n,1)
c[1]=0
for i=2:n
 l=d*(i-1)
 t=tstep
 check=0
  while check<l
   check, error = quadgk(f,0,t)
   x=a*cos(t)
   y=b*sin(t)
   t+=tstep
   c[i]=y/x^2
  end
end

return c
end
##########################################################
##########################################################
function plotsetup(nrho,nc)

t=collect(range(0,pi/2,200))
x=collect(range(0,4,200))
y=0*x
fig=plot(x,y,legend=false)

rhos=getrhos(nrho)
for i=1:nrho
 x=sqrt(rhos[i])*cos.(t)
 y=sqrt(rhos[i])*sqrt(.5)*sin.(t)
 plot!(x,y,legend=false,xlims=(0,4),ylims=(0,3))
end

c=getc(nc)
x=collect(range(0,4,200))
for i=1:nc
y=c[i].*x.^2
plot!(x,y,legend=false,xlims=(0,4),ylims=(0,3))
end

pts=getintersections(nrho,nc)
for i=1:nrho
scatter!(pts[:,2*i-1],pts[:,2*i],legend=false,markersize=3)
end

return fig
end
##########################################################
##########################################################
function getintersections(nrho,nc)

npts=nrho*nc
pts=Matrix{Float64}(undef,nc,nrho*2)
rhos=getrhos(nrho) 
c=getc(nc)  
c=reverse(c)
for i=1:nc #loop through parabolas
 for j=1:nrho #loop through ellipses
 d=c[i]
 r=rhos[j]
 f(y)=2*d*y^2+y-r*d
 y=fzero(f,1.5)
 x=sqrt(y/d)
if y==0
x=sqrt(r)
end
 pts[i,2*j]=y
 pts[i,2*j-1]=x
 end
end

pts=reverse(pts,dims=1)
return pts
end
##########################################################
##########################################################
function getparameters(nrho,nc)
pts=getintersections(nrho,nc)
t=Matrix{Float64}(undef,nc+1,nrho)
t[1,:].=0.0
t[end,:].=pi/2
rhos=getrhos(nrho)
r=sqrt.(rhos)
for i=2:nc #loop parabolas
for j=1:nrho #loop ellipses
r1=r[j]
x=pts[i,2*j-1]
t[i,j]=acos(x/r1)
end
end
return t
end
##########################################################
##########################################################
function getarclengths(nrho,nc)
arcs=Matrix{Float64}(undef,nc,nrho)
rhos=getrhos(nrho) 
t=getparameters(nrho,nc)
for i=1:nrho #loop over ellipses
for j=1:nc #loop over segments

t1=t[j,i]
t2=t[j+1,i]
f(x)=sqrt(rhos[i]*sin(x)^2+rhos[i]*(1/2)*cos(x)^2)
integral, error = quadgk(f, t1, t2)
arcs[j,i]=integral
end
end

return arcs
end
##########################################################
##########################################################
function getgplengths(nrho,nc)
pts=getintersections(nrho,nc)
ds=Matrix{Float64}(undef,nc,nrho-1)
c=getc(nc)


for i=1:nc #loop through parabolas
for j=1:nrho-1 #loop through segments 
d=c[i]
x1=pts[i,2*j-1]
x2=pts[i,2*j+1]
f(x)=sqrt(1+4*d^2*x^2)
integral, error = quadgk(f, x1, x2)
ds[i,j]=integral
end
end

return ds
end
##########################################################
##########################################################
function getcurvatures(nrho,nc)
pts=getintersections(nrho,nc)
rhos=getrhos(nrho) 
curv=Matrix{Float64}(undef,nc,nrho)

for i=1:nc #loop over parabolas
for j=1:nrho #loop over ellipse intersections
x=pts[i,2*j-1]
y=pts[i,2*j]
a=sqrt(rhos[j])
b=a*sqrt(.5)
curv[i,j]=(1/(a^2*b^2))*((x^2/a^4)+(y^2/b^4))^(-3/2)
end
end

return curv
end
##########################################################
##########################################################
function getcurvderivative(nrho,nc)
curv=getcurvatures(nrho,nc)
arcs=getarclengths(nrho,nc)
dc=Matrix{Float64}(undef,nc,nrho)
t=getparameters(nrho,nc)
tstep=pi/100000000
rhos=getrhos(nrho) 

for i=1:nrho #loop over ellipses
for j=1:nc #loop over segments
c1=curv[j,i]
t1=t[j,i]
t2=t1+tstep
a=sqrt(rhos[i])
b=a*sqrt(.5)
x2=a*cos(t2)
y2=b*sin(t2)
c2=(1/(a^2*b^2))*((x2^2/a^4)+(y2^2/b^4))^(-3/2)
f(s)=sqrt(a^2*sin(s)^2+b^2*cos(s)^2)
dr, error = quadgk(f, t1, t2)
dc[j,i]=(c2-c1)/dr
end
end

return dc
end
##########################################################
##########################################################
function getarclengthsgponly(nrho,nc)
arcs=Matrix{Float64}(undef,nc,nrho)
curv=getcurvatures(nrho,nc)
dc=getcurvderivative(nrho,nc)
dr=getgplengths(nrho,nc)
arcs1=getarclengths(nrho,nc)
arcs[:,1]=arcs1[:,1]

for i=2:nrho #loop over ellipses
for j=1:nc #loop over segments
k=curv[j,i-1]
dk=dc[j,i-1]
li=arcs[j,i-1]
l0=arcs[j,1]
ds=dr[j,i-1]
arcs[j,i]=li*(1+k*ds+dk/2*l0*ds)
end
end

return arcs
end
##########################################################
##########################################################
function getareatrap(nrho,nc)
dr=getgplengths(nrho,nc)
arcs=getarclengths(nrho,nc)
rhos=getrhos(nrho)
a=sqrt(rhos[1])
b=a*sqrt(.5)
area=pi*a*b/4

for i=1:nrho-1 #loop down gradient path segments
for j=1:nc #loop around ellipse segments
area+=dr[j,i]*(arcs[j,i]+arcs[j,i+1])/2
end
end
return area
end
##########################################################
##########################################################
function checkareatrap(nrho,nc)
area=getareatrap(nrho,nc)
exact=pi*4*sqrt(1/2)
error=abs(exact-area)/exact
return error
end
##########################################################
##########################################################
function getareagponlytrap(nrho,nc)
dr=getgplengths(nrho,nc)
arcs=getarclengthsgponly(nrho,nc)
rhos=getrhos(nrho)
a=sqrt(rhos[1])
b=a*sqrt(.5)
area=pi*a*b/4

for i=1:nrho-1 #loop down gradient path segments
for j=1:nc #loop around ellipse segments
area+=dr[j,i]*(arcs[j,i]+arcs[j,i+1])/2
end
end
return area
end
##########################################################
##########################################################
function checkareagponlytrap(nrho,nc)
area=getareagponlytrap(nrho,nc)
exact=pi*4*sqrt(1/2)
error=abs(exact-area)/exact
return error
end
##########################################################
##########################################################
function checkarclengths(nrho,nc)
err=Matrix{Float64}(undef,nrho,1)
arcs=getarclengths(nrho,nc)
rhos=getrhos(nrho) 

for i=1:nrho #loop through ellipses
a=sqrt(rhos[i])
b=a*sqrt(.5)
f(t)=a*sqrt(1-(1-b^2/a^2)*sin(t)^2)
exact, error = quadgk(f, 0, pi/2)
len=0.0
 for j=1:nc
len+=arcs[j,i]
 end
err[i]=abs(exact-len)/exact
end
return err
end
##########################################################
##########################################################
function checkarclengthsgponly(nrho,nc)
err=Matrix{Float64}(undef,nrho,1)
arcs=getarclengthsgponly(nrho,nc)
rhos=getrhos(nrho)

for i=1:nrho #loop through ellipses
a=sqrt(rhos[i])
b=a*sqrt(.5)
f(t)=a*sqrt(1-(1-b^2/a^2)*sin(t)^2)
exact, error = quadgk(f, 0, pi/2)
len=0.0
 for j=1:nc
len+=arcs[j,i]
 end
err[i]=abs(exact-len)/exact
end
return err
end
##########################################################
##########################################################
function checkarclengthserr(nrho,nc)
arcs=getarclengthsgponly(nrho,nc)
arcs2=getarclengths(nrho,nc)
err=abs.(arcs-arcs2)./abs.(arcs2)
return err
end
##########################################################
##########################################################
function maxarclengtherr(nrho,nc)
err=checkarclengthserr(nrho,nc)
max=maximum(err)
return max
end
##########################################################
##########################################################
function errortablearea()
n=(5,10,20,40,80,160,320)
err=Matrix{Float64}(undef,7,1)
eoc=Matrix{Float64}(undef,7,1)
err[1]=checkareatrap(n[1],n[1])
for i=2:7
err[i]=checkareatrap(n[i],n[i])
eoc[i]=log(err[i-1]/err[i])/log(2)
end
@printf "\n Relative error of ellipse area using reference arc lengths with n ellipses and n parabolas"
@printf "\n n      rel. error            eoc"
@printf "\n %g       %g           n/a" n[1] err[1]
for i=2:7
@printf "\n %g       %g            %g" n[i] err[i] eoc[i] 
end
end
##########################################################
##########################################################
function errortableareagponly()
n=(5,10,20,40,80,160,320)
err=Matrix{Float64}(undef,7,1)
eoc=Matrix{Float64}(undef,7,1)
err[1]=checkareatrap(n[1],n[1])
for i=2:7
err[i]=checkareagponlytrap(n[i],n[i])
eoc[i]=log(err[i-1]/err[i])/log(2)
end
@printf "\n Relative error of ellipse area using recurrence relation with n ellipses and n parabolas"
@printf "\n n      rel. error            eoc"
@printf "\n %g       %g           n/a" n[1] err[1]
for i=2:7
@printf "\n %g       %g            %g" n[i] err[i] eoc[i]
end
end
##########################################################
##########################################################
function errortablearclength()
n=(5,10,20,40,80,160,320)
err=Matrix{Float64}(undef,7,1)
eoc=Matrix{Float64}(undef,7,1)
err[1]=maxarclengtherr(n[1],n[1])
for i=2:7
err[i]=maxarclengtherr(n[i],n[i])
eoc[i]=log(err[i-1]/err[i])/log(2)
end
@printf "\n Max relative error of arc lengths with n ellipses and n parabolas"
@printf "\n n      rel. error            eoc"
@printf "\n %g       %g           n/a" n[1] err[1]
for i=2:7
@printf "\n %g       %g            %g" n[i] err[i] eoc[i]
end
end
##########################################################
##########################################################
function animategrids()
n=(5,10,20,40,80,160)

anim = @animate for i ∈ 1:6
plotsetup(n[i],n[i])
end
gif(anim, "grids_fps5.gif", fps = 5)
end
