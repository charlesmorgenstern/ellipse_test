using Plots
using QuadGK: quadgk
using Roots
using LinearAlgebra: norm, eigen, dot, cross
##########################################################
##########################################################
function plotsetup()
t=collect(range(0,pi/2,200))
x=collect(range(0,4,200))
y=0*x
fig=plot(x,y,legend=false)
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)
n=length(rhos)
for i=1:n
x=sqrt(rhos[i])*cos.(t)
y=sqrt(rhos[i])*sqrt(.5)*sin.(t)
plot!(x,y,legend=false,xlims=(0,4),ylims=(0,3))
end
c=  (9.63132, 2.36291, 1.01601, 0.543444, 0.323286, 0.201919, 0.126465, 0.0744632, 0.0345981)
n=length(c)
x=collect(range(0,4,200))
for i=1:n
y=c[i].*x.^2
plot!(x,y,legend=false,xlims=(0,4),ylims=(0,3))
end
pts=getintersections()
for i=1:10
scatter!(pts[:,2*i-1],pts[:,2*i],legend=false)
end

return fig
end
##########################################################
##########################################################
function getintersections()
pts=Matrix{Float64}(undef,10,20)
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)
c=  (9.63132, 2.36291, 1.01601, 0.543444, 0.323286, 0.201919, 0.126465, 0.0744632, 0.0345981, 0.0)
c=reverse(c)
for i=1:10 #loop through parabolas
 for j=1:10 #loop through ellipses
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

return pts
end
##########################################################
##########################################################
function getparameters()
pts=getintersections()
t=Matrix{Float64}(undef,11,10)
t[1,:].=0.0
t[11,:].=pi/2
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)
r=sqrt.(rhos)
for i=2:10 #loop parabolas
for j=1:10 #loop ellipses
r1=r[j]
x=pts[i,2*j-1]
t[i,j]=acos(x/r1)
end
end
return t
end
##########################################################
##########################################################
function getarclengths()
arcs=Matrix{Float64}(undef,10,10)
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)
t=getparameters()
for i=1:10 #loop over ellipses
for j=1:10 #loop over segments

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
function getgplengths()
pts=getintersections()
ds=Matrix{Float64}(undef,10,9)
c=  (9.63132, 2.36291, 1.01601, 0.543444, 0.323286, 0.201919, 0.126465, 0.0744632, 0.0345981, 0.0)
c=reverse(c)

for i=1:10 #loop through parabolas
for j=1:9 #loop through segments 
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
function getcurvatures()
pts=getintersections()
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)
curv=Matrix{Float64}(undef,10,10)

for i=1:10 #loop over parabolas
for j=1:10 #loop over ellipse intersections
x=pts[i,2*j-1]
y=pts[i,2*j]
a=sqrt(rhos[j])
b=a*sqrt(.5)
curv[i,j]=(1/(a^2*b^2))*((x^2/a^4)+(y^2/b^4))^(-3/2)
end
end

#t=pi/2
#for i=1:10 #final curvature along y axis
#a=sqrt(rhos[i])
#b=a*sqrt(.5)
#x=a*cos(t)
#y=b*sin(t)
#curv[11,i]=(1/(a^2*b^2))*((x^2/a^4)+(y^2/b^4))^(-3/2)
#end

return curv
end
##########################################################
##########################################################
function getcurvderivative()
curv=getcurvatures()
arcs=getarclengths()
dc=Matrix{Float64}(undef,10,10)
t=getparameters()
tstep=pi/100000000
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)

for i=1:10 #loop over ellipses
for j=1:10 #loop over segments
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
function getarclengthsgponly()
arcs=Matrix{Float64}(undef,10,10)
curv=getcurvatures()
dc=getcurvderivative()
dr=getgplengths()
arcs1=getarclengths()
arcs[:,1]=arcs1[:,1]

for i=2:10 #loop over ellipses
for j=1:10 #loop over segments
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
function getareatrap()
dr=getgplengths()
arcs=getarclengths()
area=pi*(1/25)*sqrt(1/2)

for i=1:9 #loop down gradient path segments
for j=1:10 #loop around ellipse segments
area+=dr[j,i]*(arcs[j,i]+arcs[j,i+1])/2
end
end
return area
end
##########################################################
##########################################################
function checkareatrap()
area=getareatrap()
exact=pi*4*sqrt(1/2)
error=abs(exact-area)/exact
return error
end
##########################################################
##########################################################
function getareagponlytrap()
dr=getgplengths()
arcs=getarclengthsgponly()
area=pi*(4/25)*sqrt(1/2)/4

for i=1:9 #loop down gradient path segments
for j=1:10 #loop around ellipse segments
area+=dr[j,i]*(arcs[j,i]+arcs[j,i+1])/2
end
end
return area
end
##########################################################
##########################################################
function checkareagponlytrap()
area=getareagponlytrap()
exact=pi*4*sqrt(1/2)
error=abs(exact-area)/exact
return error
end
##########################################################
##########################################################
function checkarclengths()
err=Matrix{Float64}(undef,10,1)
arcs=getarclengths()
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)

for i=1:10 #loop through ellipses
a=sqrt(rhos[i])
b=a*sqrt(.5)
f(t)=a*sqrt(1-(1-b^2/a^2)*sin(t)^2)
exact, error = quadgk(f, 0, pi/2)
len=0.0
 for j=1:10
len+=arcs[j,i]
 end
err[i]=abs(exact-len)/exact
end
return err
end
##########################################################
##########################################################
function checkarclengthsgponly()
err=Matrix{Float64}(undef,10,1)
arcs=getarclengthsgponly()
rhos = (4/25, 16/25, 36/25, 64/25, 4, 144/25, 196/25, 256/25, 324/25, 16)

for i=1:10 #loop through ellipses
a=sqrt(rhos[i])
b=a*sqrt(.5)
f(t)=a*sqrt(1-(1-b^2/a^2)*sin(t)^2)
exact, error = quadgk(f, 0, pi/2)
len=0.0
 for j=1:10
len+=arcs[j,i]
 end
err[i]=abs(exact-len)/exact
end
return err
end
##########################################################
##########################################################
function checkarclengthserr()
arcs=getarclengthsgponly()
arcs2=getarclengths()
err=abs.(arcs-arcs2)./abs.(arcs2)
return err
end
