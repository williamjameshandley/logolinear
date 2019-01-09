with(LinearAlgebra):

# Planck mass
mp:=1:

# Series indexing
step:=2/3:
jmax:=10:
pm := -1:
#b := 0:
#N_p :=0:
eta_p := 0:

# boolean conversion
bool := x -> if x then 1 else 0 end if:

`latex/special_names`[phi_p]:="\\phi_\\p":
`latex/special_names`[b0]:="b_0":
`latex/special_names`[b4]:="b_4":
`latex/special_names`[pm]:="pm":
`latex/special_names`[xx]:="(\\log t)":
coll := (expr,pow,func) -> sort(collect(expr,pow,func),pow,ascending):
simp := expr ->
    coll(expand(expr),t,
            x -> coll(x,xx,
                y -> coll(y, m,
                    z -> coll(z, phi_p,simplify)
                    )
                )
        ):

BellCoeff := proc()
if _nrest = 0  then
    1
else
    add(seq( p*BellCoeff(op([_rest][1.._nrest-p]))*[_rest][p],p=1.._nrest))/_nrest
end if
end proc:

DellCoeff := proc()
if _nrest = 0  then
    0
else
    [_rest][_nrest] - add(seq( (_nrest-p)*DellCoeff(op([_rest][1.._nrest-p]))*[_rest][p],p=1.._nrest-1))/_nrest
end if
end proc:


EellCoeff := proc(alpha)
if _nrest = 0  then
    1
else
    add(seq( ((alpha+1)*p-_nrest)*EellCoeff(alpha,op([_rest][1.._nrest-p]))*[_rest][p],p=1.._nrest))/_nrest
end if
end proc:


# Potential definition
V      := j -> 1/2 * m^2 *     add(seq(seq( phi[p]*phi[q]*bool(p+q=j-2),p=0..j-2,step),q=0..j-2,step)):
dVdphi := j ->       m^2 * phi[j-2] * bool(j-2=0):

# A_j definition
A:= j -> Matrix([
        [j,0,-1,0,0],
        [0,j,0,-1,0],
        [0,0,j-1+2*h[0],2/3/mp^2*v[0],0],
        [0,0,3*v[0],j-1+3*h[0],0],
        [0,0,0,0,j]
        ]):

# F_j definition
F := j -> Vector([
0,
0,
1/3/mp^2 * V(j)
-add(seq(seq((h[p]*h[q]+v[p]*v[q]/3/mp^2)*bool(p+q=j)*bool(q<>j)*bool(p<>j),p=0..2*j,step),q=0..2*j,step)),
-dVdphi(j)
-3*add(seq(seq(v[p]*h[q]*bool(p+q=j)*bool(p<>j)*bool(q<>j),p=0..2*j,step),q=0..2*j,step)),
exp(-N_p)*BellCoeff(seq(-N[q],q=step..j-step,step))*bool(j-step>=0)
]):


pmsubs := seq(pm^(2*j)=1,j=1..10*jmax),seq(pm^(2*j-1)=pm,j=1..10*jmax):

for j from 0 by step to jmax do
    if j = 0 then
        hjk := <N_p + 1/3*xx, phi_p+pm*sqrt(2/3)*mp*xx, 1/3, pm*sqrt(2/3)*mp, eta_p>:

    elif j = 4/3 then
        hjk := -9/14*<b, -pm*3/4*sqrt(6)*b*mp, 4/3*b, -pm*sqrt(6)*b*mp, 0>:

    else
        Fj := F(j):

        Nj := 0:
        for jj from 1 to 5 do
            ans := 'ans': lcoeff(Fj[jj],xx,`ans`):
            Nj := max(Nj,simplify(log(ans)/log(xx),symbolic)):
        end do:

        xjk := <0,0,0,0,0>:
        hjk := <0,0,0,0,0>:
        Ajinv := MatrixInverse(A(j)):
        for k from Nj+1 by -1 to 0 do
            hjk := hjk + xx^k * xjk:
            Fjk := coeff(Fj,xx,k-1):
            xjk := expand~(Multiply(Ajinv,Fjk-k*xjk)):
        end do:

    end if:
    hjk := subs(pmsubs,hjk):
    N[j], phi[j], h[j], v[j], eta[j] := op(convert(hjk,list)):
end do:

N_ser := jmax -> add(seq(N[j]*t^j, j=0..jmax,step)):
phi_ser := jmax -> add(seq(phi[j]*t^j, j=0..jmax,step)):


H_ser := jmax -> add(seq( (
collect(t*subs(log(t)=xx,diff(subs(xx=log(t),N[j]),t)) + j*N[j],xx)
)* t^(j-1), j=0..jmax,step)):

dphi_ser := jmax -> add(seq( (
collect(t*subs(log(t)=xx,diff(subs(xx=log(t),phi[j]),t)) + j*phi[j],xx)
)* t^(j-1), j=0..jmax,step)):

latex(simp(phi_ser(2)));
latex(simp(H_ser(2)));

truncate := (expr,k) -> convert(simplify(subs(x=t^(1/3),series(simplify(subs(t=x^3,expr),symbolic),x,k*3))),polynom):

log_ := (x, l) -> log(x[l]) + l * log(t) + add(seq(DellCoeff(seq(x[k+l]/x[l],k=step..j-step,step))* t^j,j=0..jmax,step)):
pow_ := (x, l, a) -> x[l]^a * add(seq(EellCoeff(a,seq(x[k+l]/x[l],k=step..j-step,step))* t^j,j=0..jmax,step)):

min_index := proc(x)
    local p,j:
    for j from 0 by step to jmax do
        if x[j] <> 0 then
            p[j]:=x[j]:
        end if:
    end do:
    return min(indices(p,`nolist`)):
end proc:

neaten := proc(x)
    local p:
    for p in [indices(x,`nolist`)] do
        x[p] := simp(x[p]);
    end do:
    x;
end proc:

log_ := proc(x)
    local j,l,ans:
    l := min_index(x):
    ans[0] := log(x[l]) + l * xx:
    for j from step by step to jmax do
        ans[j]:= DellCoeff(seq(x[k+l]/x[l],k=step..j-step,step)):
    end do:
    neaten(ans):
end proc:


exp_ := proc(a,x)
    local j,l,ans,p,q:
    l := min_index(x):
    if l = 0 then
        p := coeff(x[0],xx,1):
        q := coeff(x[0],xx,0):
        if x[0]-p*xx-q <>0 then
            print("=====================================WARNING=====================================")
        end if:
    else
        p := 0:
        q := 0:
    end if:

    for j from 0 by step to jmax-p*a do
        ans[j+p*a]:= exp(a*q)*BellCoeff(seq(a*x[k],k=step..j,step)):
    end do:
    neaten(ans):
end proc:

pow_ := proc(x,a)
    local j,l,ans:
    l := min_index(x):
    for j from l*a by step to jmax do
        ans[j]:= x[l]^a*EellCoeff(a,seq(x[k+l]/x[l],k=step..j-step,step)):
    end do:
    neaten(ans):
end proc:

mult_ := proc(x,y)
    local p,q,ans:
    for p in [indices(x,`nolist`)] do
        for q in [indices(y,`nolist`)] do
            if p+q < jmax then
                if assigned(ans[p+q]) then
                    ans[p+q] := ans[p+q] + x[p]*y[q]:
                else
                    ans[p+q] := x[p]*y[q]:
                end if:
            end if:
        end do:
    end do:
    neaten(ans):
end proc:

ser_ := (x,n) -> simp(add(seq(x[i]*t**i*bool(i<=n), i in indices(x,`nolist`)))):

eta_ := simp(ser_(eta,8/3)):
phi_ := simp(ser_(phi,2)):
dphi_ := simp(truncate(simp(subs(log(t)=xx,diff(subs(xx=log(t),ser_(phi,12/3)),t))),3)):
a_ := ser_(exp_(1,N),7/3):
h_ := simp(truncate(simp(
subs(log(t)=xx,diff(subs(xx=log(t),ser_(N,12/3)),t))
*
ser_(exp_(-1,N),5))
,3)):

1 - dphi_ser(2)^2/6/H_ser(2)^2;
series(%,t,4);
#simp(diff(phi_ser(2),t));

stop;
simp(H_ser(2));
simp(simp(a_*exp(-N_p)));
stop;
latex(simp(a_*exp(-N_p)));
latex(simp(eta_*exp(N_p)));

phi0 := 20.7:
m := 6e-6:
r := 100:

vars:={phi_p=evalf(phi0+sqrt(2/3)*log(sqrt(2/3)/(r*m*phi0))),t=sqrt(2/3)/(r*m*phi0), N_p=0}:
eqns:=subs(xx=log(t),[phi_ = 20, dphi_ = -r*m*phi0, a_=1]):

evalf(subs(vars,eqns)):
ans := fsolve(eqns,vars);

evalf(subs(xx=log(t),ans,eta_));
evalf(subs(xx=log(t),ans,1/(2*sqrt((r^2+1)/6)*m*phi0)));
evalf(subs(xx=log(t),ans,(eta_-1/(2*sqrt((r^2+1)/6)*m*phi0))/eta_));
