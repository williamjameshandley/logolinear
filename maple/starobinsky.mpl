with(LinearAlgebra):

# Planck mass
mp:=1:

# Series indexing
step:=2/3:
jmax:=10:

# Choose branch
pm:=-1:

# boolean conversion
bool := x -> if x then 1 else 0 end if:


BellCoeff := proc()
if _nrest = 0  then
    1
else
    add(seq( p*BellCoeff(op([_rest][1.._nrest-p]))*[_rest][p],p=1.._nrest))/_nrest
end if
end proc:

# Potential definition
if pm > 0 then
    V := j -> mp^2 * Lambda^4 * (
            bool(j-2=0) 
            -2*Phi_p*  BellCoeff(seq(-  sqrt(2/3)*phi[q]/mp,q=2/3..(j-4/3),step))*bool(j-4/3>=0)
            +  Phi_p^2*BellCoeff(seq(-2*sqrt(2/3)*phi[q]/mp,q=2/3..(j-2/3),step))*bool(j-2/3>=0)
            ):
    dVdphi := j -> mp * Lambda^4 * (
            +2*sqrt(2/3)*Phi_p  *BellCoeff(seq(-  sqrt(2/3)*phi[q]/mp,q=2/3..(j-4/3),step))*bool(j-4/3>=0)
            -2*sqrt(2/3)*Phi_p^2*BellCoeff(seq(-2*sqrt(2/3)*phi[q]/mp,q=2/3..(j-2/3),step))*bool(j-2/3>=0)
            ):         
else:
    V := j -> mp^2 * Lambda^4 * (
            bool(j-2=0) 
            -2*Phi_p*  BellCoeff(seq(-  sqrt(2/3)*phi[q]/mp,q=2/3..(j-8/3),step))*bool(j-8/3>=0)
            +  Phi_p^2*BellCoeff(seq(-2*sqrt(2/3)*phi[q]/mp,q=2/3..(j-10/3),step))*bool(j-10/3>=0)
            ):
    dVdphi := j -> mp * Lambda^4 * (
            +2*sqrt(2/3)*Phi_p  *BellCoeff(seq(-  sqrt(2/3)*phi[q]/mp,q=2/3..(j-8/3),step))*bool(j-8/3>=0)
            -2*sqrt(2/3)*Phi_p^2*BellCoeff(seq(-2*sqrt(2/3)*phi[q]/mp,q=2/3..(j-10/3),step))*bool(j-10/3>=0)
            ):         
end if:

# A_j definition
A:= j -> Matrix([
        [j,0,-1,0],
        [0,j,0,-1],
        [0,0,j-1+2*h[0],2/3/mp^2*v[0]],
        [0,0,3*v[0],j-1+3*h[0]]
        ]):

# F_j definition
F := j -> Vector([
0,
0,
1/3/mp^2 * V(j)
-add(seq(seq((h[p]*h[q]+v[p]*v[q]/3/mp^2)*bool(p+q=j)*bool(q<>j)*bool(p<>j),p=0..2*j,step),q=0..2*j,step)),
-dVdphi(j)
-3*add(seq(seq(v[p]*h[q]*bool(p+q=j)*bool(p<>j)*bool(q<>j),p=0..2*j,step),q=0..2*j,step))
]):


pmsubs := seq(pm^(2*j)=1,j=1..10*jmax),seq(pm^(2*j-1)=pm,j=1..10*jmax):

for j from 0 by step to jmax do
    if j = 0 then
        hjk := <N_p + 1/3*log(t), phi_p+pm*sqrt(2/3)*mp*log(t), 1/3, pm*sqrt(2/3)*mp>:

    elif j = 4/3 then
        if pm < 0 then:
            hjk := < -9*b/14, 27*sqrt(6)/56*b, -6*b/7, 9*sqrt(6)/14 * b >:
        else:
            hjk := < 
            -9*b/14 - 9 * Lambda^4*Phi_p/4 - 27*Lambda^8*Phi_p^4/50,
            27*sqrt(6)/56 * ( b + 2*Lambda^4*Phi_p/9 + 49*Lambda^8*Phi_p^4/300),
            -6*b/7 - 6 * Lambda^4*Phi_p/7 - 18*Lambda^8*Phi_p^4/25,
            9*sqrt(6)/14 * ( b + 2*Lambda^4*Phi_p/9 + 49*Lambda^8*Phi_p^4/300)
            >:
        end if:
        

    else
        Fj := F(j):

        Nj := 0:
        for jj from 1 to 4 do
            ans := 'ans': lcoeff(Fj[jj],log(t),`ans`):
            Nj := max(Nj,simplify(log(ans)/log(log(t)),symbolic)):
        end do:

        xjk := <0,0,0,0>:
        hjk := <0,0,0,0>:
        Ajinv := MatrixInverse(A(j)):
        for k from Nj+1 by -1 to 0 do
            hjk := hjk + log(t)^k * xjk:
            Fjk := coeff(Fj,log(t),k-1):
            xjk := expand~(Multiply(Ajinv,Fjk-k*xjk)):
        end do:

    end if:
    hjk := subs(pmsubs,hjk):
    N[j], phi[j], h[j], v[j] := op(convert(hjk,list)):
end do:

N_ser := jmax -> add(seq(N[j]*t^j, j=0..jmax,step)):
phi_ser := jmax -> add(seq(phi[j]*t^j, j=0..jmax,step)):

H_ser := jmax -> add(seq( (
collect(t*diff(N[j],t) + j*N[j],log(t))
)* t^(j-1), j=0..jmax,step)):


`latex/special_names`[phi_p]:="\\phi_\\p":
`latex/special_names`[Phi_p]:="\\Phi_\\p":
`latex/special_names`[b0]:="b_0":
`latex/special_names`[b4]:="b_4":
`latex/special_names`[pm]:="pm":
`latex/special_names`[xx]:="(\\log t)":
coll := (expr,pow,func) -> sort(collect(expr,pow,func),pow,ascending):
simp := expr ->
    coll(subs(log(t)=xx,expand(expr)),t,
            x -> coll(x,xx,
                y -> coll(y, [Lambda],
                    z -> coll(z, phi_p,
                        w -> coll(w, Phi_p, ww -> simplify(ww,symbolic))
                    )
                    )
                )
        ):

# Paper printing
latex(simp(H_ser(2+2*step)));
latex(simp(phi_ser(2+2*step)));
