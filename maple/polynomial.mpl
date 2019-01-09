with(LinearAlgebra):

# Planck mass
mp:=1:

# Series indexing
step:=2/3:
jmax:=10:

# boolean conversion
bool := x -> if x then 1 else 0 end if:

# Potential definition
V      := j -> 
1/2 * m^2 *     add(seq(seq(
                    phi[p]*phi[q]*bool(p+q=j-2),
                    p=0..j-2,step),q=0..j-2,step)) + 
1/24 * lambda * add(seq(seq(seq(seq(
                    phi[p]*phi[q]*phi[r]*phi[s]*bool(p+q+r+s=j-2),
                    p=0..j-2,step),q=0..j-2,step),r=0..j-2,step),s=0..j-2,step))
+ mp^2 * Lambda * bool(j-2=0):
dVdphi := j ->       m^2 * phi[j-2] * bool(j-2=0)+ 
1/6  * lambda * add(seq(seq(seq(
                    phi[p]*phi[q]*phi[r]*bool(p+q+r=j-2),
                    p=0..j-2,step),q=0..j-2,step),r=0..j-2,step)):

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
        hjk := -9/14*<b, -pm*3/4*sqrt(6)*b*mp, 4/3*b, -pm*sqrt(6)*b*mp>:

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
`latex/special_names`[b0]:="b_0":
`latex/special_names`[b4]:="b_4":
`latex/special_names`[pm]:="pm":
`latex/special_names`[xx]:="(\\log t)":
coll := (expr,pow,func) -> sort(collect(expr,pow,func),pow,ascending):
simp := expr ->
    coll(subs(log(t)=xx,expand(expr)),t,
            x -> coll(x,xx,
                y -> coll(y, [lambda,m],
                    z -> coll(z, phi_p,simplify)
                    )
                )
        ):

# Paper printing
latex(simp(H_ser(2+2*step)));
latex(simp(phi_ser(2+2*step)));

#las_simp := expr ->
#    coll(subs(b=-7/6 * 16/27/mp*sqrt(3/2)*b4, phi_p=b0,m=mu,lambda=0,log(t)=xx,expand(expr)),t,
#            x -> coll(x,xx,
#                y -> coll(y, [mu,b0,b4],simplify)
#                )
#        ):

# Lasenby doran printing
#las_simp(H_ser(2+1*step));
#las_simp(phi_ser(2+1*step));
