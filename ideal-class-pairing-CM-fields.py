#!/usr/bin/env python
# coding: utf-8

# # Ideal Class Pairing Over Totally Real Fields

# --- 
# 
# ## 1. Totally Real Fields

# The function `enumerate_totallyreal_fields_all(d, x,return_pari_objects=False)` returns a list of totally real fields of degree *deg* and discriminant <= *disc*.
# 
# **Notes:**
# 
# - Each element of the list is a pair (*D*,*p*) where *D* is the discriminant of the field and *p* is the defining polynomial.
# 
# - The argument `return_pari_objects=False` is used so that the Sage recognizes the polynomial as a polynomial.
# 
# - If the argument `return_seqs=True` is used, then the polynomial is given as a sequence of coefficients.
# 
# For example, there are 5 totally real fields of degree 2 of discriminant <= 20:

# In[62]:


enumerate_totallyreal_fields_all(2,20,return_pari_objects=False)


# The function `trf_list_bounded_deg(d,D,t)` creates a list of all totally real fields of degree up to *d* and discriminant up to *D* ; the name of the primitive element is given by 't'. The function takes the following three arguments:
# - deg_bound, a positive integer that bounds the degrees of the fields,
# - disc_bound, a positive integer that bounds the discriminant of the fields,
# - 't', a string that names the primitive element.

# In[63]:


def trf_list_bounded_deg(deg_bound,disc_bound,adjoined_root_name):
    unflat_list_of_fields= []
    for d in (1..deg_bound):
        trf_enumerated_list = enumerate_totallyreal_fields_all(d,disc_bound,return_pari_objects=False)
        trf_poly_list = [y[1] for y in trf_enumerated_list]
        unflat_list_of_fields.append([NumberField(f,adjoined_root_name) for f in trf_poly_list])
    flat_list_of_fields = []
    for deg_d_fields in unflat_list_of_fields:
        for field in deg_d_fields:
            flat_list_of_fields.append(field)
    return flat_list_of_fields


# For example

# In[64]:


deg_bound = 3
disc_bound = 50
adjoined_root_name = 't'

trf_list_bounded_deg(deg_bound,disc_bound,adjoined_root_name)


# We require our totally real field to be of class number 1. So we select those from the `trf_list_bounded_deg(deg_bound,disc_bound,adjoined_root_name)` with class number 1

# In[65]:


deg_bound = 3
disc_bound = 50
adjoined_root_name = 't'
trf_full_list = trf_list_bounded_deg(deg_bound,disc_bound,adjoined_root_name)

class_1_trf_list = [F for F in trf_full_list if F.class_number()==1];
class_1_trf_list


# To choose a class number 1 totally real field, we look at `class_1_trf_list` and write

# In[66]:


def_poly = x^2 - 7
F.<t> = NumberField(def_poly)
F


# For example, we compute several arithmetic invariants of the chosen field:

# In[67]:


def_poly = x^2 - 7
F.<t> = NumberField(def_poly)

OF = F.maximal_order()
print('ring of integers has integral basis =',OF.gens())

UF = UnitGroup(F);
print('Unit group =',UF)
print('has fundamental units =', UF.fundamental_units())


# ---
# 
# ## 2. Elliptic Curves Over Totally Real Fields

# ### 2.1. Defining elliptic curves and computing points on them

# Elliptic curves and their points can be constructed as follows:
# 
# 1. The function `EllipticCurve(F,[A,B])` defines an elliptic curve over a field *F* with short Weierstrass equation *y^2 = x^3 + a_4 x + a_6* with coefficients *a_4* and *a_6* in *F*.
# 2. We can compute the torsion subgroup tor(*E* ( *F* )) with `torsion_points()`.
# 3. With `gens()`, we can compute "some" generators of the Mordell-Weil group of *E* over the field *F*. (*Note*: The function `gens()` will not always give you a complete generating set and/or linearly independent points. However, for the purposes of ideal class pairings, we are interested in just the existence of nontorsion points.)
# 
# For example:

# In[68]:


#1. First define the totally real field
def_poly = x^2 - 7
F.<t> = NumberField(def_poly)
print('- F =',F)

#2. Define the elliptic curve E
a4 = -2
a6 = 4
E = EllipticCurve(F,[a4,a6]);
print('- E =',E)

#3. The torsion points of E are:
print('- Torsion points of E:',E.torsion_points())

#3. The generators of the nontorsion part of E are:
print('- Some independent nontorsion points of E:')

E.gens()


# ### 2.2. Manipulating the coordinates of points on elliptic curves

# The coordinates of a point *P* in *E(F)* aren't necessarily given in lowest terms (recall that *F* is taken to be of class number 1). The function `reduce_coord` takes a point on *E* and returns a list (*x,y*) where *x* and *y* are reduced fractions in *F*. More precisely, it takes as input
# 1. trf = a totally real field of class number 1
# 2. coord = a list [x,y] with two elements that are the coordinates of a point in *E(F)*
# 
# and returns a list [[*a*,*c*],[*b*,*d*]] where *x*=*a*/*c* and *y*=*b*/*d* are reduced fractions, that is gcd(*a*,*c*)=gcd(*b*,*d*)=1.

# In[69]:


def reduce_coord(trf,coord):
    O_trf = trf.maximal_order()
    
    x_num = coord[0].numerator()
    x_den = coord[0].denominator()
    y_num = coord[1].numerator()
    y_den = coord[1].denominator()
    
    x_gcd = O_trf(x_num).gcd(x_den)
    y_gcd = O_trf(y_num).gcd(y_den)
    
    return [[x_num/x_gcd,x_den/x_gcd],[y_num/y_gcd,y_den/y_gcd]]


# For example

# In[70]:


#1. First define the totally real field
def_poly = x^2 - 7
F.<t> = NumberField(def_poly)
#F = QQ
print('- F =',F)

#2. Define the elliptic curve E
a4 = -2
a6 = 4
E = EllipticCurve(F,[a4,a6]);
print('- E =',E)

#3. Find nontorsion point on E
P = E.gens()[0]
coord = list((P+P).xy())
print('- Let P be the point P =',coord)

#4. Reduce coordinates to reduced fractions
print('- Then P = (a/c,b/d), written in reduced fractions, with [[a,c],[b,d]] =')
reduce_coord(F,coord)


# Another example (over **Q**):

# In[71]:


#1. Define the elliptic curve E
a4 = -2
a6 = 5
E = EllipticCurve(QQ,[a4,a6]);
print('- E =',E)

#3. Find nontorsion point on E
P = E.gens()[0]
coord = list((P+P+P).xy())
print('- Take the nontorsion point: P =',coord)

#4. Reduce coordinates to reduced fractions
print('Then P = (a/c,b/d), written in reduced fractions, with [[a,b],[c,d]] =')
reduce_coord(F,coord)


# The coordinates of a point in *E* ( *F* ) can be written in the form ( *a* / *c*^2 , *b* / *c*^3 ) where *a,b,c* are elements of the ring of integers of *F*. The function `standard_form` takes in as input
# 1. trf = totally real field of class number 1
# 2. coord = a list [x,y] with two elements that are the coordinates of a point *P* in *E(F)*
# 
# and returns a list [a,b,c] so that *P* = ( *a* / *c*^2 , *b* / *c*^3 ).

# In[72]:


def standard_form(trf,coord):
    red = reduce_coord(F,coord)  
    x_num = red[0][0] 
    y_num = red[1][0]
    e = red[1][1] / red[0][1]
    return [x_num,y_num,e]


# For example

# In[73]:


#1. First define the totally real field
def_poly = x^2 - 7
F.<t> = NumberField(def_poly)
print('- F =',F)

#2. Define the elliptic curve E
a4 = -2
a6 = 4
E = EllipticCurve(F,[a4,a6]);
print('- E =',E)

#3. Find nontorsion point on E
P = E.gens()[0].xy()
coord = list(P)
print('- Take the nontorsion point: P =',coord)
print('- The reduced form of P is (A/C,B/D) where [[A,C],[B,D]] =',reduce_coord(F,coord))

#4. Convert coordinates to standard form
print('- Then the standard form of P is (a/c^2,b/c^3), where [a,b,c] =')
standard_form(F,coord)


# ### 2.3. Twisting elliptic curves

# We can *twist* an elliptic curve *E* over *F*. The function `twist` takes as input:
# 
# 1. trf = a totally real field of class number 1,
# 2. a4, a6 = the coefficients of the short Weierstrass equation for E,
# 3. disc = a negative fundamental discriminant,
# 
# and returns the elliptic curve over *trf* with short Weierstrass equation
# >*y*^2 = *x*^3 + *a_4* *disc*^2 *x* + *a_6* *disc*^3.

# In[74]:


def twist(trf,A,B,disc):
    return EllipticCurve(trf,[A*disc^2,B*disc^3])


# For example, we can compute (possibly not all) generators of Mordell-Weil group of a twist.

# In[75]:


#1. Define the totally real field
def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
print('- F =',F)

#2. Define the elliptic curve E
a4 = -2
a6 = 5
E = EllipticCurve(F,[a4,a6]);
print('- E =',E)

#3.Twist E by d to get E_twist
disc = -15
print('The elliptic curve E twisted by disc =',disc,'is given by:')

E_twist = twist(F,a4,a6,disc);
E_twist


# The twisted elliptic curve defined by `twist` is the "canonical model" of the twist. However, we wish to change to another non-canonical model. In particular, we are interested in how the coordinates of the points change.
# 
# The function `nc_twist_coord` takes as input:
# 1. coord = A list of coordinates [x,y] of a point P on the canonical model of the twist of E given by
# >*y*^2 = *x*^3 + *a_4* *disc*^2 *x* + *a_6* *disc*^3
# 2. disc = a negative fundamental discriminant,
# 
# and returns a list [x',y'] which are the coordinates of the corresponding point P' on the non-canonical model of the twist defined by
# > *disc* (*y'*/2)^2 = *x'*^3 + *a_4* *x'* + *a_6*.
# 
# The change of variables is *x'=x/disc* and *y'=2y/disc^2*.

# In[76]:


def nc_twist_coord(coord,disc):
    x_coord = coord[0]
    y_coord = coord[1]
    new_x_coord = x_coord/disc
    new_y_coord = (2*y_coord)/(disc^2)
    return [new_x_coord,new_y_coord]


# For example

# In[77]:


#1. Define F and the elliptic curve E
def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
a4 = -2
a6 = 5
E = EllipticCurve(F,[a4,a6]);
print('- E =',E)

#2.Twist E by disc to get E_twist
disc = -15
E_twist = twist(F,a4,a6,disc);
print('- Twist of E by',disc,' is the curve')
print('  E_twist =',E_twist)

#3. Select a point P_twist on E_twist
P_twist = E_twist.gens()[0];
print('- Take some nontorsion point on E_twist:',P_twist.xy())

#4. Find coordinates of P_twist in the non canonical model
P_nc = nc_twist_coord(P_twist,disc)
print('- The corresponding point on the noncanonical model is:',vector(P_nc))
print('- It has standard form (a/c^2,b/c^3) with [a,b,c] =',standard_form(F,P_nc))


# ---
# 
# ## 3. CM Fields

# CM fields are totally imaginary quadratic extensions of totally real fields. A CM field *K* can be written in the form
# > *K* = *F*(*a*)    where *a*^2 is in *F* and *f*(a^2)<0 for every embedding *f* of *K* into **C**
# 
# Below is a list of negative fundamental discriminants whose square roots we can add to *F*.

# In[78]:


[-n for n in (1..50) if is_fundamental_discriminant(-n)==1]


# For the special type of discrimintant $-D_E(t):=-4(t^3 + a_4 t + a_6)$, we find a $t$ for which $-D_E(t)$ is fundamental:

# In[79]:


a4 = -2
a6 = 4

n = 1
while is_fundamental_discriminant(-4*(n^3 + a4*n + a6))==0:
    n += 1
print('first t for which -D_E(t) is fundamental is t =',n)
print('which corresponds to -D_E(t) =',-4*(n^3 + a4*n + a6))
print('indeed, is_fundamental_discriminant(',-4*(n^3 + a4*n + a6),') =',is_fundamental_discriminant(-4*(n^3 + a4*n + a6)))


# Given a totally real field *F* as above and a discriminant *disc* < 0, we construct an imaginary quadratic extension *K* of *F*.
# 
# *Notation*:
# - sqrt(*disc*) is denoted by *T*.
# - The relative quadratic extension is denoted by *Kr* so *Kr* = *F* ( *T* )
# - The absolute extension is denoted by *K* so *K* = Q( *S* ), were *S* is some generator.

# In[80]:


#1. Define F
def_poly = x^2 - 7
F.<t> = NumberField(def_poly)
print('- F =',F)

#2. Choose a negative fundamental discriminant and its minimal polynomial
D = -5252
def_poly_CM = x^2 - D

#3. The associated CM field is
Kr.<T> = F.extension(def_poly_CM)
print('- The CM-field as a relative extension is: Kr =',Kr)

K = Kr.absolute_field('S');
print('- The CM-field as an absolute extension is: K =',K)

#4. Class number
print('- Notice that K is indeed a CM-field: K.is_CM() =',K.is_CM())
print('- The class number of K is:')

K.class_number()


# ---
# 
# ## 4. Quadratic Forms Over Totally Real Fields

# Given a totally real field, we determine its ring of integers. This is the ring of coefficient of our binary quadratic forms.

# In[81]:


def_poly = x^2 - 3
F.<t> = NumberField(def_poly)

OF = F.maximal_order();
OF


# We can compute an integral basis of the ring of integers with `.basis()`

# In[82]:


def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
OF = F.maximal_order();

OF.basis()


# A binary quadratic form *a x*^2 + *b x y* + *c y*^2 over the ring of integers of a totally real field *F* is constructed with the function `QuadraticForm(O,2,[a,b,c])`. The polynomial version is obtained by then using `.polynomial()`.

# In[83]:


#1. Define F and its ring of integers OF
def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
OF = F.ring_of_integers();

#2. Write the coefficients of the form ax^2+bxy+cy^2 
a = 2
b = 2
c = 1+t

#3. Define the associated form with the above coefficients
Q = QuadraticForm(OF,2,[a,b,c]);
print('The quadratic form is:',Q)
print('As a polynomial, this form is:')
Q.polynomial()


# The discriminant of a form *Q* can be computed with `.disc()`

# In[84]:


def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
OF = F.ring_of_integers();

a = 2
b = 2
c = 1+t
Q = QuadraticForm(OF,2,[a,b,c]);

print('The discriminant of ',Q.polynomial(),'is =')
Q.disc()


# The function `primitive_check(Q,O)` checks whether the binary quadratic form *Q* is primitive, i.e. its coefficients are relatively prime in the ring *O*.

# In[85]:


def primitive_check(Q,O):
    f = Q.polynomial()
    return O.ideal(Q.coefficients()) == O.ideal(1)


# For example, below is an example of a primitive form and a nonprimitive form

# In[86]:


def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
OF = F.ring_of_integers();

a1 = 3*t
b1 = 20
c1 = 1+t
Q1 = QuadraticForm(OF,2,[a1,b1,c1]);
f1 = Q1.polynomial()
print('The form ',f1,'is primitive:',primitive_check(Q1,OF))

a2 = 6
b2 = 2
c2 = 1+t
Q2 = QuadraticForm(OF,2,[a2,b2,c2]);
f2 = Q2.polynomial()
print('The form ',f2,'is primitive:',primitive_check(Q2,OF))


# ---
# 
# ## Ideal Class Pairing

# Given a point *P* on *E* defined by
# >*y*^2 = *x*^3 + a_4 *x* +a_6
# 
# and a *Q* on the nonstandard model of the *disc*-twist *E_twist* of *E*, defined by
# >d(y/2)^2 = x^3 + a_4 x + a_6,
# 
# we define a binary quadratic form with the function `point_pairing_coef`. This function takes as input:
# 
# **Input:**
# 
# 1. trf = A totally real field of class number 1,
# 2. disc = a negative fundamental discriminant,
# 3. P_coord = a list [A,B,C] for which *P* = (*A*/*C*^2,*B*/*C*^3) in *E*(*F*) is in standard form,
# 4. Q_coord = a list [u,v,w] for which *Q* = (*u*/*w*^2,*v*/*w*^3) in *E_twist*(*F*) is in standard form,
# 5. l = an algebraic integer in *F*
# 
# **Output:** A list [a,b,c] which form the coefficients of a binary quadratic form *a x*^2 + *b x y* + *c y*^2.

# In[87]:


def point_pairing_coef(trf,disc,P_coord,Q_coord,l):
    
    OF = trf.maximal_order() # ring of integers of trf
    
    A = P_coord[0]
    B = P_coord[1]
    C = P_coord[2]
    u = Q_coord[0]
    v = Q_coord[1]
    w = Q_coord[2]
    
    alpha = abs(w^2 * A - u * C^2)
    G = OF(alpha).gcd(OF(v^2 * C^6))
    
    coef_xx = alpha/G
    coef_xy = (2 * w^3 * B + l * coef_xx)/(v * C^3)
    coef_yy = ((2 * w^3 * B + l * coef_xx)^2 - disc * v^2 * C^6)/(4 * v^2 * C^6 * coef_xx)
    
    return [coef_xx,coef_xy,coef_yy]


# For example, over *F* = **Q**, the rationals, we have

# In[88]:


#1. Define the elliptic curve E
a4 = -3
a6 = 3
E = EllipticCurve(QQ,[a4,a6]);
print('- E = ',E)

#2.  Twist by disc
disc = -4
E_twist = twist(QQ,a4,a6,disc)
print('- Twist by:',disc)
print('- E_twist =',E_twist) 

#3. Find nontorsion point on E and write it in the form (a/c^2,b/c^3)
P = E.gens()[0]+E.gens()[0]+E.gens()[0]
P_xy = list(P.xy())
P_coord = standard_form(QQ,P_xy)
print('- Nontorsion point E: P =',P,'with standard form:',P_coord)

#4. Find any point on the twist Ed and use the change of variables to take the point to the noncanonical model
Q_twist = E_twist.gens()[0]
print('- Point on E_twist: Q_twist =',Q_twist)

Q = nc_twist_coord(Q_twist.xy(),disc)
Q_coord = standard_form(QQ,Q)
print('  The corresponding point on the noncanonical model is Q =',vector(Q),'with standard form:',Q_coord)

#6. Determine values of the parameter l
box_size = 1000
l_values = list((-box_size..box_size))
#l_values = [coeff[0]+coeff[1]*t for coeff in Tuples((-box_size..box_size),2).list()]

#7. Compute the quadratic forms for those values of l
print('- The quadratic forms with',-box_size-1,'< l <',box_size+1,'with integral coefficients are listed below as [l,[a,b,c]]:')
[[l,point_pairing_coef(QQ,disc,P_coord,Q_coord,l)] for l in l_values if point_pairing_coef(QQ,disc,P_coord,Q_coord,l)[2].is_integral()]


# Another example (not over **Q**)

# In[89]:


#1. Define the totally real field and its ring of integers
def_poly = x^2 - 3
F.<t> = NumberField(def_poly)
OF = F.maximal_order()
print('- F =',F)

#1. Define the elliptic curve E
a4 = -2
a6 = 5
E = EllipticCurve(F,[a4,a6]);
print('- E = ',E)

#2.  Twist by disc
disc = -3
E_twist = twist(F,a4,a6,disc)
print('- Twist by:',disc)
print('- E_twist =',E_twist) 

#3. Find nontorsion point on E and write it in the form (a/c^2,b/c^3)
P = E.gens()[0]
P_xy = list(P.xy())
P_coord = standard_form(F,P_xy)
print('- Nontorsion point E: P =',P,'with standard form:',P_coord)

#4. Find any point on the twist Ed and use the change of variables to take the point to the noncanonical model
Q_twist = E_twist.gens()[0]
print('- Point on E_twist: Q_twist =',Q_twist)

Q = nc_twist_coord(Q_twist.xy(),disc)
Q_coord = standard_form(F,Q)
print('  The corresponding point on the noncanonical model is Q =',vector(Q),'with standard form:',Q_coord)

#6. Determine values of the parameter l
box_size = 2
int_basis = vector(F,OF.gens())
basis_length = len(int_basis)
l_values = [vector(F,coef).dot_product(int_basis) for coef in Tuples((-box_size..box_size),basis_length).list()]

#7. Compute the quadratic forms for those values of l
print('- The quadratic forms with parameters l with integral coefficients bounded by',box_size,'are listed below as [l,[a,b,c]]:')
[[l,point_pairing_coef(F,disc,P_coord,Q_coord,l)] for l in l_values]
#[[l,point_pairing_coef(F,disc,P_coord,Q_coord,l)] for l in l_values if point_pairing_coef(F,disc,P_coord,Q_coord,l)[2].is_integral()]


# We consider the following elliptic curve and its twist

# In[ ]:





# ---
# 
# ## Scratch

# In[90]:


# Ideal Class Pairing over CM fields

def ICP(n,m,d,P):
    # every rational point can be rewritten as (A/C^2,B/C^3)
    A = P.xy()[0].numerator()
    B = P.xy()[1].numerator() 
    C = (P.xy()[1].denominator())/(P.xy()[0].denominator())
    # the two rational points on the twist can be rewritten as (u/w^2,v/w^3)
    u = [0,n]
    v = [2*m,2*m]
    w = [1,1]
    # Define some constants attached to the above coordinates
    alpha = [abs(A*w[0]^2-u[0]*C^2),abs(A*w[1]^2-u[1]*C^2)]
    G = [gcd(alpha[0],C^6 * v[0]^2),gcd(alpha[1],C^6 * v[1]^2)]
    H = [gcd(2*w[0]^3 * B, v[0]*C^3),gcd(2*w[1]^3 * B, v[1]*C^3)]
    var('x')
    k = [Integer(solve_mod((alpha[0]/G[0])*x==-(2*w[0]^3 * B/H[0])-(C^3 * v[0]*d/H[0]),2*C^3 * v[0]/H[0])[0][0]),
        Integer(solve_mod((alpha[1]/G[1])*x==-(2*w[1]^3 * B/H[1])-(C^3 * v[1]*d/H[1]),2*C^3 * v[1]/H[1])[0][0])]
    l = [Integer(mod(H[0]*k[0],2*C^3 * v[0])),Integer(mod(H[1]*k[1],2*C^3 * v[1]))]
    # Coefficients of binary quadratic form
    a = [alpha[0]/G[0],alpha[1]/G[1]]
    b = [(2 * w[0]^3 * B + l[0] * (alpha[0]/G[0]) )/(v[0] * C^3),
         (2 * w[1]^3 * B + l[1] * (alpha[1]/G[1]) )/(v[1] * C^3)]
    c = [((2 * w[0]^3 * B + l[0] * (alpha[0]/G[0]))^2+C^6 * v[0]^2 * d)*G[0]/(4*v[0]^2*C^6*alpha[0]),
        ((2 * w[1]^3 * B + l[1] * (alpha[1]/G[1]))^2+C^6 * v[1]^2 * d)*G[1]/(4*v[1]^2*C^6*alpha[1])]
    return [[a[0],b[0],c[0]],[a[1],b[1],c[1]]]

def ICPQF(n,m,d,P):
    icplist = ICP(n,m,d,P)
    Q = [BinaryQF(icplist[0]),BinaryQF(icplist[1])]
    return [Q[0].reduced_form(),Q[1].reduced_form()]#[0] if Q[0].is_equivalent(Q[1]) else [Q[0].reduced_form(),Q[1].reduced_form()]

# Given a parameters n,m,d and a point P on the curve E(n,m,d), checks w

