The file eqnx.gp contains the function fxj(N,x,j), which computes the modular function F_N(x,j) corresponding to the elliptic curve E_{a1}(N). The inputs are:
- N: conductor of the elliptic curve  
- x and j: variables of the output polynomial F_N(x,j)

The file eqnxJ.gp contains the function fxJ(N,x,J), which computes the modular function f_N(x,J) corresponding to E_{a1}(N). The inputs are:
- N: conductor  
- x and J: variables of the output polynomial f_N(x,J)

The file eqny.gp contains the function fyj(N,y,j), which computes the modular function G_N(y,j) for E_{a1}(N). The inputs are:
- N: conductor  
- y and j: variables of the output polynomial G_N(y,j)

The file eqnyJ.gp contains the function fyJ(N,y,J), which computes the modular function g_N(y,J) for E_{a1}(N). The inputs are:
- N: conductor  
- y and J: variables of the output polynomial g_N(y,J)
Note: All output results are saved in the C:\ directory.

The following commands compute the results in Example 2.3:
gp >read("F://eqnx.gp");
gp >fxj(11,x,j)
gp >read("F://eqnxJ.gp");
gp >fxJ(11,x,j)

The file F389xj.nb contains the modular polynomial F₃₈₉(x,j) in Mathematica format.

The file 'Computational results for Algorithm 1.zip' contains modular polynomials constructed by Algorithm 1 for elliptic curves with conductors not exceeding 100.



The file fibre.gp contains the function fibre(x0,y0,ell,N) for computing fibers of elliptic curves with conductor N at point (x0,y0).As a complete example, this file includes functions fxj(x,j) and fxJ(x,J), which correspond to F₄₆(x,j) and f₄₆(x,J) respectively.To compute fibers for elliptic curves with other conductors, you need to replace these fxj(x,j) and fxJ(x,J) functions with the corresponding modular functions for the desired conductor.
The following code sample demonstrates the computation process. Note that the output only contains three results, excluding the two cusps mapped to (4,-2) - see Example 6.1 for reference.
gp >default(parisize,"200G");
gp > read("F://fibre.gp");
gp > fibre(4,-2,"46a1",46)
%2 = [-0.50000000000000000000000000000000000000 + 0.10425720702853738133907474052527595478*I, 0.48373983739837398373983739837398373984 + 0.00084761956933770228730955073597785329086*I, 0.46774193548387096774193548387096774194 + 0.0033631357105979800431959593717830953156*I]


The file cuspv.gp contains the function cuspv(ell,N0) for computing values at cusps for an elliptic curve "ell" with conductor N0.
The output consists of several vectors. In each vector:
The first component represents the numbering of the cusp;
The second component represents the cusp s/r itself;
The third component gives the value of ϕ at the cusp s/r;
The fourth component contains the integral value 
γ(s/r) computed in the fourth substep 4 of Step 2 in Algorithm 4.
The following code demonstrates the computation process, showing how to calculate values of the modular parametrization ϕ at cusps for the elliptic curve E_{a1}(46) when N=46.
gp > read("F://cuspv.gp");
gp > cuspv("46a1",46)
%4 = [[1, 0, [4.0000000000000000000000000000000000000 + 3.3849542192879981816379771130809918687 E-75*I, -2.0000000000000000000000000000000000003 - 1.6451305057994868586671099363494727687 E-37*I], 0.66090411132397185138974328389639853101 - 5.141032832553433659 E-39*I], [2, 1/2, [0], -1.3218082226479437027794865677927970620 - 1.9834712217832027509 E-38*I], [3, 1/23, [4.0000000000000000000000000000000000000 - 2.0145356427244756343287825513810855766 E-72*I, -2.0000000000000000000000000000000000184 + 1.7544786153318167694980095844993546515 E-36*I], 1.9827123339719155541692298516891955925 + 5.482745673956431878 E-38*I], [4, 1/46, [0], 7.640713280344868802 E-38 - 1.0823589833745110963 E-37*I]].
Note that this function only implements the approximate computation in Step 2 of Algorithm 4,
not the complete algorithm. The full algorithm can be easily obtained.


The file yang.gp contains the function Yang(e,N,u,X,Y,eX,eY,ex,x,y,nx,ny,dx,dy). Below is an explanation of each parameter:

- Parameter u is the q-expansion of the x-coordinate after modular parametrization of elliptic curve E  
- Parameters X and Y are the so-called Yang-Pairs, which are modular functions of the function field of X₀(N) chosen by Yang. They are used to give the defining equation of X₀(N), and here X and Y are given by their q-expansions. The file yang.gp contains Yang-Pairs for X₀(46). To compute modular curves for other conductors, you need to replace X and Y accordingly.  
- eX and eY are the minimal degrees of X and Y respectively  
- ex is the minimal degree of u  
- x and y are variables used in the output polynomials  
- nx and ny are the degrees of x and y in the numerator polynomial of the output rational function  
- dx and dy are the degrees of x and y in the denominator polynomial of the output rational function  
- The output is a vector of two polynomials: the first polynomial is P1, the second is Q1, with x = P1/Q1  

Note: The precision of u, X and Y must be greater than max((nx+1)*(ny+1),(dx+1)*(dy+1)). By default, the output is saved in the C:\ directory.

The following commands compute the results in Example 9.2:
gp > read("F://yang.gp");
gp > e = ellinit("46a1");
gp > [u,v] = elltaniyama(e,300);
The following command computes the two polynomials representing x where x = P1/Q1:
gp > Yang(e,46,u,X,Y,6,7,2,x,y,4,4,4,3)
The following command computes the two polynomials representing y where y = P2/Q2:
gp > Yang(e,46,v,X,Y,6,7,3,x,y,4,5,4,4)
