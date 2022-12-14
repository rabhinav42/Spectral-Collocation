Short description of each matlab program:

1. a_put, eur_call, and b_call are used to get the actual results.
2. call_act and uno_act calculate the actual values of European call and Up and Out Barrier call
   for given parameters. (d12 is a function used for this computation, to get d+ and d- values)
3. put_act calculates the value of an American put with given parameters using a binomial tree
   approach.
4. xi returns the CGL points for the given value of N. cheb calculates the D1 matrix. 
   lagrangeD returns [D1, D2].
5. lagrangeN finds Lj(x) for given parameters.
6. X does the transformation Xi(x) = S
7. R returns R1,i(x) and R2,i(x).
8. du returns dui for given parameters.
9. Z returns a matrix Zi whose entries are Zi,j,k.
10. phi1 and phi2 calculate boundary conditions based on the option and psi1 the initial
   condition. g calculates the penalty term.
11. F_call, despite the name, F(tau, U(tau)) for any type of option.
12. solve_ivp actually solves the IVP dU/dtau = F(tau, U) using variable step size.
13. solve_ivp_const solves the same IVP using constant step size.