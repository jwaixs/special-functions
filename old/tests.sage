#
#f = [ base_U(replace_z_x(weights[i][i-1,0]))[-1].full_simplify() for i in range(3, 10) ]
#g = [ base_U(replace_z_x(weights[i][i-2,1]))[-1].full_simplify() for i in range(3, 10) ]
#h = [ base_U(replace_z_x(weights[i][i-3,2]))[-1].full_simplify() for i in range(4, 10) ]
#k = [ base_U(replace_z_x(weights[i][i-4,3]))[-1].full_simplify() for i in range(5, 10) ]
#k2 = [ base_U(replace_z_x(weights[i][i-5,4]))[-1].full_simplify() for i in range(6, 10) ]
#
#for i in range(len(f)): print (f[i] / (1 - q**2) * (1 - q**(2*(i + 3))) ).full_simplify()
#for i in range(len(g)): print (g[i] / qbinomial(i + 3, 1, q**2)**(-2) / (1 - q**4) * (1 - q**(2*(i + 3) - 2)) ).full_simplify()
#for i in range(len(h)): print (h[i] / qbinomial(i + 4, 2, q**2)**(-2) / ((1 - q**(2*(i+4) + 2))/(1 - q**2)) / ((1 - q**6)/(1 - q**(2*(i+4)-4))) ).full_simplify()
#for i in range(len(k)): print (k[i] / qbinomial(i + 5, 3, q**2)**(-2) / ((1 - q**(2*(i+5) + 2))/(1 - q**2)) / ((1 - q**8)/(1 - q**(2*(i+5)-6))) ).full_simplify()
#for i in range(len(k2)): print (k2[i] / qbinomial(i + 5, 4, q**2)**(-2) / ((1 - q**10)/(1 - q**(2*(i+5) - 8))) / ((1 - q**(2*(i + 5) + 2))/(1 - q**2)) ).full_simplify() 

