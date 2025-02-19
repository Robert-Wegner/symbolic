
from sympy import *
import time
import pickle
import copy

epsilon, epsiloninv, u, v, a, theta = symbols('epsilon, eta, u, v, a, \\theta', real=True)
q = symbols('q', real=False)

differentiable_symbols = [u, v, a, theta, q]

def var_deriv_name(var):
    if "_x" in var.name:
        return var.name + 'x'
    else:
        #return '(' + var.name + ')_x'
        return var.name + '_x'

    
def deriv(poly):
    res = 0
    original = differentiable_symbols.copy()
    for sym in original:
        deriv_term = Derivative(poly, sym).doit()
        if deriv_term != 0:
            #print("looking at: ", sym)
            newName = var_deriv_name(sym)

            if newName in [s.name for s in differentiable_symbols]:
                dsym = differentiable_symbols[[s.name for s in differentiable_symbols].index(newName)]
                #print("newName found: ", dsym)
            else:
                dsym = Symbol(newName, real=True)
                differentiable_symbols.append(dsym)
                #print("newName not found: ", dsym)
            #print("syms: ", syms)
            res += deriv_term * dsym
            #print("res: ", res)
    
    return res

def higher_deriv(var, n):
    if n > 0:
        return deriv(higher_deriv(var, n-1))
    else:
        return var
    
def polynomize(expr):
    return simplify(Poly(expr, epsilon, epsiloninv).subs(epsiloninv, 1/epsilon))

def substituter(expr, var, sub, magnitude=1, scale=1):
    original = differentiable_symbols.copy()

    dvar = var
    dsub = sub
    expr = expr.subs(dvar, magnitude * dsub)
    cont = True
    while cont:
        dvar = scale * deriv(dvar)
        dsub = scale * deriv(dsub)
        cont = False
        if Derivative(expr, dvar).doit() != 0:
            expr = expr.subs(dvar, magnitude * dsub)

        if var_deriv_name(dvar) in [sym.name for sym in original]:
            cont = True
            
    return expr
    
def multi_substituter(expr, data, magnitude=1, scale=1):
    for (var, sub) in data:
        expr = substituter(expr, var, sub, scale)
    return expr




#test = pickle.load(open("sigma_seq_order9.p", "rb"))
#display(test[4])


def custom_display(expr, title):
    print("--- Real: ", title, "---")
    display(polynomize(re(expr)))
    #display(polynomize(multi_substituter(re(expr), [(u, theta/2 + a), (v, theta/2 - a)])))
    print("--- Imag: ", title, "---")
    display(polynomize(im(expr)))
    #display(polynomize(multi_substituter(im(expr), [(u, theta/2 + a), (v, theta/2 - a)])))
    
r_factor = 1 #or -1?
p_factor = 1

def p_coeff(n):
    if n == -1:
        return simplify(epsilon / sqrt(2))
    if n >= 0:
        return p_factor * simplify(1 / sqrt(2) * epsiloninv**(2*n+1) * Rational((-2)**(-n) * (factorial(2*n) / (factorial(n) * factorial(n+1)))))
    return 0

def p_coeff_inv(n):
    if n == -1:
        return simplify(sqrt(2) * epsiloninv)
    if n >= 0:
        return p_factor * simplify(sqrt(2) * epsilon**(2*n+1) * Rational((-2)**n * (factorial(n) * factorial(n+1)) / (factorial(2*n))))
    return 0

r_seq = []

sigma_seq = []

N = 20

F_list = []
factors = []
combo = []
combo_tracker = []
X = 0
Y = 0
nn_is_even = True

t = t = time.time()
for n in range(-1, N):
    print(" --- nn = ", n+1, " ---")

    if n == -1:
        r_seq.append(r_factor * I * epsiloninv * sqrt(2))
        #custom_display(r_seq[-1], "r_0")

        sigma_seq.append((u + v) * epsilon**2 * p_coeff_inv(-1) / (4 * I))
        #custom_display(sigma_seq[-1], "sigma_0")
    
    elif n == 0:
        r_seq.append(r_factor * I / sqrt(2) * epsiloninv * (u - v - 2 * epsiloninv**2))
        #custom_display(r_seq[-1], "r_1")

        sigma_seq.append(simplify(deriv(sigma_seq[0]) + sigma_seq[0]**2 - 2 * I * p_coeff(0) * sigma_seq[0] + epsilon**2 / 2 * deriv(v) * r_seq[0] 
                                   + u * (1 + epsilon**2 / 2 * v)))
        #custom_display(sigma_seq[-1], "sigma_1")
    else:
        r_square_term_1 = 0
        for k in range(0, n+1):
            r_square_term_1 += r_seq[k] * r_seq[n - k]

        r_square_term_2 = 0
        for k in range(1, n+1):
            r_square_term_2 += r_seq[k] * r_seq[n + 1 - k]

        sigma_square_term = 0
        for k in range(0, n+1):
            sigma_square_term += (sigma_seq[k] - 2 * I * p_coeff(k) * (2 * I * p_coeff(-1))**k) * sigma_seq[n - k]

        r_seq.append(r_factor * simplify(I * epsiloninv / sqrt(2) 
                              * (- deriv(r_seq[n]) + (1 + epsilon**2 / 2 * v) * r_square_term_1 + epsilon**2 / 2 * r_square_term_2)))
        #custom_display(r_seq[-1], "r_" + str(n+1))

        sigma_seq.append(simplify(deriv(sigma_seq[n]) + sigma_square_term 
                                   + epsilon**2 / 2 * (2 * I * p_coeff(-1))**n * deriv(v) * r_seq[n]))
        #custom_display(sigma_seq[-1], "sigma_" + str(n+1))

    
    nn = n + 1
    factors.append([])
    F_list.append(Symbol('F_{' + str(nn) + '}'))
    if nn_is_even:
        X = im(sigma_seq[nn])
        Y = F_list[nn]
    else:
        X = re(sigma_seq[nn])
        Y = F_list[nn]
    #display(polynomize(X))
    k_is_even = False
    uu = v
    vv = u
    for k in range(1, nn+1):
        print("nn = ", nn, "k = ", k)
        if k_is_even:
            kk = k // 2 - 1
            var = higher_deriv(uu, 2*kk)
            
            prev_coeff = Poly(polynomize(combo[k-1]).subs(epsilon, 0), var).coeff_monomial(var)
            target_coeff = Poly(polynomize(epsilon**(nn-1 - 2*kk) * X).subs(epsilon, 0), var).coeff_monomial(var)
            
            factor = - target_coeff / prev_coeff * epsiloninv**(nn-1 - 2*kk)
            X += factor * combo[k-1]
            Y += factor * combo_tracker[k-1]

            #display("add to X combo at ", k-1, "with factor", polynomize(- target_coeff / prev_coeff * epsiloninv**(n-1 - 2*kk)))
            factors[-1].append(polynomize(factor))
            #display("now X = : ", polynomize(X))
            
        else:
            kk = k // 2
            var = higher_deriv(vv, 2*kk)
            
            prev_coeff = Poly(polynomize(combo[k-1] * epsiloninv).subs(epsilon, 0), var).coeff_monomial(var)
            target_coeff = Poly(polynomize(epsilon**(nn-1 - 2*kk) * X).subs(epsilon, 0), var).coeff_monomial(var)
            
            factor = - target_coeff / prev_coeff * epsiloninv**(nn - 2*kk)
            X += factor * combo[k-1]
            Y += factor * combo_tracker[k-1]

            #display("add to from X combo at ", k-1, "with factor", polynomize(- target_coeff / prev_coeff * epsiloninv**(n - 2*kk)))
            factors[-1].append(polynomize(factor))
            #display("now X = : ", polynomize(X))

        k_is_even = not k_is_even

    #display("combo: ", polynomize(X))   

    combo.append(copy.deepcopy(X))
    combo_tracker.append(copy.deepcopy(Y))
    
    output_folder = "data2/"
    pickle.dump(sigma_seq, open(output_folder + "sigma_seq_" + str(nn) + ".p", "wb"))
    pickle.dump(r_seq, open(output_folder + "r_seq_" + str(nn) + ".p", "wb"))
    pickle.dump(combo, open(output_folder + "combo_" + str(nn) + ".p", "wb"))
    pickle.dump(combo_tracker, open(output_folder + "combo_tracker_" + str(nn) + ".p", "wb"))
    pickle.dump(factors, open(output_folder + "factors_" + str(nn) + ".p", "wb"))
    
    nn_is_even = not nn_is_even
    
    print(time.time() - t)
    t = time.time()

