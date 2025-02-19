from sympy import *
from IPython.display import display
import time

epsilon, epsiloninv, u, v, a, theta = symbols('epsilon, eta, u, v, a, \\theta', real=True)
q = symbols('q', real=False)

differentiable_symbols = [u, v, a, theta, q]

R = [Symbol('r^{' + str(n) + '}') for n in range(13)]

for r in R:
    differentiable_symbols.append(r)


def var_deriv_name(var):
    if "_{x" in var.name:
        return var.name[:-1] + 'x}'
    else:
        return '(' + var.name + ')_{x}'
        #return var.name + '_x'

    
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

def depolynomize(poly):
    monoms = poly.monoms()
    coeffs = poly.coeffs()
    
    X = 0
    for (k, m) in enumerate(monoms):
        X += epsilon**m[0] * epsiloninv**m[1] * coeffs[k]
            
    return X

def poly_simplify(expr):
    return depolynomize(polynomize(expr))

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


import pickle
import datetime

output_folder = "res/"

def p_coeff(n):
    if n == -1:
        return simplify(epsilon / sqrt(2))
    if n >= 0:
        return simplify(1 / sqrt(2) * epsiloninv**(2*n+1) * Rational((-2)**(-n) * (factorial(2*n) / (factorial(n) * factorial(n+1)))))
    return 0

def q_coeff(n):
    if n % 2 == 0 or n < 0:
        return 0
    else:
        nn = (n - 1) // 2
        if nn % 2 == 0:
            return I**n * epsiloninv**n * catalan(nn) / (sqrt(2) * 2**nn)
        elif nn % 2 == 1:
            return - I**n * epsiloninv**n* catalan(nn) / (sqrt(2) * 2**nn)

def a_coeff(n):
    return catalan(n) / 2**n
    
N = int(input("Max iterations? "))

loading = input("Load from previous results? (Y/N) ")

r_sym = [Symbol('r_{' + str(n) + '}') for n in range(N+1)]
r_eq = [r_sym[n] for n in range(N+1)]
r_poly = [r_sym[n] for n in range(N+1)]

sig_sym = [Symbol('\sigma_{' + str(n) + '}') for n in range(N+1)]
sig_eq = [sig_sym[n] for n in range(N+1)]
sig_poly = [sig_sym[n] for n in range(N+1)]

E_sym = [Symbol('E_{' + str(n) + '}') for n in range(N+1)]
E_eq = [E_sym[n] for n in range(N+1)]
E_poly = [E_sym[n] for n in range(N+1)]

Et_sym = [Symbol('\\tilde{E}_{' + str(n) + '}') for n in range(N+1)]
Et_eq = [Et_sym[n] for n in range(N+1)]
Et_poly = [Et_sym[n] for n in range(N+1)]

variant_sign = -1

r_factor = -1 #Symbol('\lambda')
r_eq[0] = r_factor * (- I * sqrt(2) * epsiloninv)
r_eq[1] = r_factor * (u - epsiloninv**2 + (Rational(1, 2) + epsilon**2 * v**2 / 2) * r_eq[0]**2)
r_poly[0] = polynomize(r_eq[0])
r_poly[1] = polynomize(r_eq[1])

sig_eq[0] = epsilon / (2 * sqrt(2) * I) * (u + v)
sig_eq[1] = deriv(sig_eq[0]) - sig_eq[0]**2 + u * (1 + epsilon**2 * v / 2) - (u + v) / 2 + epsilon**2 * deriv(v) / 2 * r_eq[0]
sig_poly[0] = polynomize(sig_eq[0])
sig_poly[1] = polynomize(sig_eq[1])

E_eq[0] = re(sig_eq[0])
Et_eq[0] = im(sig_eq[0])
E_poly[0] = polynomize(E_eq[0])
Et_poly[0] = polynomize(Et_eq[0])

E_eq[1] = re(sig_eq[1]) - variant_sign * epsiloninv * sqrt(2) * Et_eq[0]
Et_eq[1] = im(sig_eq[1])
E_poly[1] = polynomize(E_eq[1])
Et_poly[1] = polynomize(Et_eq[1])

for n in range(1, N):
    
    print("n =", n)
    print(datetime.datetime.now())

    try:
        if loading == "Y":
            print("Loading next iteration.")
            sig_eq = pickle.load(open(output_folder + "sig_eq_" + str(n+1) + ".p", "rb"))
            sig_poly = pickle.load(open(output_folder + "sig_poly_" + str(n+1) + ".p", "rb"))
            r_eq = pickle.load(open(output_folder + "r_eq_" + str(n+1) + ".p", "rb"))
            r_poly = pickle.load(open(output_folder + "r_poly_" + str(n+1) + ".p", "rb"))
            E_eq = pickle.load(open(output_folder + "E_eq_" + str(n+1) + ".p", "rb"))
            E_poly = pickle.load(open(output_folder + "E_poly_" + str(n+1) + ".p", "rb"))
            Et_eq = pickle.load(open(output_folder + "Et_eq_" + str(n+1) + ".p", "rb"))
            Et_poly = pickle.load(open(output_folder + "Et_poly_" + str(n+1) + ".p", "rb"))
            print("Loading successful.")
            print("Skipping this iteration.")

        if loading == "N":
            print("Loading disabled.")
            raise Exception("Loading disabled.")
        
    except:
        print("Loading failed.")
        print("Computing this iteration.")

        X = 0
        
        X += - deriv(r_eq[n]) + q_coeff(n) * epsiloninv**2
        
        for k in range(0, n+1):
            X += (Rational(1, 2) + epsilon**2 * v**2 / 2) * r_eq[k] * r_eq[n-k]
            
        for k in range(1, n+1):
            X += epsilon / (2 * sqrt(2) * I) * r_eq[k] * r_eq[n+1-k]
            
        for k in range(0, n+1):
            for j in range(0, k+1):
                X += q_coeff(j) * r_eq[k-j] * r_eq[n-k] / 2
            
        r_eq[n+1] = r_factor * X
        r_poly[n+1] = polynomize(r_eq[n+1])
        r_eq[n+1] = depolynomize(r_poly[n+1])
        
        X = 0
        
        X += deriv(sig_eq[n]) + q_coeff(n) / 2 * (u + v)
        
        for k in range(0, n+1):
            X += - sig_eq[k] * sig_eq[n-k]
            
        X += epsilon**2 * deriv(v) / 2 * r_eq[n]
        
        sig_eq[n+1] = X
        sig_poly[n+1] = polynomize(sig_eq[n+1])
        sig_eq[n+1] = depolynomize(sig_poly[n+1])

        E_eq[n+1] = re(sig_eq[n+1]) - variant_sign * epsiloninv * sqrt(2) * Et_eq[n]
        E_poly[n+1] = polynomize(E_eq[n+1])
        E_eq[n+1] = depolynomize(E_poly[n+1])

        Et_eq[n+1] = im(sig_eq[n+1])
        for k in range(n-1, -1, -2):
            Et_eq[n+1] += Et_eq[k] * epsiloninv**(n+1-k) * a_coeff((n+1-k)//2 - 1)
        Et_poly[n+1] = polynomize(Et_eq[n+1])
        Et_eq[n+1] = depolynomize(Et_poly[n+1])


        try: 
            print("Saving results.")
            pickle.dump(sig_eq, open(output_folder + "sig_eq_" + str(n+1) + ".p", "wb"))
            pickle.dump(sig_poly, open(output_folder + "sig_poly_" + str(n+1) + ".p", "wb"))
            pickle.dump(r_eq, open(output_folder + "r_eq_" + str(n+1) + ".p", "wb"))
            pickle.dump(r_poly, open(output_folder + "r_poly_" + str(n+1) + ".p", "wb"))
            pickle.dump(E_eq, open(output_folder + "E_eq_" + str(n+1) + ".p", "wb"))
            pickle.dump(E_poly, open(output_folder + "E_poly_" + str(n+1) + ".p", "wb"))
            pickle.dump(Et_eq, open(output_folder + "Et_eq_" + str(n+1) + ".p", "wb"))
            pickle.dump(Et_poly, open(output_folder + "Et_poly_" + str(n+1) + ".p", "wb"))
            print("Saving successful.")
        except:
            print("Saving failed.")
