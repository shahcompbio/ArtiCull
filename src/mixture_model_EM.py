import functools
import numpy as np
from random import choices
from numpy.random import beta, poisson, binomial
from scipy.stats import beta as beta_rv, binom as binom_rv, gaussian_kde
import scipy.integrate as integrate
import scipy.special

from utils_io import update_progress

eps = 0.005
prog = 0

def initialize_EM(df):
    assignments = df.apply(initial_assignment, axis=1)
    df['assignment'] = assignments
    weights = [len([a for a in assignments if a == i])/len(df) for i in range(4)]
    dists = estimate_kde_hard(df)
    return df, weights, dists

def initial_assignment(mut):
    #var, tot, vaf, cn = mut
    # 0 = CLONAL
    # 1,2 = Subclonal Private
    # 3 = Shared subclonal

    if mut['ccf_1'] > 0.6 and mut['ccf_2'] > 0.6: return 0
    elif mut['ccf_2'] <.05 : return 1
    elif mut['ccf_1'] <.05 : return 2
    else: return 3

def EM(dists, weights, df):
    print("Computing Responsibilities")

    df_subset = df.sample(n=min(len(df), 1000), replace=False)
    R = df_subset.parallel_apply(lambda x: compute_responsibilities(x, dists, weights), axis=1).tolist()
    print("Updating Distributions")
    dists = estimate_kde_soft(dists, df_subset, R)
    weights = [sum([r[i] for r in R])/len(df_subset) for i in range(4)]

    return dists, weights, R

def compute_responsibilities(mut, dists, weights):
    #global prog
    #if prog % 200 == 0: print(prog)
    #prog += 1
    Z_1, Z_2, Z_a1 = dists

    v1, t1 = mut['var_1'], mut['tot_1']
    v2, t2 = mut['var_2'], mut['tot_2']
    c1, c2 = mut['cn_1'], mut['cn_2']

    if c1 == 0 or c2 == 0:
        return [0,0.5,0.5,0]

    r = []
    # Public cluster
    ## This does not properly account for error
    a = binom_rv(t1, 1/c1).pmf(v1)
    b = binom_rv(t2, 1/c2).pmf(v2)
    r.append(weights[0]*a*b)

    # Private Clone 1
    a = compute_integral(Z_1, v1, t1, c1)
    b = binom_rv(t2, eps).pmf(v2)
    r.append(weights[1]*a*b)

    # Private Clone 2
    a= binom_rv(t1, eps).pmf(v1)
    b = compute_integral(Z_2, v2, t2, c2)
    r.append(weights[2]*a*b)

    # Artifact
    a = compute_integral(Z_a1, v1, t1, c1)
    b = compute_integral(Z_a1, v2, t2, c2)
    r.append(weights[3]*a*b)

    return r

@functools.lru_cache
def binom_coeff(n, x):
    # Memoized to drastically speed up repeated calls
    return scipy.special.binom(n, x)

def binom_pmf(n,p,x):
    return binom_coeff(n,x)*p**x*(1-p)**(n-x)

def compute_integral(Z, v, t, cn):
    def integrand(x):
        x_adj = x/cn - (2*x*eps)/cn + eps
        return Z.pdf(x)[0] * binom_pmf(t, x_adj, v)
    return integrate.quad(integrand, a = 0, b = 1)[0]

def estimate_kde_soft(prev_dists, df, R):

    # Private Clone 1
    points = df['ccf_1']
    weights = [r[1] for r in R]
    Z_1 = gaussian_kde(points, weights = weights)

    # Private Clone 2
    points = df['ccf_1']
    weights = [r[2] for r in R]
    Z_2 = gaussian_kde(points, weights = weights)

    # Artifact
    points_1 = df['ccf_1']
    points_2 = df['ccf_2']
    weights = [r[3] for r in R]
    Z_a1 = gaussian_kde(list(points_1)+list(points_2), weights = weights+weights)

    # TODO: Remove unused distributions
    return Z_1, Z_2, Z_a1

def estimate_kde_hard(df):

    # Private Clone 1
    points = df[df['assignment'] == 1]['ccf_1']
    Z_1 = gaussian_kde(points)

    # Private Clone 2
    points = df[df['assignment'] == 2]['ccf_2']
    Z_2 = gaussian_kde(points)

    # Artifact
    points_1 = df[df['assignment'] == 3]['ccf_1']
    points_2 = df[df['assignment'] == 3]['ccf_2']
    Z_a1 = gaussian_kde(list(points_1) + list(points_2))

    # TODO: Remove unused distributions
    return Z_1, Z_2, Z_a1
