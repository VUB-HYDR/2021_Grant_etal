# IPython log file

runfile('/home/luke/documents/python/isimip/final_data_code_env/part3_OF/main_p3.py', wdir='/home/luke/documents/python/isimip/final_data_code_env/part3_OF')
y = np.matrix(y).T
X = np.matrix(X).T
Z = np.transpose(np.matrix(ctl))
nb_runs_x = np.matrix(nb_runs_x)
import numpy as np
y = np.matrix(y).T
X = np.matrix(X).T
Z = np.transpose(np.matrix(ctl))
nb_runs_x = np.matrix(nb_runs_x)
nbts = y.shape[0]
# Spatial dimension
n_spa = (trunc+1)**2
# Spatio_temporal dimension (ie dimension of y)
n_st = n_spa * nbts
# number of different forcings
I = X.shape[1]
nb_runs_ctl = np.shape(Z)[1]
half_1_end = int(np.floor(nb_runs_ctl / 2))
Z1 = Z[:,:half_1_end]
if nb_runs_ctl % 2 == 0:
    Z2 = Z[:,half_1_end]
else:
    Z2 = Z[:,half_1_end+1:]
n_red = n_st - n_spa
U = projfullrank(nbts, n_spa)
runfile('/home/luke/documents/python/isimip/final_data_code_env/part3_OF/proc/PyDnA.py', wdir='/home/luke/documents/python/isimip/final_data_code_env/part3_OF/proc')
n_red = n_st - n_spa
U = projfullrank(nbts, n_spa)
yc = np.dot(U, y)
Z1c = np.dot(U, Z1)
Z2c = np.dot(U, Z2)
Xc = np.dot(U, X)
proj = np.identity(X.shape[1])
Cf = regC(Z1c.T)
Cf1 = np.real(spla.inv(Cf))
#Matrix is singular and may not have a square root. can be ignored
Cf12 = np.real(spla.inv(spla.sqrtm(Cf)))
#Matrix is singular and may not have a square root. can be ignored
Cfp12 = np.real(spla.sqrtm(Cf))
c0, c1, c2, d_cons, x_tilde_white, y_tilde_white = tls(np.dot(Xc.T, Cf12), np.dot(yc.T, Cf12), np.dot(Z2c.T, Cf12), nb_runs_x, proj, formule_ic_tls)
x_tilde = np.dot(Cfp12, x_tilde_white.T)
y_tilde = np.dot(Cfp12, y_tilde_white.T)
beta_hat = c0.T
beta_hat_inf = c1.T
beta_hat_sup = c2.T
X_test = np.dot(Xc.T,Cf12)
Y_test = np.dot(yc.T, Cf12)
Y_test = np.dot(yc.T, Cf12)
Z2_test = np.dot(Z2c.T, Cf12)
nX_test = nb_runs_x
PROJ_test = proj
n = Y.shape[1]
m = X.shape[0]
n = Y_test.shape[1]
m = X_test.shape[0]
n
#[Out]# 32
m
#[Out]# 1
X = np.multiply(X_test, (np.sqrt(nX_test).T * np.ones((1, n))))
if X.shape[0] == 1: # adjusted
    DnX = np.sqrt(nX).squeeze()
else:
    DnX = np.diag(np.sqrt(nX).A1)
X_test = np.multiply(X_test, (np.sqrt(nX_test).T * np.ones((1, n))))
if X_test.shape[0] == 1: # adjusted
    DnX = np.sqrt(nX_test).squeeze()
else:
    DnX = np.diag(np.sqrt(nX_Test).A1)
M = np.vstack([X_test, Y_test])
U, D, V = np.linalg.svd(M)
V = V.T
Uk = U[:, -1]
Uk_proj = np.vstack([PROJ * DnX * Uk[:-1], Uk[-1]])
Uk_proj = np.vstack([PROJ_test * DnX * Uk[:-1], Uk[-1]])
beta_hat = - Uk_proj[:-1] / Uk_proj[-1]
Uk_proj[:-1]
#[Out]# matrix([[-0.8280026]])
Uk_proj[-1]
#[Out]# matrix([[0.9783408]])
beta_hat
#[Out]# matrix([[0.8463335]])
beta_hat_inf = np.zeros(beta_hat.shape)
beta_hat_sup = np.zeros(beta_hat.shape)
D_tilde = np.matrix(np.zeros(M.shape))
np.fill_diagonal(D_tilde, D)
D_tilde[m, m] = 0
Z_tilde = U * D_tilde * V.T
X_tilde = Z_tilde[0:m, :] / (np.dot(np.sqrt(nX).T, np.ones((1, n))))
Y_tilde = Z_tilde[m, :]
D_tilde = np.matrix(np.zeros(M.shape))
np.fill_diagonal(D_tilde, D)
D_tilde[m, m] = 0
Z_tilde = U * D_tilde * V.T
X_tilde = Z_tilde[0:m, :] / (np.dot(np.sqrt(nX_test).T, np.ones((1, n))))
Y_tilde = Z_tilde[m, :]
# The square singular values (denoted by lambda in AS03)
d = D**2
# Computation of corrected singular value (cf Eq 34 in AS03)
d_hat = np.zeros(d.shape)
NAZv = Z2_test.shape[0]
for i in range(len(d)):
    vi = V[:, i].T
    if Formule_IC_TLS == "AS03":
        # Formule Allen & Stott (2003)
        d_hat[i] = d[i] / np.dot(np.dot(np.dot(vi, Z2_test.T), Z2_test / NAZv), vi.T)
    elif Formule_IC_TLS == "ODP": 
         # Formule ODP (Allen, Stone, etc)
        d_hat[i] = d[i] / np.dot(np.power(vi, 2), np.sum(np.power(Z2_test, 2), axis=0).T / NAZv)
Formule_IC_TLS = formule_ic_tls
# The square singular values (denoted by lambda in AS03)
d = D**2
# Computation of corrected singular value (cf Eq 34 in AS03)
d_hat = np.zeros(d.shape)
NAZv = Z2_test.shape[0]
for i in range(len(d)):
    vi = V[:, i].T
    if Formule_IC_TLS == "AS03":
        # Formule Allen & Stott (2003)
        d_hat[i] = d[i] / np.dot(np.dot(np.dot(vi, Z2_test.T), Z2_test / NAZv), vi.T)
    elif Formule_IC_TLS == "ODP": 
         # Formule ODP (Allen, Stone, etc)
        d_hat[i] = d[i] / np.dot(np.power(vi, 2), np.sum(np.power(Z2_test, 2), axis=0).T / NAZv)
d_cons = d_hat[-1]
d_cons
#[Out]# 38.559428751531804
d_H0_cons = consist_mc_tls(Cf, Xc, nb_runs_x, NZ1, NZ2, N_cons_mc, formule_ic_tls)
NZ1 = Z1c.shape[1]
NZ2 = Z2c.shape[1]
N_cons_mc = 1000
d_H0_cons = consist_mc_tls(Cf, Xc, nb_runs_x, NZ1, NZ2, N_cons_mc, formule_ic_tls)
N = len(d_H0)
d_H0 = d_H0_cons
d = d_conds
d = d_cons
N = len(d_H0)
h = 1.06 * np.std(d_H0, ddof=1) * N ** (-1./5)
onem = sps.norm.cdf(d, d_H0, h)
h
#[Out]# 2.4557668248965667
d
#[Out]# 38.559428751531804
onem_test = sps.norm.cdf(d, loc=d_H0[1], scale=h)
onem_test = sps.norm.cdf(d, loc=d_H0[2], scale=h)
