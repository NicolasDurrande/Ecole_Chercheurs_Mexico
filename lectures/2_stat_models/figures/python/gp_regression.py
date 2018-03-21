import numpy as np
import gpflow

import matplotlib.pyplot as plt
from  matplotlib import cm
#import matplotlib2tikz as m2t
plt.close('all')
plt.ion()


colors500 = [ cm.coolwarm(x) for x in np.linspace(0, 1, 500) ]
x = np.linspace(-1, 1, 100)[:, np.newaxis]

##############################
## GP samples
num_samples = 200

k = gpflow.kernels.RBF(1, lengthscales=0.3)
K = k.compute_K(x, x)
Z = np.random.multivariate_normal(np.zeros(x.shape[0]), K, num_samples+10000).T

colors = [ cm.coolwarm_r(x) for x in np.linspace(0, 1, num_samples) ]
#colors = [ cm.bone_r(x) for x in np.linspace(0, 1, num_samples) ]

plt.figure(figsize=(9, 5))
for i in range(num_samples):
    plt.plot(x, Z[:,i], color=colors[i], linewidth=.8, alpha=.8)
plt.xlim((-1, 1))
plt.ylim((-4, 4))
plt.xticks([-1,1], fontsize=14)
plt.yticks([-4, 0, 4], fontsize=14)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$Y(x)$', fontsize=18)

plt.tight_layout()
plt.savefig('gp_samples.pdf', transparent=True)

##############################
## GP definition
n = 10
idx = np.random.choice(np.arange(x.shape[0]), n)

for id in idx:
    plt.axvline(x[id], 0, 1, color='k', linestyle='--', linewidth=.5)

plt.savefig('gp_definition_a.pdf', transparent=True)


A = np.random.uniform(-2, 2, (1, n))

Y = A @ Z[idx, :]
Y = Y.flatten()

plt.figure(figsize=(5,5))
plt.hist(Y, bins=21, density=True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('$\\alpha^T Y(X)$', fontsize=18)
plt.ylabel('', fontsize=18)
plt.tight_layout()
plt.savefig('gp_definition_b.pdf', transparent=True)

##############################
## GP distribution
x = np.linspace(-1, 1, 200)[:, np.newaxis]

num_samples = 20
colors = [ cm.coolwarm_r(x) for x in np.linspace(0, 1, num_samples) ]

means = [3*x, 0*x]
kerns = [gpflow.kernels.RBF, gpflow.kernels.Matern12]

for j in range(2):
    mean = means[j]
    ki = kerns[j](1, lengthscales=0.3)
    
    kx = ki.compute_K(x, np.zeros((1,1)))

    K = ki.compute_K(x, x)
    Z = mean + np.random.multivariate_normal(np.zeros(x.shape[0]), K, num_samples).T

    z_range = (np.min(Z), np.max(Z))
    plt.figure(figsize=(11, 3.5))

    plt.subplot(1, 3, 1)
    plt.plot(x, mean)
    plt.xticks([-1, 1], fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('$m(x)$', fontsize=18)
    plt.ylim(z_range)

    plt.subplot(1, 3, 2)
    plt.plot(x, kx)
    plt.xticks([-1, 0, .5, 1], fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('$k(0, x)$', fontsize=18)
    plt.axvline(0, 0, 1, color='k', linestyle='--', linewidth=.5)
    plt.axvline(0.5, 0, 1, color='k', linestyle='--', linewidth=.5)

    plt.subplot(1, 3, 3)
    for i in range(num_samples):
        plt.plot(x, Z[:,i], color=colors[i], linewidth=.8, alpha=.8)
    plt.xlim((-1, 1))
    plt.ylim((-4, 4))
    plt.xticks([-1, 0, .5, 1], fontsize=14)
    plt.yticks(fontsize=14)
    plt.title('$Y(x)$', fontsize=18)
    plt.ylim(z_range)
    plt.axvline(0, 0, 1, color='k', linestyle='--', linewidth=.5)
    plt.axvline(0.5, 0, 1, color='k', linestyle='--', linewidth=.5)

    plt.tight_layout()
    plt.savefig('gp_distribution_{}.pdf'.format(j), transparent=True)


##############################
## observations
def f(x):
    y = 2*x + 2*np.sin(2*np.pi*x)
    return y

n = 7
X = np.linspace(-.7, .7, n)[:, np.newaxis]
Y = f(X)

plt.figure(figsize=(9, 5))
plt.scatter(X, Y, marker='x', color='k', s=25, linewidth=2)
plt.xlim((-1, 1))
plt.ylim((-5, 5))
plt.xticks([-1,1], fontsize=14)
plt.yticks([-5, 0, 4], fontsize=14)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$f$', fontsize=18)

plt.tight_layout()
plt.savefig('gp_obs.pdf', transparent=True)

##############################
## GP prior samples
num_samples = 200

k = gpflow.kernels.RBF(1, variance=2., lengthscales=0.15)
K = k.compute_K(x, x)
Z = np.random.multivariate_normal(np.zeros(x.shape[0]), K, num_samples+10000).T

colors = [ cm.coolwarm_r(x) for x in np.linspace(0, 1, num_samples) ]
#colors = [ cm.bone_r(x) for x in np.linspace(0, 1, num_samples) ]

plt.figure(figsize=(9, 5))
for i in range(num_samples):
    plt.plot(x, Z[:,i], color=colors[i], linewidth=.8, alpha=.8)
plt.xlim((-1, 1))
plt.ylim((-5, 5))
plt.xticks([-1,1], fontsize=14)
plt.yticks([-4, 0, 4], fontsize=14)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$Y(x)$', fontsize=18)

plt.tight_layout()
plt.savefig('gp_samples2.pdf', transparent=True)


#############################
## conditional samples

m = gpflow.models.GPR(X, Y, kern=k)
m.likelihood.variance = 0.001

num_samples = 50

Z = m.predict_f_samples(x, num_samples)[:,:,0].T

colors = [ cm.coolwarm_r(x) for x in np.linspace(0, 1, num_samples) ]
#colors = [ cm.bone_r(x) for x in np.linspace(0, 1, num_samples) ]

plt.figure(figsize=(9, 5))
for i in range(num_samples):
    plt.plot(x, Z[:,i], color=colors[i], linewidth=.8, alpha=.8)
plt.xlim((-1, 1))
plt.ylim((-5, 5))
plt.xticks([-1,1], fontsize=14)
plt.yticks([-4, 0, 4], fontsize=14)
plt.xlabel('$x$', fontsize=18)
plt.ylabel('$Y(x) | Y(X)=F$', fontsize=18)
plt.plot(X, Y, 'kx', mew=2)

plt.tight_layout()
plt.savefig('gp_cond_samples.pdf', transparent=True)

def plot(x, m):
    mean, var = m.predict_f(x)
    plt.fill_between(x[:,0],
                     mean[:,0] - 2*np.sqrt(var[:,0]),
                     mean[:,0] + 2*np.sqrt(var[:,0]),
                     color='C0', alpha=0.2)
    plt.plot(x, mean, 'C0', lw=2)
    plt.plot(x, mean[:,0]-2*np.sqrt(var[:,0]), 'C0', lw=.5)
    plt.plot(x, mean[:,0]+2*np.sqrt(var[:,0]), 'C0', lw=.5)
    plt.xlim((-1, 1))
    plt.ylim((-5, 5))
    plt.xticks([-1,1], fontsize=14)
    plt.yticks([-4, 0, 4], fontsize=14)
    plt.xlabel('$x$', fontsize=18)
    plt.ylabel('$Y(x) | Y(X)=F$', fontsize=18)
    plt.plot(X, Y, 'kx', mew=2)

plt.figure(figsize=(9, 5))
plot(x, m)
plt.tight_layout()
plt.savefig('gp_model.pdf', transparent=True)

