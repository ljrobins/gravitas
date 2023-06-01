import scipy
import gravitas
import pyspaceaware as ps
import numpy as np
import pyvista as pv

def pairwise_norm(rv1: np.ndarray, rv2: np.ndarray, return_deep_grids: bool = False) -> np.ndarray:
    r1s, r2s = rv1.shape[0], rv2.shape[0]
    rv1_deep = rv1.reshape((r1s, 1, 3))
    rv2_deep = rv2.reshape((r2s, 1, 3))
    rv1_grid = np.swapaxes(np.tile(rv1_deep, (r2s, 1)), 0, 1)
    rv2_grid = np.repeat(rv2_deep, repeats=r1s, axis=1)
    rv_norm = np.linalg.norm(rv1_grid-rv2_grid, axis=2)
    if not return_deep_grids:
        return rv_norm
    else:
        return rv_norm, rv1_grid, rv2_grid

def fun_with_legacy(x: np.ndarray, 
                    acc_truth: np.ndarray, 
                    r_sample: np.ndarray,
                    x_legacy) -> np.ndarray:
    x_all = np.concatenate((x, x_legacy)) if len(x_legacy) else x
    return fun(x_all, acc_truth, r_sample)

def fun(x: np.ndarray, 
        acc_truth: np.ndarray, 
        r_sample: np.ndarray) -> float:
    x = x.reshape((-1,4))
    r_mascon, mu = x[:,:3], x[:,[-1]]
    mug = np.repeat(mu.reshape((-1,1)), repeats=r_sample.shape[0], axis=1)[...,np.newaxis]
    r_diff_norm, r_sample_grid, r_mascon_grid = pairwise_norm(r_sample, r_mascon, return_deep_grids=True)
    r_mascon_to_sample = r_sample_grid - r_mascon_grid
    acc_rec = np.sum(-mug * r_mascon_to_sample / r_diff_norm[:,:,np.newaxis]**3, axis=0)
    fx = np.sum(ps.vecnorm(acc_rec - acc_truth))
    return fx

ntest, alt_samples = 200, 1
ptest = np.empty((0,3))
for alt in np.linspace(300, 1000, alt_samples):
    alt_pts = (ps.AstroConstants.earth_r_eq + alt) * ps.spiral_sample_sphere(ntest // alt_samples)
    ptest = np.vstack((ptest, alt_pts))

x_legacy = []

max_order = 100
acc_truth = gravitas.acceleration(ptest, max_order, "EGM96")
f = []
mui = 10
for i in range(1000):
    print(f'Optimizing {i} mascons!')
    nmascon = 1
    mu0 = ps.AstroConstants.earth_mu / nmascon if i == 0 else mui
    x0 = np.hstack((ps.rand_point_in_ball(6378, nmascon), mu0 * np.ones((nmascon,1))))

    opt = scipy.optimize.minimize(lambda x: fun_with_legacy(x, acc_truth, ptest, x_legacy), x0.flatten(), method='BFGS')
    f.append(opt.fun)
    print(f'{opt.fun:.2e}')
    if i > 0 and f[-1] > f[-2]:
        mui *= 1
    x_legacy = np.concatenate((x_legacy, opt.x)) if len(x_legacy) else opt.x

optx = x_legacy.reshape((-1,4))
optr, optmu = optx[:,:3], optx[:,[-1]]

pl = pv.Plotter()
ps.scatter3(pl, optr, scalars=optmu)
ps.scatter3(pl, ptest, color='r')
ps.two_sphere(pl, radius=6378)
pl.camera.focal_point = (0.0, 0.0, 0.0)
pl.show()