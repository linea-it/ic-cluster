import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import fitting, models
from scipy.interpolate import interp1d
import copy

##Catalog Functions:
def reduced_catalog(catalog, seed, p_0, del_catalog = False):
    '''
    Reduce the number of rows in a catalog by randomly choosing rows.
    
    Args:
        catalog (array): the catalog to be reduced.
        seed (int)
    '''
    sq1 = np.random.SeedSequence(seed)
    rng = np.random.default_rng(sq1)
    random_indexs = rng.choice(2, size=len(catalog), p=[p_0,1-p_0]).astype('bool')
    catalog_reduced = catalog[random_indexs]
    
    if del_catalog:
        del catalog
    
    return catalog_reduced

## Astrophysics Functions:
def look_dist(z,size,option, sigma_m = 0.308, sigma_k = 0.0, sigma_l = 0.692, H_0 = 67.8, c = 299792.458):
    
    def definite_integral(f,lim_inf,lim_sup):
        return integrate.quad(lambda x: f(x), lim_inf, lim_sup)[0]

    def inv_E(z):
        return (sigma_m*((1+z)**3.0) + sigma_k*((1+z)**2.0) + sigma_l)**(-0.5)

    def com_dist_lof(z):
        return (c/H_0) * definite_integral(inv_E, 0, z)
    
    def com_dist_trans(z):
        if sigma_k == 0.0:
            return com_dist_lof(z)

        elif sigma_k > 0.0:
            return (c/H_0)*(sigma_k**(-0.5))*np.sinh((sigma_k**0.5)*com_dist_lof(z)/(c/H_0))

        elif sigma_k < 0.0:
            return (c/H_0)*(np.abs(sigma_k)**(-0.5))*np.sinh((np.abs(sigma_k)**0.5)*com_dist_lof(z)/(c/H_0))

    def ang_diam_dist(z):
        return com_dist_trans(z)/(1+z)
    
    D_A = ang_diam_dist(z)
    
    if option == 'degree':
        return np.degrees(size/D_A)
    
    elif option == 'arcmin':
        return 60*np.degrees(size/D_A)
    
    elif option == 'arcsec':
        return 3600*np.degrees(size/D_A)
    
    elif option == 'mpc':
        return np.radians(size/60)*D_A
    
def look_dist_w0wa(z,size,option, w_0 = -1.0 , w_a = 0.0, sigma_m = 0.3156, sigma_k = 0.0, H_0 = 67.27, c = 299792.458):
    
    def definite_integral(f,lim_inf,lim_sup):
        return integrate.quad(lambda x: f(x), lim_inf, lim_sup)[0]
    
    def w(z):
        return w_0 + (z/(1+z))*w_a
    
    def int_w(z):
        return (1+w(z))/(1+z)
    
    def inv_E(z):
        return (sigma_m*((1+z)**3.0) + sigma_k*((1+z)**2.0) + (1-sigma_m-sigma_k)*np.exp(3*definite_integral(int_w,0,z)))**(-0.5)

    def com_dist_lof(z):
        return (c/H_0) * definite_integral(inv_E, 0, z)
    
    def com_dist_trans(z):
        if sigma_k == 0.0:
            return com_dist_lof(z)

        elif sigma_k > 0.0:
            return (c/H_0)*(sigma_k**(-0.5))*np.sinh((sigma_k**0.5)*com_dist_lof(z)/(c/H_0))

        elif sigma_k < 0.0:
            return (c/H_0)*(np.abs(sigma_k)**(-0.5))*np.sin((np.abs(sigma_k)**0.5)*com_dist_lof(z)/(c/H_0))

    def ang_diam_dist(z):
        return com_dist_trans(z)/(1+z)
    
    D_A = ang_diam_dist(z)
    
    if option == 'degree':
        return np.degrees(size/D_A)
    
    elif option == 'arcmin':
        return 60*np.degrees(size/D_A)
    
    elif option == 'arcsec':
        return 3600*np.degrees(size/D_A)
    
    elif option == 'mpc':
        return np.radians(size/60)*D_A
    
##Plot Functions
def make_fig(nrows, ncols, figsize, titulo, show_up = True):
	if show_up == True:
		fig, axs = plt.subplots(nrows=nrows, ncols=ncols,constrained_layout=True, figsize = (figsize[0],figsize[1]))
	else:
		fig, axs = plt.subplots(nrows=nrows, ncols=ncols,constrained_layout=False,figsize = (figsize[0],figsize[1]))
	plt.suptitle(titulo, fontsize = 18, fontweight = 'bold')
	return fig, axs

## Ajustes
def continuum_black_body(x, y, output='polynomial', degree=3, n_iterate=7, lower_threshold=4, upper_threshold=4, verbose=False, weights=None):

    """
    Builds a polynomial continuum from segments of a spectrum,
    given in the form of wl and flux arrays.

    Parameters
    ----------
    x : array-like
        Independent variable
    y : array-like
        y = f(x)
    output: string
        Specifies what will be returned by the function

        'ratio'      = ratio between fitted continuum and the spectrum
        'difference' = difference between fitted continuum and the spectrum
        'function'   = continuum function evaluated at x

    degree : integer
        Degree of polynomial for the fit
    n_iterate : integer
        Number of rejection iterations
    lower_threshold : float
        Lower threshold for point rejection in units of standard
        deviation of the residuals
    upper_threshold : float
        Upper threshold for point rejection in units of standard
        deviation of the residuals
    verbose : boolean
        Prints information about the fitting
    weights : array-like
        Weights for continuum fitting. Must be the shape of x and y.

    Returns
    -------
    c : tuple

        c[0]: numpy.ndarray
            Input x coordinates
        c[1]: numpy.ndarray
            See parameter "output".

    """

    assert not np.isnan(x).all(), 'All x values are NaN.'
    assert not np.isnan(y).all(), 'All y values are NaN.'

    x_full = copy.deepcopy(x)
    # NOTE: For now, interp1d skips interpolation of NaNs.
    s = interp1d(x, y)

    if weights is None:
        weights = np.ones_like(x)

    if np.isnan(y).any():
        nan_mask = np.isnan(s(x))
        x = x[~nan_mask]
        weights = copy.deepcopy(weights)[~nan_mask]
        warnings.warn(
            'NaN values found in data! Removed {:d} out of {:d} data points.'.format(
                np.count_nonzero(nan_mask), len(x_full)),
            category=RuntimeWarning,
        )

    model = models.Legendre1D(degree=degree)
    fitter = fitting.LinearLSQFitter()

    for i in range(n_iterate):

        f = fitter(model, x, s(x), weights=weights)
        res = s(x) - f(x)
        sig = np.std(res)
        rej_cond = ((res < upper_threshold * sig) & (res > -lower_threshold * sig))

        if np.sum(rej_cond) <= degree:
            if verbose:
                warnings.warn('Not enough fitting points. Stopped at iteration {:d}. sig={:.2e}'.format(i, sig))
            break

        if np.sum(weights != 0.0) <= degree:
            if verbose:
                warnings.warn(
                    'Number of non-zero values in weights vector is lower than the polynomial degree. '
                    'Stopped at iteration {:d}. sig={:.2e}'.format(i, sig))
            break

        x = x[rej_cond]
        weights = weights[rej_cond]

    if verbose:
        print('Final number of points used in the fit: {:d}'.format(len(x)))
        print('Rejection ratio: {:.2f}'.format(1. - float(len(x)) / float(len(x_full))))

    p = fitter(model, x, s(x), weights=weights)

    out = {'difference': (x_full, s(x_full) - p(x_full)), 'function': (x_full, p(x_full)), 'polynomial': p}
    if all(p(x_full) == 0.0):
        #warnings.warn('Continuum is identically zero. Setting ratio to NaN.', stacklevel=2)
        nan_spec = np.empty_like(x_full)
        nan_spec[:] = np.nan
        out['ratio'] = (x_full, nan_spec)
    else:
        out['ratio'] = (x_full, s(x_full) / p(x_full))

    return out[output]
