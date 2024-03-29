#!/usr/bin/env python3

import argparse
import numpy as np
import collections
from scipy import stats
from statsmodels.stats._knockoff import RegressionFDR
import pandas

__all__ = ['fdrcorrection', 'fdrcorrection_twostage', 'local_fdr',
		   'multipletests', 'NullDistribution', 'RegressionFDR']
# ==============================================
#
# Part 1: Multiple Tests and P-Value Correction
#
# ==============================================


def _ecdf(x):
	'''no frills empirical cdf used in fdrcorrection
	'''
	nobs = len(x)
	return np.arange(1,nobs+1)/float(nobs)

multitest_methods_names = {'b': 'Bonferroni',
						   's': 'Sidak',
						   'h': 'Holm',
						   'hs': 'Holm-Sidak',
						   'sh': 'Simes-Hochberg',
						   'ho': 'Hommel',
						   'fdr_bh': 'FDR Benjamini-Hochberg',
						   'fdr_by': 'FDR Benjamini-Yekutieli',
						   'fdr_tsbh': 'FDR 2-stage Benjamini-Hochberg',
						   'fdr_tsbky': 'FDR 2-stage Benjamini-Krieger-Yekutieli',
						   'fdr_gbs': 'FDR adaptive Gavrilov-Benjamini-Sarkar'
						   }

_alias_list = [['b', 'bonf', 'bonferroni'],
			   ['s', 'sidak'],
			   ['h', 'holm'],
			   ['hs', 'holm-sidak'],
			   ['sh', 'simes-hochberg'],
			   ['ho', 'hommel'],
			   ['fdr_bh', 'fdr_i', 'fdr_p', 'fdri', 'fdrp'],
			   ['fdr_by', 'fdr_n', 'fdr_c', 'fdrn', 'fdrcorr'],
			   ['fdr_tsbh', 'fdr_2sbh'],
			   ['fdr_tsbky', 'fdr_2sbky', 'fdr_twostage'],
			   ['fdr_gbs']
			   ]


multitest_alias = {}
for m in _alias_list:
	multitest_alias[m[0]] = m[0]
	for a in m[1:]:
		multitest_alias[a] = m[0]

def multipletests(pvals, alpha=0.05, method='hs', is_sorted=False,
				  returnsorted=False):
	"""
	Test results and p-value correction for multiple tests

	Parameters
	----------
	pvals : array_like, 1-d
		uncorrected p-values.   Must be 1-dimensional.
	alpha : float
		FWER, family-wise error rate, e.g. 0.1
	method : str
		Method used for testing and adjustment of pvalues. Can be either the
		full name or initial letters. Available methods are:

		- `bonferroni` : one-step correction
		- `sidak` : one-step correction
		- `holm-sidak` : step down method using Sidak adjustments
		- `holm` : step-down method using Bonferroni adjustments
		- `simes-hochberg` : step-up method  (independent)
		- `hommel` : closed method based on Simes tests (non-negative)
		- `fdr_bh` : Benjamini/Hochberg  (non-negative)
		- `fdr_by` : Benjamini/Yekutieli (negative)
		- `fdr_tsbh` : two stage fdr correction (non-negative)
		- `fdr_tsbky` : two stage fdr correction (non-negative)

	is_sorted : bool
		If False (default), the p_values will be sorted, but the corrected
		pvalues are in the original order. If True, then it assumed that the
		pvalues are already sorted in ascending order.
	returnsorted : bool
		 not tested, return sorted p-values instead of original sequence

	Returns
	-------
	reject : ndarray, boolean
		true for hypothesis that can be rejected for given alpha
	pvals_corrected : ndarray
		p-values corrected for multiple tests
	alphacSidak : float
		corrected alpha for Sidak method
	alphacBonf : float
		corrected alpha for Bonferroni method

	Notes
	-----
	There may be API changes for this function in the future.

	Except for 'fdr_twostage', the p-value correction is independent of the
	alpha specified as argument. In these cases the corrected p-values
	can also be compared with a different alpha. In the case of 'fdr_twostage',
	the corrected p-values are specific to the given alpha, see
	``fdrcorrection_twostage``.

	The 'fdr_gbs' procedure is not verified against another package, p-values
	are derived from scratch and are not derived in the reference. In Monte
	Carlo experiments the method worked correctly and maintained the false
	discovery rate.

	All procedures that are included, control FWER or FDR in the independent
	case, and most are robust in the positively correlated case.

	`fdr_gbs`: high power, fdr control for independent case and only small
	violation in positively correlated case

	**Timing**:

	Most of the time with large arrays is spent in `argsort`. When
	we want to calculate the p-value for several methods, then it is more
	efficient to presort the pvalues, and put the results back into the
	original order outside of the function.

	Method='hommel' is very slow for large arrays, since it requires the
	evaluation of n partitions, where n is the number of p-values.
	"""
	import gc
	pvals = np.asarray(pvals)
	alphaf = alpha  # Notation ?

	if not is_sorted:
		sortind = np.argsort(pvals)
		pvals = np.take(pvals, sortind)

	ntests = len(pvals)
	alphacSidak = 1 - np.power((1. - alphaf), 1./ntests)
	alphacBonf = alphaf / float(ntests)
	if method.lower() in ['b', 'bonf', 'bonferroni']:
		reject = pvals <= alphacBonf
		pvals_corrected = pvals * float(ntests)

	elif method.lower() in ['s', 'sidak']:
		reject = pvals <= alphacSidak
		pvals_corrected = 1 - np.power((1. - pvals), ntests)

	elif method.lower() in ['hs', 'holm-sidak']:
		alphacSidak_all = 1 - np.power((1. - alphaf),
									   1./np.arange(ntests, 0, -1))
		notreject = pvals > alphacSidak_all
		del alphacSidak_all

		nr_index = np.nonzero(notreject)[0]
		if nr_index.size == 0:
			# nonreject is empty, all rejected
			notrejectmin = len(pvals)
		else:
			notrejectmin = np.min(nr_index)
		notreject[notrejectmin:] = True
		reject = ~notreject
		del notreject

		pvals_corrected_raw = 1 - np.power((1. - pvals),
										   np.arange(ntests, 0, -1))
		pvals_corrected = np.maximum.accumulate(pvals_corrected_raw)
		del pvals_corrected_raw

	elif method.lower() in ['h', 'holm']:
		notreject = pvals > alphaf / np.arange(ntests, 0, -1)
		nr_index = np.nonzero(notreject)[0]
		if nr_index.size == 0:
			# nonreject is empty, all rejected
			notrejectmin = len(pvals)
		else:
			notrejectmin = np.min(nr_index)
		notreject[notrejectmin:] = True
		reject = ~notreject
		pvals_corrected_raw = pvals * np.arange(ntests, 0, -1)
		pvals_corrected = np.maximum.accumulate(pvals_corrected_raw)
		del pvals_corrected_raw
		gc.collect()

	elif method.lower() in ['sh', 'simes-hochberg']:
		alphash = alphaf / np.arange(ntests, 0, -1)
		reject = pvals <= alphash
		rejind = np.nonzero(reject)
		if rejind[0].size > 0:
			rejectmax = np.max(np.nonzero(reject))
			reject[:rejectmax] = True
		pvals_corrected_raw = np.arange(ntests, 0, -1) * pvals
		pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
		del pvals_corrected_raw

	elif method.lower() in ['ho', 'hommel']:
		# we need a copy because we overwrite it in a loop
		a = pvals.copy()
		for m in range(ntests, 1, -1):
			cim = np.min(m * pvals[-m:] / np.arange(1,m+1.))
			a[-m:] = np.maximum(a[-m:], cim)
			a[:-m] = np.maximum(a[:-m], np.minimum(m * pvals[:-m], cim))
		pvals_corrected = a
		reject = a <= alphaf

	elif method.lower() in ['fdr_bh', 'fdr_i', 'fdr_p', 'fdri', 'fdrp']:
		# delegate, call with sorted pvals
		reject, pvals_corrected = fdrcorrection(pvals, alpha=alpha,
												 method='indep',
												 is_sorted=True)
	elif method.lower() in ['fdr_by', 'fdr_n', 'fdr_c', 'fdrn', 'fdrcorr']:
		# delegate, call with sorted pvals
		reject, pvals_corrected = fdrcorrection(pvals, alpha=alpha,
												 method='n',
												 is_sorted=True)
	elif method.lower() in ['fdr_tsbky', 'fdr_2sbky', 'fdr_twostage']:
		# delegate, call with sorted pvals
		reject, pvals_corrected = fdrcorrection_twostage(pvals, alpha=alpha,
														 method='bky',
														 is_sorted=True)[:2]
	elif method.lower() in ['fdr_tsbh', 'fdr_2sbh']:
		# delegate, call with sorted pvals
		reject, pvals_corrected = fdrcorrection_twostage(pvals, alpha=alpha,
														 method='bh',
														 is_sorted=True)[:2]

	elif method.lower() in ['fdr_gbs']:
		#adaptive stepdown in Gavrilov, Benjamini, Sarkar, Annals of Statistics 2009
##        notreject = pvals > alphaf / np.arange(ntests, 0, -1) #alphacSidak
##        notrejectmin = np.min(np.nonzero(notreject))
##        notreject[notrejectmin:] = True
##        reject = ~notreject

		ii = np.arange(1, ntests + 1)
		q = (ntests + 1. - ii)/ii * pvals / (1. - pvals)
		pvals_corrected_raw = np.maximum.accumulate(q) #up requirementd

		pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
		del pvals_corrected_raw
		reject = pvals_corrected <= alpha

	else:
		raise ValueError('method not recognized')

	if pvals_corrected is not None: #not necessary anymore
		pvals_corrected[pvals_corrected>1] = 1
	if is_sorted or returnsorted:
		return reject, pvals_corrected, alphacSidak, alphacBonf
	else:
		pvals_corrected_ = np.empty_like(pvals_corrected)
		pvals_corrected_[sortind] = pvals_corrected
		del pvals_corrected
		reject_ = np.empty_like(reject)
		reject_[sortind] = reject
		return reject_, pvals_corrected_, alphacSidak, alphacBonf



def fdrcorrection(pvals, alpha=0.05, method='indep', is_sorted=False):
	'''pvalue correction for false discovery rate

	This covers Benjamini/Hochberg for independent or positively correlated and
	Benjamini/Yekutieli for general or negatively correlated tests. Both are
	available in the function multipletests, as method=`fdr_bh`, resp. `fdr_by`.

	Parameters
	----------
	pvals : array_like
		set of p-values of the individual tests.
	alpha : float
		error rate
	method : {'indep', 'negcorr'}
	is_sorted : bool
		If False (default), the p_values will be sorted, but the corrected
		pvalues are in the original order. If True, then it assumed that the
		pvalues are already sorted in ascending order.

	Returns
	-------
	rejected : ndarray, bool
		True if a hypothesis is rejected, False if not
	pvalue-corrected : ndarray
		pvalues adjusted for multiple hypothesis testing to limit FDR

	Notes
	-----

	If there is prior information on the fraction of true hypothesis, then alpha
	should be set to alpha * m/m_0 where m is the number of tests,
	given by the p-values, and m_0 is an estimate of the true hypothesis.
	(see Benjamini, Krieger and Yekuteli)

	The two-step method of Benjamini, Krieger and Yekutiel that estimates the number
	of false hypotheses will be available (soon).

	Method names can be abbreviated to first letter, 'i' or 'p' for fdr_bh and 'n' for
	fdr_by.



	'''
	pvals = np.asarray(pvals)

	if not is_sorted:
		pvals_sortind = np.argsort(pvals)
		pvals_sorted = np.take(pvals, pvals_sortind)
	else:
		pvals_sorted = pvals  # alias

	if method in ['i', 'indep', 'p', 'poscorr']:
		ecdffactor = _ecdf(pvals_sorted)
	elif method in ['n', 'negcorr']:
		cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))   #corrected this
		ecdffactor = _ecdf(pvals_sorted) / cm
##    elif method in ['n', 'negcorr']:
##        cm = np.sum(np.arange(len(pvals)))
##        ecdffactor = ecdf(pvals_sorted)/cm
	else:
		raise ValueError('only indep and negcorr implemented')
	reject = pvals_sorted <= ecdffactor*alpha
	if reject.any():
		rejectmax = max(np.nonzero(reject)[0])
		reject[:rejectmax] = True

	pvals_corrected_raw = pvals_sorted / ecdffactor
	pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
	del pvals_corrected_raw
	pvals_corrected[pvals_corrected>1] = 1
	if not is_sorted:
		pvals_corrected_ = np.empty_like(pvals_corrected)
		pvals_corrected_[pvals_sortind] = pvals_corrected
		del pvals_corrected
		reject_ = np.empty_like(reject)
		reject_[pvals_sortind] = reject
		return reject_, pvals_corrected_
	else:
		return reject, pvals_corrected



def fdrcorrection_twostage(pvals, alpha=0.05, method='bky', iter=False,
						   is_sorted=False):
	'''(iterated) two stage linear step-up procedure with estimation of number of true
	hypotheses

	Benjamini, Krieger and Yekuteli, procedure in Definition 6

	Parameters
	----------
	pvals : array_like
		set of p-values of the individual tests.
	alpha : float
		error rate
	method : {'bky', 'bh')
		see Notes for details

		* 'bky' - implements the procedure in Definition 6 of Benjamini, Krieger
		   and Yekuteli 2006
		* 'bh' - the two stage method of Benjamini and Hochberg

	iter : bool

	Returns
	-------
	rejected : ndarray, bool
		True if a hypothesis is rejected, False if not
	pvalue-corrected : ndarray
		pvalues adjusted for multiple hypotheses testing to limit FDR
	m0 : int
		ntest - rej, estimated number of true hypotheses
	alpha_stages : list of floats
		A list of alphas that have been used at each stage

	Notes
	-----
	The returned corrected p-values are specific to the given alpha, they
	cannot be used for a different alpha.

	The returned corrected p-values are from the last stage of the fdr_bh
	linear step-up procedure (fdrcorrection0 with method='indep') corrected
	for the estimated fraction of true hypotheses.
	This means that the rejection decision can be obtained with
	``pval_corrected <= alpha``, where ``alpha`` is the original significance
	level.
	(Note: This has changed from earlier versions (<0.5.0) of statsmodels.)

	BKY described several other multi-stage methods, which would be easy to implement.
	However, in their simulation the simple two-stage method (with iter=False) was the
	most robust to the presence of positive correlation

	TODO: What should be returned?

	'''
	pvals = np.asarray(pvals)

	if not is_sorted:
		pvals_sortind = np.argsort(pvals)
		pvals = np.take(pvals, pvals_sortind)

	ntests = len(pvals)
	if method == 'bky':
		fact = (1.+alpha)
		alpha_prime = alpha / fact
	elif method == 'bh':
		fact = 1.
		alpha_prime = alpha
	else:
		raise ValueError("only 'bky' and 'bh' are available as method")

	alpha_stages = [alpha_prime]
	rej, pvalscorr = fdrcorrection(pvals, alpha=alpha_prime, method='indep',
								   is_sorted=True)
	r1 = rej.sum()
	if (r1 == 0) or (r1 == ntests):
		return rej, pvalscorr * fact, ntests - r1, alpha_stages
	ri_old = r1

	while True:
		ntests0 = 1.0 * ntests - ri_old
		alpha_star = alpha_prime * ntests / ntests0
		alpha_stages.append(alpha_star)
		#print ntests0, alpha_star
		rej, pvalscorr = fdrcorrection(pvals, alpha=alpha_star, method='indep',
									   is_sorted=True)
		ri = rej.sum()
		if (not iter) or ri == ri_old:
			break
		elif ri < ri_old:
			# prevent cycles and endless loops
			raise RuntimeError(" oops - should not be here")
		ri_old = ri

	# make adjustment to pvalscorr to reflect estimated number of Non-Null cases
	# decision is then pvalscorr < alpha  (or <=)
	pvalscorr *= ntests0 * 1.0 /  ntests
	if method == 'bky':
		pvalscorr *= (1. + alpha)

	if not is_sorted:
		pvalscorr_ = np.empty_like(pvalscorr)
		pvalscorr_[pvals_sortind] = pvalscorr
		del pvalscorr
		reject = np.empty_like(rej)
		reject[pvals_sortind] = rej
		return reject, pvalscorr_, ntests - ri, alpha_stages
	else:
		return rej, pvalscorr, ntests - ri, alpha_stages



def local_fdr(zscores, null_proportion=1.0, null_pdf=None, deg=7,
			  nbins=30, alpha=0):
	"""
	Calculate local FDR values for a list of Z-scores.

	Parameters
	----------
	zscores : array_like
		A vector of Z-scores
	null_proportion : float
		The assumed proportion of true null hypotheses
	null_pdf : function mapping reals to positive reals
		The density of null Z-scores; if None, use standard normal
	deg : int
		The maximum exponent in the polynomial expansion of the
		density of non-null Z-scores
	nbins : int
		The number of bins for estimating the marginal density
		of Z-scores.
	alpha : float
		Use Poisson ridge regression with parameter alpha to estimate
		the density of non-null Z-scores.

	Returns
	-------
	fdr : array_like
		A vector of FDR values

	References
	----------
	B Efron (2008).  Microarrays, Empirical Bayes, and the Two-Groups
	Model.  Statistical Science 23:1, 1-22.

	Examples
	--------
	Basic use (the null Z-scores are taken to be standard normal):

	>>> from statsmodels.stats.multitest import local_fdr
	>>> import numpy as np
	>>> zscores = np.random.randn(30)
	>>> fdr = local_fdr(zscores)

	Use a Gaussian null distribution estimated from the data:

	>>> null = EmpiricalNull(zscores)
	>>> fdr = local_fdr(zscores, null_pdf=null.pdf)
	"""

	from statsmodels.genmod.generalized_linear_model import GLM
	from statsmodels.genmod.generalized_linear_model import families
	from statsmodels.regression.linear_model import OLS

	# Bins for Poisson modeling of the marginal Z-score density
	minz = min(zscores)
	maxz = max(zscores)
	bins = np.linspace(minz, maxz, nbins)

	# Bin counts
	zhist = np.histogram(zscores, bins)[0]

	# Bin centers
	zbins = (bins[:-1] + bins[1:]) / 2

	# The design matrix at bin centers
	dmat = np.vander(zbins, deg + 1)

	# Rescale the design matrix
	sd = dmat.std(0)
	ii = sd >1e-8
	dmat[:, ii] /= sd[ii]

	start = OLS(np.log(1 + zhist), dmat).fit().params

	# Poisson regression
	if alpha > 0:
		md = GLM(zhist, dmat, family=families.Poisson()).fit_regularized(L1_wt=0, alpha=alpha, start_params=start)
	else:
		md = GLM(zhist, dmat, family=families.Poisson()).fit(start_params=start)

	# The design matrix for all Z-scores
	dmat_full = np.vander(zscores, deg + 1)
	dmat_full[:, ii] /= sd[ii]

	# The height of the estimated marginal density of Z-scores,
	# evaluated at every observed Z-score.
	fz = md.predict(dmat_full) / (len(zscores) * (bins[1] - bins[0]))

	# The null density.
	if null_pdf is None:
		f0 = np.exp(-0.5 * zscores**2) / np.sqrt(2 * np.pi)
	else:
		f0 = null_pdf(zscores)

	# The local FDR values
	fdr = null_proportion * f0 / fz

	fdr = np.clip(fdr, 0, 1)

	return fdr



class NullDistribution(object):
	"""
	Estimate a Gaussian distribution for the null Z-scores.

	The observed Z-scores consist of both null and non-null values.
	The fitted distribution of null Z-scores is Gaussian, but may have
	non-zero mean and/or non-unit scale.

	Parameters
	----------
	zscores : array_like
		The observed Z-scores.
	null_lb : float
		Z-scores between `null_lb` and `null_ub` are all considered to be
		true null hypotheses.
	null_ub : float
		See `null_lb`.
	estimate_mean : bool
		If True, estimate the mean of the distribution.  If False, the
		mean is fixed at zero.
	estimate_scale : bool
		If True, estimate the scale of the distribution.  If False, the
		scale parameter is fixed at 1.
	estimate_null_proportion : bool
		If True, estimate the proportion of true null hypotheses (i.e.
		the proportion of z-scores with expected value zero).  If False,
		this parameter is fixed at 1.

	Attributes
	----------
	mean : float
		The estimated mean of the empirical null distribution
	sd : float
		The estimated standard deviation of the empirical null distribution
	null_proportion : float
		The estimated proportion of true null hypotheses among all hypotheses

	References
	----------
	B Efron (2008).  Microarrays, Empirical Bayes, and the Two-Groups
	Model.  Statistical Science 23:1, 1-22.

	Notes
	-----
	See also:

	http://nipy.org/nipy/labs/enn.html#nipy.algorithms.statistics.empirical_pvalue.NormalEmpiricalNull.fdr
	"""

	def __init__(self, zscores, null_lb=-1, null_ub=1, estimate_mean=True,
				 estimate_scale=True, estimate_null_proportion=False):

		# Extract the null z-scores
		ii = np.flatnonzero((zscores >= null_lb) & (zscores <= null_ub))
		if len(ii) == 0:
			raise RuntimeError("No Z-scores fall between null_lb and null_ub")
		zscores0 = zscores[ii]

		# Number of Z-scores, and null Z-scores
		n_zs, n_zs0 = len(zscores), len(zscores0)

		# Unpack and transform the parameters to the natural scale, hold
		# parameters fixed as specified.
		def xform(params):

			mean = 0.
			sd = 1.
			prob = 1.

			ii = 0
			if estimate_mean:
				mean = params[ii]
				ii += 1
			if estimate_scale:
				sd = np.exp(params[ii])
				ii += 1
			if estimate_null_proportion:
				prob = 1 / (1 + np.exp(-params[ii]))

			return mean, sd, prob


		from scipy.stats.distributions import norm


		def fun(params):
			"""
			Negative log-likelihood of z-scores.

			The function has three arguments, packed into a vector:

			mean : location parameter
			logscale : log of the scale parameter
			logitprop : logit of the proportion of true nulls

			The implementation follows section 4 from Efron 2008.
			"""

			d, s, p = xform(params)

			# Mass within the central region
			central_mass = (norm.cdf((null_ub - d) / s) -
							norm.cdf((null_lb - d) / s))

			# Probability that a Z-score is null and is in the central region
			cp = p * central_mass

			# Binomial term
			rval = n_zs0 * np.log(cp) + (n_zs - n_zs0) * np.log(1 - cp)

			# Truncated Gaussian term for null Z-scores
			zv = (zscores0 - d) / s
			rval += np.sum(-zv**2 / 2) - n_zs0 * np.log(s)
			rval -= n_zs0 * np.log(central_mass)

			return -rval


		# Estimate the parameters
		from scipy.optimize import minimize
		# starting values are mean = 0, scale = 1, p0 ~ 1
		mz = minimize(fun, np.r_[0., 0, 3], method="Nelder-Mead")
		mean, sd, prob = xform(mz['x'])

		self.mean = mean
		self.sd = sd
		self.null_proportion = prob


	# The fitted null density function
def pdf(self, zscores):
		"""
		Evaluates the fitted empirical null Z-score density.

		Parameters
		----------
		zscores : scalar or array_like
			The point or points at which the density is to be
			evaluated.

		Returns
		-------
		The empirical null Z-score density evaluated at the given
		points.
		"""

		zval = (zscores - self.mean) / self.sd
		return np.exp(-0.5*zval**2 - np.log(self.sd) - 0.5*np.log(2*np.pi))

def read_wig(wig_fyle):
	#create list of coverage values, position dependent
	return_dict = {}
	with open(wig_fyle,"r") as inp:
		for line in inp:
			return_dict[int(line.strip().split('\t')[0])] = [int(line.strip().split('\t')[1]),float(line.strip().split('\t')[2])]
	return return_dict

def read_pauses(fyle):
	return_dict = {}
	with open(fyle,'r') as inp:
		firstline = inp.readline()
		for line in inp:
			
			peak = int(line.strip().split(',')[0])
			strand = line.strip().split(',')[1]
	
			return_dict[peak] = strand

	return return_dict

def pause_compile(pauses,wig,mode):
	return_dict = {}
	if mode == 'Score':
		for key,value in pauses.items():
			count = wig[key][0]
			score = wig[key][1]

			if value == '+':
				if score > 0:
					return_dict[key] = abs(score)
				else:
					return_dict[key] = 0
				

			else:
				if score < 0:
					return_dict[key] = abs(score)
				else:
					return_dict[key] = 0
			
	elif mode == 'Count':
		for key,value in pauses.items():
			count = wig[key][0]
			score = wig[key][1]
			
			if value == '+':
				if score > 0:
					return_dict[key] = abs(count)
				else:
					return_dict[key] = 0
				

			else:
				if score < 0:
					return_dict[key] = abs(count)
				else:
					return_dict[key] = 0
			

	return return_dict

def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('pauses',type=str,help='')
	parser.add_argument('WT_1',type=str,help='')
	parser.add_argument('WT_2',type=str,help='')
	parser.add_argument('WT_3',type=str,help='')
	parser.add_argument('Mutant_1',type=str,help='')
	parser.add_argument('Mutant_2',type=str,help='')
	parser.add_argument('Mutant_3',type=str,help='')
	parser.add_argument('mode',type=str,help='Count or Score')
	parser.add_argument('output',type=str,help='')
	
	args = parser.parse_args()

	pauses = read_pauses(args.pauses)
	WT_1 = read_wig(args.WT_1)
	WT_2 = read_wig(args.WT_2)
	WT_3 = read_wig(args.WT_3)
	Mutant_1 = read_wig(args.Mutant_1)
	Mutant_2 = read_wig(args.Mutant_2)
	Mutant_3 = read_wig(args.Mutant_3)

	if args.mode == 'Count':
		WT_1_Count = pause_compile(pauses,WT_1,'Count')
		WT_2_Count = pause_compile(pauses,WT_2,'Count')
		WT_3_Count = pause_compile(pauses,WT_3,'Count')
		Mutant_1_Count = pause_compile(pauses,Mutant_1,'Count')
		Mutant_2_Count = pause_compile(pauses,Mutant_2,'Count')
		Mutant_3_Count = pause_compile(pauses,Mutant_3,'Count')

		WT_1,WT_2,WT_3,Mutant_1,Mutant_2,Mutant_3 = {},{},{},{},{},{}

		compiled = {}

		for key in pauses.keys():
			temp = [WT_1_Count[key],WT_2_Count[key],WT_3_Count[key],Mutant_1_Count[key],Mutant_2_Count[key],Mutant_3_Count[key]]
			if 'NA' not in temp:
				compiled[key] = temp

		compiled_sorted = collections.OrderedDict(sorted(compiled.items()))

		with open('Count_'+args.output,'w') as outp:
			outp.write('ID,WT_1,WT_2,WT_3,Mutant_1,Mutant_2,Mutant_3\n')
			for key,value in compiled_sorted.items():
				outp.write(str(key)+','+','.join(list(map(str,value)))+'\n')

	elif args.mode == 'Score':
		WT_1_Score = pause_compile(pauses,WT_1,'Score')
		WT_2_Score = pause_compile(pauses,WT_2,'Score')
		WT_3_Score = pause_compile(pauses,WT_3,'Score')
		Mutant_1_Score = pause_compile(pauses,Mutant_1,'Score')
		Mutant_2_Score = pause_compile(pauses,Mutant_2,'Score')
		Mutant_3_Score = pause_compile(pauses,Mutant_3,'Score')

		WT_1,WT_2,WT_3,Mutant_1,Mutant_2,Mutant_3 = {},{},{},{},{},{}

		compiled = {}

		for key in pauses.keys():
			temp = [WT_1_Score[key],WT_2_Score[key],WT_3_Score[key],Mutant_1_Score[key],Mutant_2_Score[key],Mutant_3_Score[key]]
			if 'NA' not in temp:
				compiled[key] = temp

		compiled_sorted = collections.OrderedDict(sorted(compiled.items()))

		p_val,FC = [],[]

		for key,value in compiled_sorted.items():
			WT = [value[0],value[1],value[2]]
			Mutant = [value[3],value[4],value[5]]
			p_val.append(stats.ttest_ind(WT,Mutant)[1])
			FC.append(np.log2(np.median(Mutant)/np.median(WT)))

		p_corrected = multipletests(p_val, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
		p_corrected = list(p_corrected)

		q = 0
		with open('Score_'+args.output,'w') as outp:
			outp.write('coord,Log2FC,pval,FDR\n')
			for key in compiled_sorted.keys():
				outp.write(str(key)+','+str(FC[q])+','+str(p_val[q])+','+str(p_corrected[q])+'\n')
				q += 1 

		with open('Score_raw_'+args.output,'w') as outp:
			outp.write('coord,WT_1,WT_2,WT_3,Mutant_1,Mutant_2,Mutant3\n')
			for key,value in compiled_sorted.items():
				outp.write(str(key)+','+','.join(list(map(str,value)))+'\n')

if __name__ == '__main__':
	main()