from Ranger import Range
from functools import partial
import numpy as np
from numba import jit
import numba as nb
import sys

if sys.version_info[0] > 2:
    xrange = range

@jit(argtypes=(nb.types.pyobject,
               nb.double,
               nb.double,
               nb.int64,
               nb.double,
               nb.boolean), restype=nb.double)
def _step_left_endpoint(obj_func, y, L, J, w, m):
    while y < obj_func(L):
        L -= w
        if m:
            J -= 1
            if J <= 0:
                break
    return L

@jit(argtypes=(nb.types.pyobject,
               nb.double,
               nb.double,
               nb.int64,
               nb.double,
               nb.boolean), restype=nb.double)
def _step_right_endpoint(obj_func, y, R, K, w, m):
    while y < obj_func(R):
        R += w
        if m:
            K -= 1
            if K <= 0:
                break
    return R

class slice_sampler(object):
    """ Generic slice-sampler.
    Base class for both univariate and multivariate methods

    See http://www.cs.toronto.edu/~radford/ftp/slc-samp.pdf for theory
    and details
    """
    def __init__(self, log_f):
        """ Instantiates a slice sampler

        Parameters
        ----------
        log_f : callable
            A function that returns the log of the function you want to
            sample and accepts a numpy array as an argument (the x)
        """
        if not hasattr(log_f, '__call__'):
            raise TypeError("log_f is not callable")
        ## The function to sample
        self._g = log_f
    def set_function(self, log_f):
        """ Sets the function being sampled

        Parameters
        ----------
        log_f : callable
            A function that returns the log of the function you want to
            sample and accepts a numpy array (or scalar if univariate function)
            as an argument (the x)
        """
        self._g = log_f

class univariate_slice_sampler(slice_sampler):
    """ Does slice sampling on univariate distributions
    """
    def __init__(self, log_f, adapt_w = False, start_w = 0.1):
        """ Instantiates a slice sampler

        Parameters
        ----------
        log_f : callable
            A function that returns the log of the function you want to
            sample and accepts a scalar as an argument (the x)
        adapt_w : boolean
            Whether to adapt w during sampling. Will work in between samplings
        start_w : float
            Starting value for w, necessary if adapting w during sampling
        """        
        super(univariate_slice_sampler, self).__init__(log_f)
        self._w = start_w
        self._adapt_w = adapt_w
        ## Parameter for number of data points going into w
        self._w_n = 1.
    def accept_doubling(self, x0, x1, y, w, interval):
        """ Test for whether a new point, x1, is an acceptable next state,
        when the interval was found by the doubling procedure
        (See Radford paper at http://www.cs.toronto.edu/~radford/ftp/slc-samp.pdf)

        Parameters
        ----------
        x0 : float
            The current point
        x1 : float
            The possible next point
        y : float
            The vertical level defining the slice
        w : float
            Estimte of typical slice size
        interval : Ranger.Range
            The interval found be the doubling procedure

        Returns
        -------
        Whether or not x1 is an acceptable next state
        """
        L_hat = interval.lowerEndpoint()
        R_hat = interval.upperEndpoint()
        D = False
        while (R_hat - L_hat) > 1.1*w:
            M = (L_hat+R_hat)/2.
            if (x0 < M and x1 >= M) or (x0 >= M and x1 < M):
                D = True
            if x1 < M:
                R_hat = M
            else:
                L_hat = M
            if D and (y > self._g(L_hat) and y >= self._g(R_hat)):
                return False
        return True
    def find_interval_doubling(self, x0, y, w, p=None):
        """ Finds an interval around a current point, x0, using a doubling procedure
        (See Radford paper at http://www.cs.toronto.edu/~radford/ftp/slc-samp.pdf)

        Parameters
        ----------
        xo : float
            The current point
        y : float
            The vertical level defining the slice
        w : float
            Estimate of typical slice size
        p : int, optional
            Integer limiting the size of a slice to (2^p)w. If None,
            then interval can grow without bound

        Returns
        -------
        Range containing the interval found
        """
        # Sample initial interval around x0
        U = np.random.random()
        L = x0-w*U
        R = L + w
        # Set limiter
        K = None
        if p:
            K = p
        # Refine interval
        while (y < self._g(L) or y < self._g(R)):
            V = np.random.random()
            if V < 0.5:
                L -= (R-L)
            else:
                R = R + (R-L)
            if p:
                K -= 1
                if K <= 0:
                    break
        # Return the inteval
        return Range.closed(L,R)
    def find_interval_step_out(self, x0, y, w, m=None):
        """ Finds an interval around a current point, x0, using the stepping-out
        procedure (See Radford paper at http://www.cs.toronto.edu/~radford/ftp/slc-samp.pdf)

        Parameters
        ----------
        xo : float
            The current point
        y : float
            The vertical level defining the slice
        w : float
            Estimate of typical slice size
        m : int, optional
            Integer, where maximum size of slice should be mw. If None,
            then interval can grow without bound

        Returns
        -------
        Range containing the interval found
        """
        # Sample initial interval around x0
        U = np.random.random()
        L = x0 - w*U
        R = L + w
        # Place limitations on interval if necessary
        V = None
        J = None
        K = None
        if m:
            V = np.random.random()
            J = np.floor(m*V)
            K = (m-1)-J
        # Get left endpoint
        L = _step_left_endpoint(self._g, y, L, J, w, m)
        # Get right endpoint
        R = _step_right_endpoint(self._g, y, R, K, w, m)
        # Return interval
        return Range.closed(L,R)
    def run_sampler(self, x0_start = 0., n_samp = 10000, interval_method='doubling',
                    w=0.1, m=None, p=None):
        """ Runs the slice sampler

        Parameters
        ----------
        x0_start : float
            An initial value for x
        n_samp : int
            The number of samples to take
        interval_method : str
            The method for determining the interval at each stage of sampling. Possible values
            are 'doubling', 'stepping'.
        w : float
            Estimate of typical slice size. If adapt_w is true, then this is overriden
        m : int, optional (Only relevant for stepping interval procedure)
            Integer, where maximum size of slice should be mw. If None,
            then interval can grow without bound.
        p : int, optional (Only relevant for doubling interval procedure)
            Integer limiting the size of a slice to (2^p)w. If None,
            then interval can grow without bound

        Returns
        -------
        Generator of samples from the distribution

        Examples
        --------
        >>> from quantgen.stats_utils.slice_sampler import univariate_slice_sampler
        >>> from scipy.stats import norm
        >>> sampler = univariate_slice_sampler(lambda x: norm.logpdf(x, loc=0, scale=5.))
        >>> samples = [x for x in sampler.run_sampler(n_samp=1000, w=2.5)]
        """
        x0 = x0_start
        interval = None
        doubling_used = True
        if interval_method != 'doubling':
            doubling_used = False
        if self._adapt_w:
            w = self._w
        for i in xrange(n_samp):
            # Draw vertical value, y, that defines the horizontal slice
            y = self._g(x0) - np.random.exponential()
            # Find interval around x0 that contains at least a big part of the slice
            if interval_method == 'doubling':
                interval = self.find_interval_doubling(x0,y,w=w,p=p)
            elif interval_method == 'stepping':
                interval = self.find_interval_step_out(x0,y,w=w,m=m)
            else:
                raise ValueError("%s is not an interval method" % interval_method)
            # Draw new point
            x0 = self.sample_by_shrinkage(x0, y, w, interval, doubling_used=doubling_used)
            # Update w if necessary
            if self._adapt_w:
                interval_length = (interval.upperEndpoint()-interval.lowerEndpoint())
                self._w = np.power(self._w,self._w_n/(self._w_n+1.))*\
                    np.power(interval_length/2.,1./(self._w_n+1.))
                self._w_n += 1.
                w = self._w            
            yield float(x0)
    def sample_by_shrinkage(self, x0, y, w, interval, doubling_used = False):
        """ Samples a point from an interval using the shrinkage procedure
        (See Radford paper at http://www.cs.toronto.edu/~radford/ftp/slc-samp.pdf)

        Parameters
        ----------
        x0 : float
            The current point
        y : float
            The vertical level defining the slice
        w : float
            Estimate of typical slice size
        interval : Ranger.Range
            The interval found be the doubling procedure
        doubling_used : boolean
            Whether the doubling procedure was used when defining the interval

        Returns
        -------
        The new point
        """
        # Set up the accept function, based on whether interval was found
        # using doubling
        accept_func = None
        if doubling_used:
            accept_func = lambda x0, x1, y, w, interval: \
              self.accept_doubling(x0, x1, y, w, interval)
        else:
            accept_func = lambda x0, x1, y, w, interval: True
        # Set initial interval
        L_bar = interval.lowerEndpoint()
        R_bar = interval.upperEndpoint()
        x1 = L_bar + 0.5*(R_bar-L_bar)
        # Run through shrinkage
        while 1:
            U = np.random.random()
            x1 = L_bar + U*(R_bar-L_bar)
            if y < self._g(x1) and accept_func(x0, x1, y, w, interval):
                break
            if x1 < x0:
                L_bar = x1
            else:
                R_bar = x1
        return x1

class multivariate_slice_sampler(slice_sampler):
    """ Run slice sampling on multivariate functions
    """
    def __init__(self, log_f, dim):
        """ Instantiates a slice sampler

        Parameters
        ----------
        log_f : callable
            A function that returns the log of the function you want to
            sample and accepts a numpy array as an argument (the x)
        dim : int
            The dimensionality of the x (e.g. 2 for a bivariate normal)
        """        
        super(multivariate_slice_sampler, self).__init__(log_f)
        self._dim = dim
    def get_univariates_mutable_sampler(self, x0_start, interval_method='doubling', w=0.1,
                                        m=None, p=None, adapt_w=False):
        """ Gets an instance of multivariate_mutable_slice_sampler, which can have its objective
        function easily changed between sampling steps

        Parameters
        ----------
        x0_start : np.ndarray
            An initial value for x
        n_samp : int
            The number of samples to take. This refers to the number of times
            the sampler cycles through ALL x_i
        interval_method : str
            The method for determining the interval at each stage of sampling. Possible values
            are 'doubling', 'stepping'.
        w : float
            Estimate of typical slice size
        m : int, optional (Only relevant for stepping interval procedure)
            Integer, where maximum size of slice should be mw. If None,
            then interval can grow without bound.
        p : int, optional (Only relevant for doubling interval procedure)
            Integer limiting the size of a slice to (2^p)w. If None,
            then interval can grow without bound
        adapt_w : boolean
            Whether to adapt w during sampling.

        Returns
        -------
        An instance of multivariate_mutable_slice_sampler
        """
        return multivariate_mutable_slice_sampler(self._g, self._dim, x0_start,
                                                  mode='univariates',
                                                  interval_method=interval_method,
                                                  w=w, m=m, p=p, adapt_w=adapt_w)
    def run_sampler_univariates(self, x0_start, n_samp = 10000,
                                interval_method='doubling', w=0.1, m=None,
                                p=None, adapt_w = False):
        """ Runs the slice sampler by cycling through univariate samplers
        of each x_i with the other x_i's held fixed (like a Gibbs sampler)

        Parameters
        ----------
        x0_start : np.ndarray
            An initial value for x
        n_samp : int
            The number of samples to take. This refers to the number of times
            the sampler cycles through ALL x_i
        interval_method : str
            The method for determining the interval at each stage of sampling. Possible values
            are 'doubling', 'stepping'.
        w : float
            Estimate of typical slice size
        m : int, optional (Only relevant for stepping interval procedure)
            Integer, where maximum size of slice should be mw. If None,
            then interval can grow without bound.
        p : int, optional (Only relevant for doubling interval procedure)
            Integer limiting the size of a slice to (2^p)w. If None,
            then interval can grow without bound
        adapt_w : boolean
            Whether to adapt w during sampling.

        Returns
        -------
        Generator of samples from the distribution

        Examples
        --------
        >>> from quantgen.stats_utils.slice_sampler import multivariate_slice_sampler
        >>> from scipy.stats import multivariate_normal
        >>> import matplotlib.pyplot as plt
        >>> import numpy as np
        >>>
        >>> cov = np.eye(2)
        >>> cov[0,1]=0.9
        >>> cov[1,0]=0.9
        >>> rv = multivariate_normal(np.zeros(2),cov)
        >>> sampler = multivariate_slice_sampler(rv.logpdf, 2)
        >>> samples = np.array([x for x in sampler.run_sampler_univariates(np.zeros(2))])        
        """
        x0 = np.array(x0_start)
        if len(x0) != self._dim:
            raise ValueError("x0_start is not of the correct dimension")
        ## Create a univariate sampler for each x_i, storing in a list
        samplers = []
        def func(x_i, adj_ind):
            mult_arr = np.ones(self._dim)
            mult_arr[adj_ind] = 0
            add_arr = np.zeros(self._dim)
            add_arr[adj_ind] = x_i
            return self._g(np.multiply(x0,mult_arr)+add_arr)        
        for i in xrange(self._dim):
            samplers.append(univariate_slice_sampler(partial(func, adj_ind=int(i)),
                                                     adapt_w=adapt_w))
        ## Run through sampling
        for i in xrange(n_samp):
            # Go through each of the univariate samplers in turn
            for j,sampler in enumerate(samplers):
                x0[j] = next(sampler.run_sampler(x0_start=x0[j], n_samp=1,
                        interval_method=interval_method, w=w, m=m, p=p))
            # Yield the current sample
            yield np.array(x0)
            

@jit(argtypes=(nb.double,
               nb.double[:],
               nb.int64,
               nb.types.pyobject,
               nb.double[:],
               nb.double[:]), restype=nb.double)
def _univariate_sampling_func(x_i, x0, adj_ind, obj_func, mult_arr, add_arr):
    add_arr[adj_ind] = x_i
    return obj_func(np.multiply(x0, mult_arr)+add_arr)

class multivariate_mutable_slice_sampler(slice_sampler):
    """ Multivariate slice sampler that allows updates to function in between
    calls
    """
    def __init__(self, log_f, dim, x0_start,  mode='univariates', interval_method='doubling',
                 w=0.1, m=None, p=None, adapt_w=False):
        """ Instantiates mutable a slice sampler

        Parameters
        ----------
        log_f : callable
            A function that returns the log of the function you want to
            sample and accepts a numpy array as an argument (the x). Note that this
            can be changed in between calls
        dim : int
            The dimensionality of the x (e.g. 2 for a bivariate normal)
        x0_start : np.ndarray
            An initial value for x
        mode : str
            The type of multivariate slice sampler being used. Only choice at the moment
            is 'univariates', where univariate sampling for each parameter is done
            in turn
        interval_method : str
            The method for determining the interval at each stage of sampling. Possible values
            are 'doubling', 'stepping'.
        w : float
            Estimate of typical slice size. Only initializes if using adapt_w
        m : int, optional (Only relevant for stepping interval procedure)
            Integer, where maximum size of slice should be mw. If None,
            then interval can grow without bound.
        p : int, optional (Only relevant for doubling interval procedure)
            Integer limiting the size of a slice to (2^p)w. If None,
            then interval can grow without bound
        adapt_w : boolean
            Whether to adapt w during sampling.        
        """
        super(multivariate_mutable_slice_sampler, self).__init__(log_f)
        self._mode = mode
        self._dim = dim
        self._mult_arr = np.ones(self._dim)
        self._add_arr = np.zeros(self._dim)
        self._interval_method=interval_method
        self._w=w
        self._m=m
        self._p=p
        self._x0 = x0_start
        ## Holds individual samplers for univariate sampling
        self._univariate_samplers = []
        ## Set up univariate samplers
        if mode == 'univariates':
            for i in xrange(self._dim):
                self._univariate_samplers.append(univariate_slice_sampler(partial(_univariate_sampling_func,
                                                                                  x0=x0_start, obj_func=log_f,
                                                                                  mult_arr=self._mult_arr,
                                                                                  add_arr=self._add_arr,
                                                                                  adj_ind=int(i)),
                                                                          adapt_w=adapt_w,
                                                                          start_w=w))
        else:
            raise ValueError("%s is not a valid mode" % str(mode))
    def sample_one(self):
        """ Runs through a single iteration of sampling

        Returns
        -------
        Array of sampled values
        """
        if self._mode == 'univariates':
            for j, sampler in enumerate(self._univariate_samplers):
                # Make the mult_arr and add_arr the correct values
                self._mult_arr[j] = 0.
                # Put in the correct objective function
                sampler.set_function(partial(_univariate_sampling_func,
                                             x0=self._x0,
                                             adj_ind=int(j),
                                             obj_func=self._g,
                                             mult_arr=self._mult_arr,
                                             add_arr=self._add_arr))
                self._x0[j] = next(sampler.run_sampler(x0_start=self._x0[j], n_samp=1,
                                                       interval_method=self._interval_method,
                                                       w=self._w, m=self._m, p=self._p))
                
                # Put mult_arr and add_arr back
                self._mult_arr[j] = 1.
                self._add_arr[j] = 0.
            # Return the sampled value
            return np.array(self._x0)
