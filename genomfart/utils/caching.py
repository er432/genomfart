from llist import dllist, dllistnode

class LRUcache(dict):
    """ A least recently used cache

    Examples
    --------
    >>> from genomfart.utils.caching import LRUcache
    >>>
    >>> retreiveFunc = lambda k: k.upper()
    >>> def disposeFunc(k,v):
    ...     print('%s Test' % v)
    ...     return
    >>>
    >>> a = LRUcache(retreiveFunc, delFunc=disposeFunc, maxsize=3)
    """
    def __init__(self, retrieveFunc, delFunc=None, maxsize=1000, *args, **kwargs):
        """ Instantiates the LRUcache

        Parameters
        ----------
        retrieveFunc : callable
            Function that takes a key as an argument and returns the value
            corresponding to that key
        delFunc : callable
            Function that takes key, value as an arguments and does something before
            deletion
        maxsize : int
            The maximum size the cache can keep before kicking keys out
        """
        if not hasattr(retrieveFunc, '__call__'):
            raise TypeError("retrieveFunc not callable")
        # Call the dictionary constructor
        super(LRUcache, self).__init__(*args, **kwargs)
        self._retrieveFunc = retrieveFunc
        self._delFunc = None
        if delFunc:
            self._delFunc = delFunc
        self._maxsize = maxsize
        # Create dictionary of key -> node object and doublelinkedList
        # of keys
        self._keyNodeDict = {}
        self._keyQueue = dllist()
        for key in self:
            self._keyNodeDict[key] = self._keyQueue.append(key)
    def __delitem__(self, key):
        """ Removes item from the cache, and performs any logic necessary prior
        to removal

        Parameters
        ----------
        key : hashable
            The key to use for indexing
        """
        if self._delFunc:
            self._delFunc(key, super(LRUcache,self).__getitem__(key))
        v = super(LRUcache, self).__getitem__(key)
        del v
        super(LRUcache, self).__setitem__(key, None)
        super(LRUcache, self).__delitem__(key)
    def __setitem__(self, key, val):
        """ Adds a new item to the cache

        Parameters
        ----------
        key : hashable
            The key to use for indexing the val
        val : object
            The value
        """
        if len(self._keyQueue) == self._maxsize:
            # Pop out the least recently used key/value pair if at the
            # maximum size
            pop_key = self._keyQueue.popleft()
            del self[pop_key]
            del self._keyNodeDict[pop_key]
        # Check if the key already in the cache, and then
        # determine how to act
        if key in self:
            # If key already in the cache, put it at the end of the queue
            removeNode = self._keyNodeDict[key]
            self._keyQueue.remove(removeNode)
            self._keyNodeDict[key] = self._keyQueue.append(removeNode)
            # Put it in the regular dictionary
            dict.__setitem__(self, key, val)
        else:
            self._keyNodeDict[key] = self._keyQueue.append(key)
            dict.__setitem__(self, key, val)
    def __getitem__(self, key):
        """ Gets an item from the cache

        Parameters
        ----------
        key : Hashable
            Key value
        """
        if key in self:
            # If key already present, move to the end of the queue
            removeNode = self._keyNodeDict[key]
            self._keyQueue.remove(removeNode)
            self._keyNodeDict[key] = self._keyQueue.append(removeNode)
        # Use the superclass function to return the value
        return dict.__getitem__(self, key)
    def __missing__(self, key):
        """ Called when there is a cache miss

        Finds and sets the item using the retrieveFunc supplied

        Parameters
        ----------
        key : Hashable

        Returns
        -------
        The value associated with the key
        """
        val = self._retrieveFunc(key)
        self[key] = val
        return dict.__getitem__(self, key)
