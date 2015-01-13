import mclatticecellobject as site
import numpy as np

class Lattice:
    """This class will create a lattice object with size
    (size)**dimension.
    It is able to apply a keyword to all sites in the lattice.
    Additionally, because it's needed by moving in the cells,
    it can take returned arguments after applying a keyword, apply them to
    a specified other cell, then give input back to the original mover."""

    def __init__(self, arguments, size = 3, dimension = 2):
        """This will create a lattice with edge length
        size, (and dimension dimension?).
        First argument in inputs to each lattice site.
        Other arguments are optional.
        syntax: Lattice(arguments, size, dimension)
        Returns None."""
        # input checks:
        if type(size) != int:
            raise TypeError("size must be an integer, given %s" % type(size))
        if size < 1:
            raise ValueError('size must be greater than 1')
        if type(dimension) != int:
            raise TypeError('dimension must be an integer')
        if dimension < 1:
            raise ValueError('dimension must be >= 1')
        if isinstance(arguments, basestring):
            raise TypeError('arguments must be a list-like, given %s' \
                            % type(arguments))
        if len(arguments) != 7:
            raise ValueError('arguments must be len 7')
        self.dimension = dimension
        self.size      = size
        #
        # declares self.object lattice as an empty numpy array
        # with object-type sites
        #
        if dimension == 2:
            self.ob_lattice = np.empty((size, size), dtype = object)
            for i in xrange(size):
                for j in xrange(size):
                    self.ob_lattice[i, j] = self.mapable_site(arguments)
            # moves is the list for picking move directions.
            # It will be referenced when a site tries to make a move,
            # and the iterator needs to interpret where it is trying
            # to make a move to.
            # It is the list of changes to the currently referenced site
            # to which it wants to move.
            # AKA, if it wants to move up, the indices will be changed by
            # (delta row = -1, delta column = 0).
            # The ordering here is arbitrary as it is just picked
            # randomly anyway and probabilities are checked afterwards.
            # To get the first part of the first item, reference as
            # self.moves[0][0]
            self.moves = ((-1, 0), (0, -1), (1, 0), (0, 1))
            # This sample_templ is a template for output.
            # It is an array of strings with length defined by the number
            # in the 'a10' at the end.
            # The length should be no greater than 10
            # because the very longest possible state returned should
            # be '1,-1,-1,-1' (though hopefully that won't ever happen).
            # If for any reason the length needs to be greater than that,
            # this length needs to be changed, or the output will be
            # truncated to 10 characters.
            self.sample_templ = np.empty((size, size), dtype = 'a10')
        else:
            raise ValueError('Lattice has not be taught to handle ' + \
                             'dimension != 2, given %i' % dimension)
        

    def mapable_site(self, arguments):
        """This function returns a cell object initialized with
        arguments as a list (tuple). It is used because it's neater
        to pass a single argument to create each lattice object.
        It was originally intended to be 'np.vectorize'd, but
        I couldn't figure out how to make an array of arguments
        and have this apply to each site. (Possibly need to define
        custom dtype for each site of the argument lattice.)
        Syntax mapable_site(arguments)
        Returns initialized Lattice Cell Object."""
        # get arguments from arguments tuple to pass to
        # site.Lattice_Cell_Object to initialize each lattice site
        molecprob           = float(arguments[0])
        reaction_energy     = float(arguments[1])
        reaction_favoritism = float(arguments[2])
        assoc_stabilization = float(arguments[3])
        assoc_favoritism    = float(arguments[4])
        excit_prob          = float(arguments[5])
        beta                = float(arguments[6])
        return site.Lattice_Cell_Object(molecprob, reaction_energy,
                                 reaction_favoritism, assoc_stabilization,
                                 assoc_favoritism, excit_prob,
                                 beta, self.dimension)

    def over_sites(self, kw):
        """This function will take a keyword and apply that keyword
        over all the lattice site. The keyword must belong to this
        list: (excite, react, move, sample).
        Syntax over_sites(keyword)
        Returns None, except if keyword == sample,
        then it will return the occupation state of each site"""
        #
        #
        # There might be a better way to do this: three things
        # I saw on the Python page on classes:
        # 1. define __iter__ and next in the class definition of the
        # lattice, then it should be easily iterable.
        # 2. use a generator
        # 3. use a generator expression
        # The main problem I see with these right now is that it might
        # be difficult to then reference wherever it is trying to
        # move to.
        # I don't know if there are any benefits. Might help
        # readability a little bit, but probably not a ton.
        #
        #
        if type(kw) != str:
            raise TypeError('keyword must be a string')
        # Simple for excite and react:
        # just apply the keyword to each lattice site.
        # We expect no output or anything. Each site just acts internally.
        #
        # Look into nditer for better iteration over arrays?
        #
        if   kw == 'excite':
            for i in xrange(self.size):
                for j in xrange(self.size):
                    self.ob_lattice[i, j].excite()
        elif kw == 'react':
            for i in xrange(self.size):
                for j in xrange(self.size):
                    self.ob_lattice[i, j].react()
        #
        # move is more difficult. If the site has something there,
        # it will attempt to move by returning its current state
        # and a direction it is trying to move. We need to catch that
        # information, pass it to where ever it is trying to move,
        # then take whatever that returns, and pass it back to the
        # original cell that attempted the move.
        #
        elif kw == 'move':
            for i in xrange(self.size):
                for j in xrange(self.size):
                    move_out = self.ob_lattice[i, j].move()
                    # If nothing to move, will return move_out = None, 
                    # which will evaluate to False.
                    # Otherwise, this will evaluate to True.
                    if move_out:
                        move_direc = self.moves[move_out[0]]
                        # Now, takes what is trying to be moved and
                        # gives it to the accept_move function of the
                        # site determined by move_direc(tion).
                        # The accept_move function will return whatever
                        # it does not accept, which needs to be then
                        # passed back to the move_result function
                        # of the original lattice site.
                        #
                        # Had some trouble indexing past edges when it wanted
                        # to move past the ends, so np.take wraps around edges
                        # but it also indexes like a flattened array.
                        # That's why the math is needed for the index.
                        # Not certain this is right, but as long as it is 
                        # unique it should work fine to converge to
                        # uniform sampling.
                        #
                        move_to = self.ob_lattice.take(\
                            self.size * (i + move_direc[0]) +\
                            (j + move_direc[1]), mode = 'wrap')
                        amove_out = move_to.accept_move(
                                            move_out[1], move_out[2],
                                            move_out[3], move_out[4],
                                            move_out[5], move_out[6])
                        self.ob_lattice[i, j].move_result(
                            amove_out[0], amove_out[1],
                            amove_out[2], amove_out[3],
                            amove_out[4], amove_out[5])
        elif kw == 'sample':
            sample_out = self.sample_templ
            for i in xrange(self.size):
                for j in xrange(self.size):
                    # This will call __repr__ of each site which will
                    # output a string of the current occupancy.
                    # This will then be added to the correct spot in
                    # the output file.
                    sample_out[i, j] = str(self.ob_lattice[i, j])
            return sample_out
        else:
            raise IOError('keyword not recognized. Given: %s' % kw)


