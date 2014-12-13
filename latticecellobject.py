#! /usr/bin/env python
from random import random
# random() returns random number between 0 and 1
from random import choice
# choice(seq) psuedorandomly returns a value of in seq
from random import randint
# randint(a,b) psuedorandomly returns an integer in range (a,b) inclusive


def pm():
    """this function returns + or - 1 psuedorandomly"""
    return choice((-1,1))

class lattice_cell_object:
    """This class contains a state that is its current state, and can process several messages.
    It is meant to simulate a lattice site in a lattice model of the reaction of 3-HTMF
    with methyl cinnamate catalyzed by TADDOL. It uses Monte Carlo probabilities for
    pretty much everything
    init will set the state of the instance with some chance of having molecules present
    checkstate will return the state of the instance
    move will attempt to move anything to another randomly chosen adjacent instance
    react will attempt to react if all molecules are present
    acceptmove will see if this instance can accept a molecule or molecules
    excite will possibly excite some molecules that are present"""

    def __init__(self, moleprob=0.01, reaction_rate=0.1,
                 reaction_favoritism=2, assoc_stabilization=0.1,
                 assoc_favoritism=2,
                 dimensions=2):
        """This initializes an instance of lattice_cell_object with a chance of population
        of molecules determined by moleprob.
        Syntax: __init__(moleprob, reaction_rate,
        reaction_favoritism, assoc_stabilization,
        association_favoritism, dimensions)
        moleprob is an optional argument that must be between 0 and 1.
        It is the probability of a molecule starting in this lattice site.
        reaction_favoritism is the relative amount the
        positive product is favored in the reaction.
        assoc_stabilization is the amount association is favored
        over molecules by themselves.
        assoc_favoritism is the relative amount the
        positive product is favored in staying associated.
        dimensions is the dimensionality of the lattice. Probably going to stay at 2d.
        Returns None."""
        # Check values of inputs
        if type(moleprob) != float:
            raise TypeError('moleprob must be a real number (float)')
        elif moleprob >= 1:
            raise ValueError('moleprob must be less than 1')
        elif moleprob <= 0:
            raise ValueError('moleprob must be greater than 0')
        elif reaction_rate <= 0:
            raise ValueError('reaction rate must be greater than 0')
        elif (reaction_favoritism <= 0
              or reaction_favoritism >= 1):
            raise ValueError('reaction favoritism must be between 0 and 1')
        elif assoc_stabilization <= 0:
            raise ValueError('assoc stabilization must be greater than 0')
        elif (assoc_favoritism <= 0
              or assoc_favoritism >= 1):
            raise ValueError('assoc favoritism must be between 0 and 1')
        # With Monte Carlo-like probability sets occupancy for catalyst.
        if random() > moleprob:
            self.catalyst = 1
        else:
            self.catalyst = 0
        # With Monte Carlo-like probability sets occupancy for htmf.
        # If occupied, can be in pm 1 orientation
        if random() > moleprob:
            self.htmf = pm()
        else:
            self.htmf = 0
        # With Monte Carlo-like probability sets occupancy for cinnamate.
        # If occupied, can be in pm 1 orientation
        if random() > moleprob:
            self.cinna = pm()
        else:
            self.cinna = 0
        # sets product occupancy to 0, as well as excitation state
        self.product = 0
        self.htmf_excitation_state = 0
        self.catalyst_excitation_state = 0
        # Need to define relevant variables for this!!
        # They should be taken as input to init and set as internal variables here
        self.p_react_pos = reaction_rate * reaction_favoritism
        self.p_react_neg = reaction_rate / reaction_favoritism
        self.p_assoc_pos = assoc_stabilization * assoc_favoritism
        self.p_assoc_neg = assoc_stabilization / assoc_favoritism
        self.max_move    = 2 * dimensions

    def __repr__(self):
        """Returns current occupation of this instance"""
        return self.catalyst,self.htmf,self.cinna,self.product

    def move(self):
        """If there is anything in this lattice site,
        then it will return which direction it wants to move as an integer
        and its current state as:
        direction, catalyst, catalyst excitation, htmf, htmf excitation, cinna,
        product"""
        # square each occupation (to avoid problems with negatives) and sum
        # this will be used to see if anything is here
        sum_state_squared = self.catalyst**2 + self.htmf**2 + self.cinna**2 + self.product**2
        if sum_state_squared > 0:
            move_direc = randint(1,self.max_move)
            return (move_direc, self.catalyst, self.catalyst_excitation_state,
                    self.htmf, self.htmf_excitation_state, self.cinna, self.product)

    def excite(self):
