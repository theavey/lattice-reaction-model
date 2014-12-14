#! /usr/bin/env python
from random import random
# random() returns random number between 0 and 1
from random import choice
# choice(seq) psuedorandomly returns a value of in seq
from random import randint
# randint(a, b) psuedorandomly returns an integer in range (a, b) inclusive


def pm():
    """this function returns + or - 1 psuedorandomly"""
    return choice((-1, 1))

class lattice_cell_object:
    """This class contains a state that is its current state,
    and can process several messages.
    It is meant to simulate a lattice site in a lattice
    model of the reaction of 3-HTMF
    with methyl cinnamate catalyzed by TADDOL.
    It uses Monte Carlo probabilities for
    pretty much everything
    init will set the state of the instance with some chance of
    having molecules present
    checkstate will return the state of the instance
    move will attempt to move anything to another randomly
    chosen adjacent instance
    react will attempt to react if all molecules are present
    acceptmove will see if this instance can accept a molecule
    or molecules
    excite will possibly excite some molecules that are present"""

    def __init__(self, molecprob = 0.01, reaction_rate = 0.1,
                 reaction_favoritism = 2.0, assoc_stabilization = 0.1,
                 assoc_favoritism = 2.0, excit_prob = 0.01,
                 move_prob = 0.5
                 dimensions = 2):
        """This initializes an instance of lattice_cell_object with a
        chance of population
        of molecules determined by molecprob.
        All arguments are optional.
        Syntax: __init__(molecprob, reaction_rate,
        reaction_favoritism, assoc_stabilization,
        association_favoritism, dimensions)
        molecprob is an optional argument that must be between 0 and 1.
        It is the probability of a molecule starting in this lattice site.
        reaction_favoritism is the relative amount the
        positive product is favored in the reaction.
        assoc_stabilization is the amount association is favored
        over molecules by themselves.
        assoc_favoritism is the relative amount the
        positive product is favored in staying associated.
        excit_prob is the probability of an excitation at any step
        dimensions is the dimensionality of the lattice. Probably
        going to stay at 2d.
        Returns None."""
        # Check values of inputs
        if type(molecprob) != float:
            raise TypeError('molecprob must be a real number (float)')
        elif moleprob >= 1:
            raise ValueError('moleprob must be less than 1')
        elif moleprob <= 0:
            raise ValueError('moleprob must be greater than 0')
        elif reaction_rate <= 0:
            raise ValueError('reaction rate must be greater than 0')
        elif type(reaction_favoritism) != float:
            raise TypeError('reaction favoritism must be a float')
        elif assoc_stabilization <= 0:
            raise ValueError('assoc stabilization must be greater than 0')
        elif type(assoc_favoritism) != float:
            raise TypeError('assoc favoritism must be a float')
        elif (move_prob > 1
              or move_prob <= 0):
            raise ValueError('move prob must be between 0 and 1')
        # With Monte Carlo-like probability sets occupancy for catalyst.
        if random() > molecprob:
            self.catalyst = 1
        else:
            self.catalyst = 0
        # With Monte Carlo-like probability sets occupancy for htmf.
        # If occupied, can be in pm 1 orientation
        if random() > molecprob:
            self.htmf = pm()
        else:
            self.htmf = 0
        # With Monte Carlo-like probability sets occupancy for cinnamate.
        # If occupied, can be in pm 1 orientation
        if random() > molecprob:
            self.cinna = pm()
        else:
            self.cinna = 0
        # sets product occupancy to 0, as well as excitation state
        self.product = 0
        self.htmf_excitation_state = 0
        self.catalyst_excitation_state = 0
        # Need to define relevant variables for this!!
        # They should be taken as input to init and set as internal
        # variables here
        self.p_react_pos = reaction_rate * reaction_favoritism
        self.p_react_neg = reaction_rate / reaction_favoritism
        self.p_assoc     = assoc_stabilization
        self.p_assoc_pos = assoc_stabilization * assoc_favoritism
        self.p_assoc_neg = assoc_stabilization / assoc_favoritism
        self.excit_prob  = excit_prob
        self.move_prob   = move_prob
        self.max_move    = 2 * dimensions

    def __repr__(self):
        """Returns current occupation of this instance"""
        return self.catalyst, self.htmf, self.cinna, self.product

    def move(self):
        """If there is anything in this lattice site,
        then it will return which direction it wants to move as an integer
        and its current state.
        Syntax: move()
        Returns:
        (direction, catalyst, catalyst excitation, htmf, htmf
        excitation, cinna, product)"""
        # square each occupation (to avoid problems with negatives) and sum
        # this will be used to see if anything is here
        sum_state_squared = self.catalyst**2 + self.htmf**2 +
        self.cinna**2 + self.product**2
        if sum_state_squared > 0:
            move_direc = randint(1, self.max_move)
            return (move_direc, self.catalyst, self.catalyst_excitation_state,
                    self.htmf, self.htmf_excitation_state, self.cinna,
                    self.product)

    def move_result(self, rcatalyst, rcatalyst_excitation_state,
                    rhtmf, rhtmf_excitation_state, rcinna, rproduct)
        """This will be called after an attempted move to re-set
        the state of this instance based on what the other instance
        didn't accept.
        syntax: move_result(rcatalyst, rcatalyst_excitation_state,
        rhtmf, rhtmf_excitation_state, rcinna, rproduct)
        Returns None"""
        self.catalyst = rcatalyst
        self.htmf     = rhtmf
        self.cinna    = rcinna
        self.product  = rproduct
        self.catalyst_excited_state = self.catalyst_excited_state
        self.htmf_excited_state     = self.htmf_excited_state

    def excite(self):
        """This will check if there is anything to be excited,
        then it will try to excite it with probability excit_prob.
        syntax: excite()
        Returns None"""
        if self.catalyst != 0:
            if random < self.excit_prob:
                self.catalyst_excitation_state = 1
        if self.htmf != 0:
            if random < self.excit_prob:
                self.htmf_excitation_state = 1

    def react(self):
        """react takes no arguments. It will check the state of this instance,
        and if all three molecules are present, will react with probability
        given by p_react_pos or p_react_neg (depending on sign of
        the product of the molecule occupancies) times the amount of
        excitation.
        If it reacts, it will update occupancies
        Syntax: react()
        Returns: none"""
        occupancy_product = self.catalyst * self.htmf * self.cinna
        total_excit = max(self.htmf_excitation_state,
                          self.catalyst_excitation_state)
        if occupancy_product == -1:
            weighted_react_prob = total_excit * p_react_neg
            if random < weighted_react_prob:
                self.product = -1
                self.htmf    =  0
                self.cinna   =  0
        elif occupancy_product == 1:
            weighted_react_prob = total_excit * p_react_pos
            if random < weighted_react_prob:
                self.product = 1
                self.htmf    = 0
                self.cinna   = 0
        # Reduce the excitation state
        self.htmf_excitation_state /= 2
        self.catalyst_excitation_state /= 2

    def accept_move(self, pcatalyst, pcatalyst_excitation_state,
                    phtmf, phtmf_excitation_state, pcinna, pproduct):
        """This will be called when another instance tries to move here.
        Calculates probability of accepting with Monte Carlo probabilities.
        Inputs interpretted as the new proposed states. (phtmf, etc.)
        Sets returned values that will be returned at end. (rhtmf, etc.)
        syntax: accept_move(catalyst, catalyst_excitation_state,
        htmf, htmf_excitation_state, cinna, product)
        returns whatever it doesn't accept as:
        (catalyst, catalyst_excitation_state, htmf, htmf_excitation_state
        cinna, product)"""
        #
        #
        # No where here does it change the signs of the cinnamate or htmf
        # (doesn't switch orientation), but it wouldn't be hard to do
        # by just adding * pm() in several places
        #
        #
        # First things first: if there is a product already here,
        # this site will not accept anything:
        if self.product != 0;
            return (pcatalyst, pcatalyst_excitation_state,
                    phtmf, phtmf_excitation_state, pcinna, pproduct)
        # product of signed occupancies to determine which association
        # stabilization to use
        occupancy_p_product = phtmf * pcinna
        # sum of abs of occupancies to determine magnitude of stabilization
        occupancy_sum = abs(pcatalyst) + abs(phtmf) + abs(pcinna)
        if occupancy_p_product == -1:
            stabil = self.p_assoc_neg * occupancy_sum
        elif occupancy_p_product == 1:
            stabil = self.p_assoc_pos * occupancy_sum
        else:
            stabil = self.p_assoc     * occupancy_sum
        # check to see if it will break up based on stabil
        if random < stabil:
            # stay together
            if random < self.move_prob:
                # accepts move of everything if not occupied
                # now need to check occupancies
                # if empty, moves
                # We know self.product is empty from above check,
                # so we can accept that immediately
                self.product = pproduct
                rproduct     = 0
                if self.catalyst == 0:
                    self.catalyst = pcatalyst
                    self.catalyst_excitation_state =
                        pcatalyst_excitation_state
                    rcatalyst     = 0
                    rcatalyst_excitation_state = 0
                else:
                    rcatalyst = pcatalyst
                    rcatalyst_excitation_state = phtmf_excitation_state
                if self.htmf == 0:
                    self.htmf = phtmf
                    rhtmf     = 0
                    self.htmf_excitation_state = phtmf_excitation_state
                    rhtmf_excitation_state     = 0
                else:
                    rhtmf = phtmf
                    rhtmf_excitation_state = phtmf_excitation_state
                if self.cinna == 0:
                    self.cinna = pcinna
                    rcinna     = 0
                else:
                    rcinna = pcinna
            else:
                # still stay together, but don't accept move
                # set all return variables to proposed variables
                return (pcatalyst, pcatalyst_excitation_state,
                        phtmf, phtmf_excitation_state, pcinna, pproduct)
        else:
            # separate
            # will move each independently with own move prob
            if random < (self.move_prob - abs(self.catalyst)):
                # Moves catalyst with some probability give by move_prob.
                # If catalyst site is occupied, because move_prob < 1,
                # then it is not possible to move the catalyst
                self.catalyst = pcatalyst
                rcatalyst     = 0
                self.catalyst_excitation_state = pcatalyst_excitation_state
                rcatalyst     = 0
                rcatalyst_excitation_state = 0
            else:
                rcatalyst = pcatalyst
                rcatalyst_excitation_state = phtmf_excitation_state
            if random < (self.move_prob - abs(self.htmf)):
                self.htmf = phtmf
                rhtmf     = 0
                self.htmf_excitation_state = phtmf_excitation_state
                rhtmf_excitation_state     = 0
            else:
                rhtmf = phtmf
                rhtmf_excitation_state = phtmf_excitation_state
            if random < (self.move_prob - abs(self.cinna)):
                self.cinna = pcinna
                rcinna     = 0
            else:
                rcinna = pcinna
            if random < (self.move_prob - abs(self.product)):
            # self.product is 0, but this shouldn't affect anything
                self.product = pproduct
                rproduct     = 0
            else:
                rproduct = pproduct
        return (rcatalyst, rcatalyst_excitation_state,
                rhtmf, rhtmf_excitation_state, rcinna, rproduct)
        
