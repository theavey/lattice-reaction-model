#! /usr/bin/env python

from random import random
# random() returns random number between 0 and 1
from random import choice
# choice(seq) psuedorandomly returns a value of in seq
from random import randint
# randint(a, b) psuedorandomly returns an integer in range (a, b) inclusive
from random import shuffle
# shuffle(list) will create a permutation of list in place with no return
from math import exp
# exponential function
from math import sqrt
# square root function


def pm():
    """this function returns + or - 1 psuedorandomly"""
    return choice((-1, 1))


class Lattice_Cell_Object:
    """This class contains a state that is its current state,
    and can process several messages.
    It is meant to simulate a lattice site in a lattice
    model of the reaction of 3-HTMF
    with methyl cinnamate catalyzed by TADDOL.
    init will set the state of the instance with some chance of
    having molecules present.
    checkstate will return the state of the instance.
    move will attempt to move anything to another randomly
    chosen adjacent instance.
    react will attempt to react if all molecules are present.
    acceptmove will see if this instance can accept a molecule
    or molecules.
    excite will possibly excite some molecules that are present."""

    def __init__(self, molecprob = 0.01, reaction_energy = 0.1,
                 reaction_favoritism = 2.0, assoc_stabilization = 0.1,
                 assoc_favoritism = 2.0, excit_prob = 0.01,
                 beta = 1.0,
                 dimension = 2):
        """This initializes an instance of lattice_cell_object with a
        chance of population
        of molecules determined by molecprob.
        All arguments are optional.
        Syntax: __init__(molecprob, reaction_energy,
        reaction_favoritism, assoc_stabilization,
        association_favoritism, excitation_probability, beta, dimension)
        molecprob is an optional argument that must be between 0 and 1.
        It is the probability of a molecule starting in this lattice site.
        reaction_energy is essentially the activation energy of the reaction.
        It is in units of k_B T.
        reaction_favoritism is the relative amount the
        positive product is favored in the reaction.
        assoc_stabilization is the amount association is favored
        over molecules by themselves in units of k_B T.
        assoc_favoritism is the relative amount the
        positive product is favored in staying associated.
        excit_prob is the probability of an excitation at any step.
        beta is the unitless inverse temperature.
        dimension is the dimensionality of the lattice. Probably
        going to stay at 2D.
        Returns None."""
        # Check values of inputs
        if type(molecprob) != float:
            raise TypeError('molecprob must be a real number (float), %s'
                            % type(molecprob))
        elif molecprob >= 1:
            raise ValueError('moleprob must be less than 1')
        elif molecprob <= 0:
            raise ValueError('moleprob must be greater than 0')
        elif reaction_energy <= 0:
            raise ValueError('reaction energy must be greater than 0')
        elif type(reaction_favoritism) != float:
            raise TypeError('reaction favoritism must be a float')
        elif assoc_stabilization <= 0:
            raise ValueError('assoc stabilization must be greater than 0')
        elif type(assoc_favoritism) != float:
            raise TypeError('assoc favoritism must be a float')
        elif beta <= 0.0:
            raise ValueError('beta (inv temp) must be greater than 0.0')
        # With molecprob probability sets occupancy for catalyst.
        if random() < molecprob:
            self.catalyst = 1
        else:
            self.catalyst = 0
        # With molecprob probability sets occupancy for htmf.
        # If occupied, can be in pm 1 orientation
        if random() < molecprob:
            self.htmf = pm()
        else:
            self.htmf = 0
        # With Monte Carlo-like probability sets occupancy for cinnamate.
        # If occupied, can be in pm 1 orientation
        if random() < molecprob:
            self.cinna = pm()
        else:
            self.cinna = 0
        # sets product occupancy to 0, as well as excitation state
        self.product = 0
        self.htmf_excitation_state = 0.0
        self.catalyst_excitation_state = 0.0
        # Need to define relevant variables for this!!
        # They should be taken as input to init and set as internal
        # variables here
        self.e_react_pos = reaction_energy / reaction_favoritism
        self.e_react_neg = reaction_energy
        self.e_assoc_pos = assoc_stabilization * assoc_favoritism
        self.e_assoc_neg = assoc_stabilization
        # This can support negative association favoritism now by just using
        # the arithemetic mean instead.
        if assoc_favoritism < 0.0:
            self.e_assoc = (self.e_assoc_neg + self.e_assoc_pos) / 2
        else:
            self.e_assoc = assoc_stabilization * sqrt(assoc_favoritism)
        self.excit_prob  = excit_prob
        self.max_move    = 2 * dimension
        self.beta        = beta
        # This number is used to determine the product repulsion in e_of_state
        # It depends on beta such that at high beta, it shouldn't overflow.
        self.large_num   = 700 / (2.5 * beta)
        #self.counter     = 1

    def __repr__(self):
        """Returns current occupation of this instance as a comma separated
        list."""
        return ','.join(map(str, (self.catalyst, self.htmf,
                                   self.cinna, self.product)))
        # !!! This needs to return a string.
        # !!! cannot return a list
        # !!! see testenv for possible implementation
        #
        # This has been fixed. Now returns in a way that's good
        # for Mathematica to interpret. Also, VERY IMPORTANT to
        # note that this cannot be any longer without changing
        # sampl_templ in initializelat, otherwise characters
        # will be dropped from the end sometimes.
        #

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
        sum_state_squared = self.catalyst**2 + self.htmf**2 + \
                                self.cinna**2 + self.product**2
        # Only if the total occupancy is non-zero will this be true.
        # If there is something here, it will attempt to move it
        # by returning a direction and the current state.
        # If it's empty, it will return None.
        if sum_state_squared:
            move_direc = randint(0, (self.max_move - 1))
            return (move_direc, self.catalyst, self.catalyst_excitation_state,
                    self.htmf, self.htmf_excitation_state, self.cinna,
                    self.product)
        else:
            return None

    def move_result(self, rcatalyst, rcatalyst_excitation_state,
                    rhtmf, rhtmf_excitation_state, rcinna, rproduct):
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
        self.catalyst_excitation_state = rcatalyst_excitation_state
        self.htmf_excitation_state     = rhtmf_excitation_state

    def excite(self):
        """This will check if there is anything to be excited,
        then it will try to excite it with probability excit_prob.
        syntax: excite()
        Returns None"""
        #randomnumber = random
        #print randomnumber, self.catalyst, self.excit_prob
        if self.catalyst != 0:
            if random() < self.excit_prob:
                self.catalyst_excitation_state = 1.0
        if self.htmf != 0:
            if random() < self.excit_prob:
                self.htmf_excitation_state = 1.0
        #print self.catalyst_excitation_state

    def react(self):
        """react takes no arguments. It will check the state of this instance,
        and if all three molecules are present, will react with probability
        given by p_react_pos or p_react_neg (depending on sign of
        the product of the molecule occupancies) times the amount of
        excitation.
        If it reacts, it will update occupancies
        Syntax: react()
        Returns: none"""
        #
        # If there is already a product here, don't want another reaction
        # occuring here and replacing it.
        #
        if self.product:
            return
        # This will be 0 if not fully occupied.
        occupancy_product = self.catalyst * self.htmf * self.cinna
        # This assumes only one thing needs to be excited and that both the
        # HTMF and TADDOL catalyst can be excited
        total_excit = max(self.htmf_excitation_state,
                          self.catalyst_excitation_state)
        # If not all occupied, skip this, move to reducing excitation
        # state.
        # This will produce the negative product with some prob.
        # If it does react, it will create product,
        # and remove both reactants and all excitation.
        if occupancy_product == -1:
            reaction_energy = self.e_react_neg
            # with psuedo-MMC probability, this will react. The probability of
            # the reaction is attenuated by the excitation of the molecules.
            if random() < total_excit * exp(-self.beta * reaction_energy):
                # Create product.
                self.product = -1
                # Remove reactants and excitation.
                self.htmf    =  0
                self.htmf_excitation_state     = 0
                self.cinna   =  0
                self.catalyst_excitation_state = 0
        # Does the same except for the positive product.
        elif occupancy_product == 1:
            reaction_energy = self.e_react_pos
            if random() < total_excit * exp(-self.beta * reaction_energy):
                self.product = 1
                self.htmf    = 0
                self.htmf_excitation_state     = 0
                self.cinna   =  0
                self.catalyst_excitation_state = 0
        # Reduce the excitation state (need to be floats, not integers).
        self.htmf_excitation_state /= 2.
        self.catalyst_excitation_state /= 2.

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
        # Define current and proposed states as vectors (cleaner, I think)
        p_state = [pcatalyst, pcatalyst_excitation_state, phtmf,
                   phtmf_excitation_state, pcinna, pproduct]
        c_state = [self.catalyst, self.catalyst_excitation_state,
                   self.htmf, self.htmf_excitation_state, self.cinna,
                   self.product]
        # List of the locations of the molecules in the state vectors
        #
        # If changing the state to a vector as opposed to separate
        # variables, this will obviously need to be changed here.
        #
        molecule_list = [0, 2, 4]
        # randomize this list to check each sequentially, but in a
        # random order (shuffle does this in place, aka changes variable
        # definition)
        shuffle(molecule_list)
        # proposeto will be a vector that will be sent to check_proposed
        # after certain parts are replaced.
        # proposeto is the site being moved into, and proposefrom is being
        # moved out of.
        #
        # The colon is absolutely necessary here, otherwise c_state will
        # change every time propose is changed. I imagine this is slower
        # unfortunately, because it probably has to go item-wise.
        #
        proposeto   = c_state[:]
        proposefrom = p_state[:]
        # The [5] is for the product, which I think I want moved first.
        # Might get rid of this later if I can get it working.
        for mol in [5] + molecule_list:
            # Try to switch each molecule between the two states
            proposeto[mol]   = p_state[mol]
            proposefrom[mol] = c_state[mol]
            # Check proposal. If it returns true, move was accepted,
            # otherwise it was rejected. Set states accordingly.
            if self.check_proposed(proposeto, proposefrom, c_state, p_state):
                c_state[mol] = proposeto[mol]
                p_state[mol] = proposefrom[mol]
                # If it is an "excitable" molecule, need to include
                # excitation state as well. Will check that here.
                # The only excitable molecules as implemented are the HTMF
                # and TADDOL (catalyst) at locations 0 and 2 in the vector.
                #
                # If changing the state to a vector as opposed to separate
                # variables, this will obviously need to be changed here.
                #
                if mol in (0, 2):
                    proposeto[mol + 1]   = p_state[mol + 1]
                    proposefrom[mol + 1] = c_state[mol + 1]
                    c_state[mol + 1] = proposeto[mol + 1]
                    p_state[mol + 1] = proposefrom[mol + 1]
            else:
                # reset proposes for next time around the loop if not accepted
                proposeto[mol]   = c_state[mol]
                proposefrom[mol] = p_state[mol]
        # End of accept move sequence:
        # Sets the state of this site to be whatever was decided above
        # and consequently set as the current state.
        #
        # If changing the state to a vector as opposed to separate
        # variables, this will obviously need to be changed here.
        #
        self.catalyst                  = c_state[0]
        self.catalyst_excitation_state = c_state[1]
        self.htmf                      = c_state[2]
        self.htmf_excitation_state     = c_state[3]
        self.cinna                     = c_state[4]
        self.product                   = c_state[5]
        # Then, it returns whatever was not accepted as a vector
        return p_state

    def e_of_state(self, state):
        """Calculates and returns the energy of the input
        state * beta (inverse temperature). Input should be a
        list-like object of length 6 of the form:
        catalyst, catalyst ex. state, htmf, htmf ex. state,
        cinnamate, product.
        Returns beta * E (a float)."""
        # Interesting extension note: this does not take into account the
        # excitation state of the molecules, which really might be related
        # to the association energy.
        #
        # A few terms to help make things cleaner later:
        # Product of occupancies to determine which association stabilization
        # to use:
        occ_prod = state[2] * state[4]
        # Sum of abs of occupancies to determine magnitude of stabilization
        # (or destabilization if there is a product present)
        occ_sum = abs(state[0]) + abs(state[2]) + abs(state[4])
        # Term to help ensure nothing overlaps with a product.
        # It will be 0 if product is 0, large negative if it has product and
        # nothing else, otherwise it will be quite large and positive.
        product_repulsion = abs(state[5]) * (occ_sum - 0.5) * self.large_num
        # Stabilization due to occupancy. This depends on the sign of
        # occ_prod and the magnitude of occ_sum.
        if occ_prod == -1:
            stabil = self.e_assoc_neg * occ_sum
        elif occ_prod == 1:
            stabil = self.e_assoc_pos * occ_sum
        else:
            stabil = self.e_assoc     * occ_sum
        # Finally, calculate the energy. The stabilization
        # will lower the energy (making it more stable).
        energy = product_repulsion - stabil
        return self.beta * energy

    def check_proposed(self, proposedto, proposedfrom, currentfrom,
                       currentto):
        """Takes two possible states and performs a metropolis MC
        comparison. It depends on e_of_state for energy calculations.
        takes two states as vectors and returns two vectors.
        Usage check_proposed(proposed, current)
        Returns 1 if accepted, 0 if not accepted."""
        # Calculate energies of current and proposed states
        #
        # These energies include the temperature.
        #
        # p_energy is energy of proposed states after move from 'from' to 'to'
        # state. c_energy is the energy of both states before any move.
        p_energy = self.e_of_state(proposedto) + self.e_of_state(proposedfrom)
        c_energy = self.e_of_state(currentto)  + self.e_of_state(currentfrom)
        ## if exp(-(p_energy - c_energy)) != 1.0:
        ##     print ""
        ##     print "look, it's not 1!!"
        ##     print p_energy, " , ", c_energy
        ##     print exp(-(p_energy - c_energy))
        if random() < exp(-(p_energy - c_energy)):
            # Negative of energy of proposed state minus the energy
            # of the current state. If it is going down in energy (or
            # the energy of the proposed state is the same), this
            # will always be true because the random number cannot be
            # > 1. Otherwise, there is some random probability that it
            # may happen depending on the magnitude of the energy
            # difference.
            # Note, the temperature is included in this energy term.
            return 1
        else:
            # move not accepted:
            ## print "move not accepted"
            return 0
