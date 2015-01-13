#! /usr/bin/env python

import initializelat as initlat
from time import strftime
from random import randint
nowtime = strftime("%Y%m%d%H%M")
# Output file is named lat_react_out followed by the time, then a random
# three digit integer to help ensure no name overlap if running multiple at
# the same time. BTW, str.zfill makes the string at least 3 characters long
# by left filling with zeros.
output_file_name = ("lat_react_out" + nowtime + str(randint(1,999)).zfill(3))

def run_react(arguments = (0.1, 0.1, 2.0, 0.5, 2.0, 0.5, 1.0),
                           steps = 2000, size = 10, sample_int = 20,
                           dimension = 2):
    """This function will run a reaction based on the lattice implementation
    defined in initializelat (which depends on laticecellobject).
    This will take arguments of parameters for each cell, number of steps,
    length along a side of the lattice, sampling interval, and dimensionality.
    The parameters for each cells are
    molecprob, reaction_energy,
    reaction_favoritism, assoc_stabilization,
    assoc_favoritism, excit_prob, beta
    All arguments are optional. It will return an array of arrays of strings
    that are the occupancies of the sites"""
    # Initialize lattice:
    lattice = initlat.Lattice(arguments, size, dimension)
    # Header for the output file. First human readable format then
    # in machine readable format.
    header  = '{{"size:%(size)i","steps:%(steps)i","sample_int:%(samp_int)i"\
            ,"args:%(argus)s"},\n{%(size)i,%(steps)i,%(samp_int)i}\
            \n' % {'size': size, 'steps': steps, \
                    'samp_int': sample_int, 'argus': arguments}
    # Open output file. Using "with" will automatically close the file
    # if an exception occurs and it has to quit.
    with open(output_file_name, 'w') as out_file:
        out_file.write(header)
        for step in xrange(steps + 1):
            lattice.over_sites('excite')
            lattice.over_sites('react')
            lattice.over_sites('move')
            # sample every sample_int
            # goes to steps + 1 so that it samples the last run
            if step in xrange(0, steps + 1, sample_int):
                # This is is a numpy array
                a_sample = lattice.over_sites('sample')
                # This comma separates the previous entry and then starts the
                # next. This formatting is done so that it can be easily read
                # in Mathematica, but the braces are matched so other programs
                # should be able to read it as well.
                out_file.write(',{{')
                # Because a_sample is an np array, the tofile will output the
                # file to out_file separated by sep
                a_sample.tofile(out_file, sep = '},{')
                out_file.write('}}')
            # Print the current number of steps for every tenth, just to keep
            # track and give an estimate of how long it may take to finish.
            if step in xrange(0, steps, steps / 10):
                print step
        # Final brace matches the open brace remaining.
        out_file.write('}')
    print output_file_name
        
    


# I'm not really sure how this works (don't remember what __name__ or
# __main__ are, but I took this from 1Dautomaton I believe that I
# wrote for that CS 11 class.
if __name__ == "__main__":
    import sys
    import os
    prog_name = os.path.basename(sys.argv[0])
    usage = "usage: %s [steps, size, " % prog_name + \
      "sample_int, dimension, arguments]"
    #try:
    if   len(sys.argv) == 1:
        # No arguments given
        print "No arguments provided, running with default values"
        run_react()
    elif len(sys.argv) == 8:
        print "Running with the one set of arguments given"
        argums = sys.argv[1:]
        run_react(arumgs[0:])
    elif len(sys.argv) == 11:
        print "Running with the 4 arguments given"
        argums = sys.argv[1:]
        run_react(argums[3:], int(argums[0]), int(argums[1]),
                  int(argums[2]))
    elif len(sys.argv) == 12:
        print "Running with the 5 arguments given"
        argums = sys.argv[1:]
        run_react(argums[4:],  int(argums[0]), int(argums[1]),
                  int(argums[2]), int(argums[3]))
    else:
        print "Weird number of args given (%i), " % len(sys.argv)
        print "running with default values. Usage: "
        print usage
        run_react()
    #except(TypeError), args:
    #    print usage
    #    print args
      
