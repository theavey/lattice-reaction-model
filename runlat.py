#! /usr/bin/env python

import initializelat as initlat
import time
nowtime = time.strftime("%Y%m%d%H%M")
output_file_name = ("lat_react_out" + nowtime)

def run_react(arguments = (0.5, 2.0, 2.0, 2., 2.0, 0.9, 0.75),
                           steps = 20000, size = 40, sample_int = 1000,
                           dimension = 2):
    """This function will run a reaction based on the lattice implementation
    defined in initializelat (which depends on laticecellobject).
    This will take arguments of parameters for each cell, number of steps,
    length along a side of the lattice, sampling interval, and dimensionality.
    The parameters for each cells are
    molecprob, reaction_rate,
    reaction_favoritism, assoc_stabilization,
    assoc_favoritism, excit_prob,
    move_prob
    All arguments are optional. It will return an array of arrays of strings
    that are the occupancies of the sites"""
    # Initialize lattice:
    lattice = initlat.Lattice(arguments, size, dimension)
    header  = '{{"size:%(size)i","steps:%(steps)i","sample_int:%(samp_int)i"\
            ,"args:%(argus)s"},\n{%(size)i,%(steps)i,%(samp_int)i}\
            \n' % {'size': size, 'steps': steps, \
                    'samp_int': sample_int, 'argus': arguments}
    # Open output file
    with open(output_file_name, 'w') as out_file:
        out_file.write(header)
        for step in xrange(steps + 1):
            lattice.over_sites('excite')
            lattice.over_sites('react')
            lattice.over_sites('move')
            if step in xrange(0, steps + 1, sample_int):
                a_sample = lattice.over_sites('sample')
                out_file.write(',{{')
                a_sample.tofile(out_file, sep = '},{')
                out_file.write('}}')
            if step in xrange(0, steps, steps / 10):
                print step
        out_file.write('}')
        
    



if __name__ == "__main__":
    import sys
    import os
    prog_name = os.path.basename(sys.argv[0])
    usage = "usage: %s [steps, size, " % prog_name + \
      "sample_int, dimension, arguments]"
    try:
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
    except(TypeError), args:
        print usage
        print args
      
