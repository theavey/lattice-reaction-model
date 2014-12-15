import numpy as np

class Xfun:

     def __init__(self, first, second):
             self.ans = first + second
             self.nans = first - second
             print self.nans

     def __repr__(self):
         return ' '.join(map(str, (self.ans, self.nans)))
     # then use whatever.split(', ') to separate string
     # back into list

     def apply(self):
         print ('sum of sum and difference is %i' % (self.ans + self.nans))


def mapable_xfun(args):
    return Xfun(args[0], args[1])

testtup = (1,2)
testarray = np.empty((3,3), dtype = object)
# testarray[:] = testtup

for i in xrange(3):
    for j in xrange(3):
        testarray[i,j] = mapable_xfun(testtup)

vxfun = np.vectorize(mapable_xfun)


class Tester:

    def __init__(self, word):
        self.word = word

    def pw(self):
        print 'the word is', self.word

def iftester():
    if 1:
        pass
    elif 0:
        print 'problem'
    print 'good'

def rettest():
    return (1, 2, 3)

def fof3(a, b, c):
    print a, b, c
