

def write_hilapw(skp, nq):
    #
    print(skp)
    # rx.in : Variation cell optimization
    #
    with open("sets.in", 'w') as f:
        print(" %d %d %d" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # rx.in : Variation cell optimization
    #
    with open("lapw.inSCF", 'w') as f:
        print(" %d %d %d" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # rx.in : Variation cell optimization
    #
    with open("lapw.inBAND", 'w') as f:
        print(" %d %d %d" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
    #
    # rx.in : Variation cell optimization
    #
    with open("lapw.inBAND", 'w') as f:
        print(" %d %d %d" % (nq[0]*2, nq[1]*2, nq[2]*2), file=f)
