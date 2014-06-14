#
# Push gaps forward in homopolymers
#
# Rewrite rule 1: XX    ===>  XX
#                 X-          -X
#
# Rewrite rule 2: X-    ===>  -X
#                 XX          XX
#
# Iterate until convergence.
#
# Also could consider alternative, more aggressive (rewritePrime)
#
# Rewrite rule 1': XX    ===>  XX
#                  Y-          -Y
#
# Rewrite rule 2': Y-    ===>  -Y
#                  XX          XX

def toString(lst):
    return "".join(lst)

def rewrite(t, q):
    assert len(t) == len(q) == 2
    t0, t1 = t
    q0, q1 = q
    if (t0 == t1 == q0) and (q1 == "-"):
        return (t, q1 + q0)
    elif (q0 == q1 == t0) and (t1 == "-"):
        return (t1 + t0, q)
    else:
        return (t, q)

def rewritePrime(t, q):
    # Allow realign around base mismatch in t/q
    assert len(t) == len(q) == 2
    t0, t1 = t
    q0, q1 = q
    if (t0 == t1) and (q1 == "-"):
        return (t, q1 + q0)
    elif (q0 == q1) and (t1 == "-"):
        return (t1 + t0, q)
    else:
        return (t, q)


def realign(t, q, rewrite=rewrite):
    assert len(t) == len(q)
    t = list(t)
    q = list(q)

    def rewriteScan(t, q):
        for i in xrange(len(t)-2):
            t[i:i+2], q[i:i+2] = rewrite(t[i:i+2], q[i:i+2])

    oldT, oldQ = None, None
    while (oldT, oldQ) != (t, q):
        oldT, oldQ = t[:], q[:]
        #print toString(t), toString(q)
        rewriteScan(t, q)

    return (toString(t), toString(q))


def realignPrime(t, q):
    return realign(t, q, rewrite=rewritePrime)


def test_realignment():
    #from nose import assert_equals as EQ

    assert realign("AAAAAA",
                   "AAA-AA") == \
                  ('AAAAAA',
                   '-AAAAA')

    assert realign("A-AAAA",
                   "AAAAAA") == \
                  ('-AAAAA',
                   'AAAAAA')

    # Doesn't realign over base not matching the homopolymer, unless we
    # use the prime rules

    assert realign("GATTTACA",
                   "GAG-TACA") == \
                  ('GATTTACA',
                   'GAG-TACA')

    assert realignPrime("GATTTACA",
                        "GAG-TACA") == \
                       ('GATTTACA',
                        'GA-GTACA')


if __name__ == '__main__':
    print "Running tests"
    test_realignment()
