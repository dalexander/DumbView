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

def toString(lst):
    return "".join(lst)

def rewrite(t, q, x):
    assert len(t) == len(q) == len(x) == 2
    t0, t1 = t
    q0, q1 = q
    x0, x1 = x
    if (t0 == t1 == q0) and (q1 == "-"):
        return (t, q1 + q0, x1 + x0)
    elif (q0 == q1 == t0) and (t1 == "-"):
        return (t1 + t0, q, x1 + x0)
    else:
        return (t, q, x)

def realign(t, q, x):
    assert len(t) == len(q) == len(x)
    t = list(t)
    q = list(q)
    x = list(x)

    def rewriteScan(t, q, x):
        origX = x[:]
        for i in xrange(len(t)-2,-1,-1):
            t[i:i+2], q[i:i+2], x[i:i+2] = rewrite(t[i:i+2], q[i:i+2], x[i:i+2])
        changed = (x != origX)
        return changed

    goAgain = True
    while goAgain:
        goAgain = rewriteScan(t, q, x)

    return (toString(t), toString(q), toString(x))



def test_realignment():
    #from nose import assert_equals as EQ

    assert realign("AAAAAA",
                   "AAA-AA",
                   "MMMDMM") == \
                  ('AAAAAA',
                   '-AAAAA',
                   'DMMMMM')

    assert realign("A-AAAA",
                   "AAAAAA",
                   "MIMMMM") == \
                  ('-AAAAA',
                   'AAAAAA',
                   'IMMMMM')

    # Doesn't realign over base not matching the homopolymer
    assert realign("GATTTACA",
                   "GAG-TACA",
                   "MMRIMMMM") == \
                  ('GATTTACA',
                   'GAG-TACA',
                   'MMRIMMMM')

    assert realign("AAAAAA",
                   "AAA--A",
                   "MMMIIM") == \
                  ('AAAAAA',
                   '--AAAA',
                   'IIMMMM')

    assert realign("AAAAAA",
                   "A-AA-A",
                   "MIMMIM") == \
                  ('AAAAAA',
                   '--AAAA',
                   'IIMMMM')



if __name__ == '__main__':
    print "Running tests"
    test_realignment()
