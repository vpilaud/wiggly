from collections import defaultdict
from functools import cmp_to_key

### WIGGLY COMPLEX ON ARCS

@cached_function
def side_points(P, i, j, k):
    r"""
    Return the relative position of P[i], P[j], P[k].
    This is 1, -1, or 0.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for (i,j,k) in Subsets(range(5), 3):
        ....:     print(i, j, k, side_points(P, i, j, k))
        ....:
        0 1 2 -1
        0 1 3 1
        0 1 4 0
        0 2 3 1
        0 2 4 1
        0 3 4 -1
        1 2 3 1
        1 2 4 1
        1 3 4 -1
        2 3 4 -1
    """
    return sign(Matrix([[P[i][0], P[i][1], 1], [P[j][0], P[j][1], 1], [P[k][0], P[k][1], 1]]).determinant())

@cached_function
def position(P, i, j, k):
    r"""
    Return the relative position of P[i], P[j], P[k].
    This is 1 or -1 if the points are not aligned, and True or False if the points are aligned (False means that k lies in the interior of the segment joining i to j).

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for (i,j,k) in Subsets(range(5), 3):
        ....:     print(i, j, k, position(P, i, j, k))
        ....:
        0 1 2 -1
        0 1 3 1
        0 1 4 True
        0 2 3 1
        0 2 4 1
        0 3 4 -1
        1 2 3 1
        1 2 4 1
        1 3 4 -1
        2 3 4 -1
    """
    d = Matrix([[P[i][0], P[i][1], 1], [P[j][0], P[j][1], 1], [P[k][0], P[k][1], 1]]).determinant()
    if d != 0: return sign(d)
    return i == k or j == k or any((P[i][c] - P[k][c])*(P[j][c] - P[k][c]) > 0 for c in [0,1])

@cached_function
def intermediate_points(P, i ,j):
    r"""
    Return the indices k such that P[k] lies in the open segment joining P[i] to P[j].

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for (i,j) in Subsets(range(5), 2):
        ....:     print(i, j, intermediate_points(P, i, j))
        ....:
        0 1 {}
        0 2 {}
        0 3 {}
        0 4 {1}
        1 2 {}
        1 3 {}
        1 4 {}
        2 3 {}
        2 4 {}
        3 4 {}
    """
    return Set([k for k in range(len(P)) if not position(P, i, j, k)])

@cached_function
def is_coherently_labeled(P):
    r"""
    Tests if the point set P is coherently labeled.
    This means that aligned points are labeled in increasing/decreasing order along each line.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: is_coherently_labeled(P)
        True

        sage: P = ((0,0), (2,0), (0,2), (2,2), (1,1))
        sage: is_coherently_labeled(P)
        False
    """
    return all(position(P, i, j, k) for (i,j,k) in [sorted(s) for s in Subsets(range(len(P)), 3)])

def arc(i, j, A, B):
    r"""
    Creates an arc.

    EXAMPLES::
        sage: arc(3, 6, [2,5], [4])
        (3, 6, {2, 5}, {4})
    """
    return (i, j, Set(A), Set(B))

def reverse_arc(a):
    r"""
    Reverses an arc a.

    EXAMPLES::
        sage: reverse_arc(arc(3, 6, [2,5], [4]))
        (6, 3, {4}, {2, 5})
    """
    return arc(a[1], a[0], a[3], a[2])

@cached_function
def side_arc_point(P, a, k):
    r"""
    Given an arc a and a point k, tell on which side of the arc a is k.
    This is 1, -1 or 0.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for k in range(5):
        ....:     print(k, side_arc_point(P, arc(0, 4, [1], []), k))
        ....:
        0 0
        1 1
        2 -1
        3 1
        4 0
    """
    (i, j, A, B) = a
    if position(P, i, j, k) is False: return 1 if k in A else -1
    if position(P, i, j, k) is True: return 0
    return position(P, i, j, k)

@cached_function
def is_boundary_arc(P, a):
    r"""
    Test if the arc a is a boundary arc of P.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: is_boundary_arc(P, arc(0, 2, [], []))
        True
        sage: is_boundary_arc(P, arc(2, 3, [], []))
        False
        sage: is_boundary_arc(P, arc(0, 4, [1], []))
        False
        sage: P = ((0,0), (1,0), (3,0), (0,3), (3,3))
        sage: is_boundary_arc(P, arc(0, 2, [1], []))
        True
        sage: is_boundary_arc(P, arc(0, 2, [], [1]))
        False
    """
    sides = [side_arc_point(P, a, k) for k in range(len(P)) if k != a[0] and k != a[1]]
    return all(side != 0 for side in sides) and len(Set(sides)) == 1

@cached_function
def arcs(P):
    r"""
    Return all arcs on a point set P.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: arcs(P)
        [(0, 1, {}, {}),
         (0, 2, {}, {}),
         (1, 2, {}, {}),
         (0, 3, {}, {}),
         (1, 3, {}, {}),
         (2, 3, {}, {}),
         (0, 4, {}, {1}),
         (0, 4, {1}, {}),
         (1, 4, {}, {}),
         (2, 4, {}, {}),
         (3, 4, {}, {})]
    """
    return [(i, j, A, intermediate_points(P, i ,j).difference(A)) for j in range(len(P)) for i in range(j) for A in Subsets(intermediate_points(P, i ,j))]

@cached_function
def relevant_arcs(P):
    r"""
    Return all relevant arcs on a point set P.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: relevant_arcs(P)
        [(0, 1, {}, {}),
         (1, 2, {}, {}),
         (1, 3, {}, {}),
         (2, 3, {}, {}),
         (0, 4, {}, {1}),
         (0, 4, {1}, {}),
         (1, 4, {}, {})]
    """
    return [a for a in arcs(P) if not is_boundary_arc(P, a)]

@cached_function
def boundary_arcs(P):
    r"""
    Return all boundary arcs on a point set P.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: boundary_arcs(P)
        [(0, 2, {}, {}), (0, 3, {}, {}), (2, 4, {}, {}), (3, 4, {}, {})]
    """
    return [a for a in arcs(P) if is_boundary_arc(P, a)]

def are_crossing(P, a1, a2):
    r"""
    Test if the two arcs a1 and a2 are crossing.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for (a1, a2) in Subsets(arcs(P), 2):
        ....:     if are_crossing(P, a1, a2):
        ....:         print(a1, a2)
        ....:
        (0, 4, {1}, {}) (1, 2, {}, {})
        (0, 4, {}, {1}) (1, 3, {}, {})
        (0, 4, {}, {1}) (2, 3, {}, {})
        (0, 4, {1}, {}) (2, 3, {}, {})
        (1, 4, {}, {}) (2, 3, {}, {})
    """
    assert(is_coherently_labeled(P), 'coherently relabel your point set')
    (i1, j1, A1, B1) = a1
    (i2, j2, A2, B2) = a2
    if A1.intersection(B2).union(A1.intersection(Set([i2,j2]))).union(B2.intersection(Set([i1,j1]))) and B1.intersection(A2).union(B1.intersection(Set([i2,j2]))).union(A2.intersection(Set([i1,j1]))): return True
    return side_arc_point(P, a1, i2) * side_arc_point(P, a1, j2) < 0 and side_arc_point(P, a2, i1) * side_arc_point(P, a2, j1) < 0

def is_crossing(P, b, aa):
    r"""
    Test if the arc b is crossing any of the arcs in aa.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: b = arc(0, 4, [1], [])
        sage: is_crossing(P, b, [arc(0, 1, [], []), arc(1, 3, [], [])])
        False
        sage: is_crossing(P, b, [arc(0, 1, [], []), arc(1, 2, [], [])])
        True
    """
    return any(are_crossing(P, a, b) for a in aa)

def is_pointed(P, b, aa):
    r"""
    Test if the arc b is pointed with the set of arcs aa.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: b = arc(0, 1, [], [])
        sage: is_pointed(P, b, [arc(1, 2, [], []), arc(1, 3, [], [])])
        False
        sage: is_pointed(P, b, [arc(1, 2, [], []), arc(1, 4, [], [])])
        False
        sage: is_pointed(P, b, [arc(1, 2, [], []), arc(0, 4, [], [])])
        True
    """
    both_direction_aa = list(aa) + [reverse_arc(a) for a in aa]
    for i in [0,1]:
        aai = [a for a in both_direction_aa if a[0] == b[i]]
        for a in aai:
            if not position(P, a[1], b[1-i], b[i]): return False
        for (a1, a2) in Subsets(aai, 2):
            if side_points(P, b[i], a1[1], a2[1]) * side_points(P, b[i], a1[1], b[1-i]) < 0 and side_points(P, b[i], a2[1], a1[1]) * side_points(P, b[i], a2[1], b[1-i]) < 0 and side_points(P, b[i], b[1-i], a1[1]) * side_points(P, b[i], b[1-i], a2[1]) < 0: return False
    return True

def are_pointed(P, aa):
    r"""
    Test if the set aa of arcs is pointed.

    EXAMPLES::
        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for aa in Subsets(arcs(P), 3):
        ....:     if not are_pointed(P, aa):
        ....:         print(aa)
        ....:
        {(1, 4, {}, {}), (0, 1, {}, {}), (0, 2, {}, {})}
        {(1, 3, {}, {}), (0, 1, {}, {}), (1, 2, {}, {})}
        {(1, 4, {}, {}), (0, 1, {}, {}), (1, 2, {}, {})}
        {(0, 3, {}, {}), (1, 4, {}, {}), (0, 1, {}, {})}
        {(1, 3, {}, {}), (0, 1, {}, {}), (1, 4, {}, {})}
        {(1, 4, {}, {}), (0, 1, {}, {}), (2, 3, {}, {})}
        {(1, 4, {}, {}), (0, 4, {}, {1}), (0, 1, {}, {})}
        {(1, 4, {}, {}), (0, 4, {1}, {}), (0, 1, {}, {})}
        {(1, 4, {}, {}), (0, 1, {}, {}), (2, 4, {}, {})}
        {(1, 4, {}, {}), (3, 4, {}, {}), (0, 1, {}, {})}
    """
    both_direction_arcs = list(aa) + [reverse_arc(a) for a in aa]
    for i in range(len(P)):
        aai = [a for a in both_direction_arcs if a[0] == i]
        for (a1, a2) in Subsets(aai, 2):
            if not position(P, a1[1], a2[1], i): return False
        for (a1, a2, a3) in Subsets(aai, 3):
            if side_points(P, i, a1[1], a2[1]) * side_points(P, i, a1[1], a3[1]) < 0 and side_points(P, i, a2[1], a1[1]) * side_points(P, i, a2[1], a3[1]) < 0 and side_points(P, i, a3[1], a1[1]) * side_points(P, i, a3[1], a2[1]) < 0: return False
    return True

@cached_function
def wiggly_complex_enum(P):
    r"""
    Enumerates the wiggly complex on the point set P.
    """
    seeds = [([], relevant_arcs(P))]
    succ = lambda aabb: [(aabb[0]+[b], aabb[1][i+1:]) for (i,b) in enumerate(aabb[1]) if is_pointed(P, b, aabb[0]) and not is_crossing(P, b, aabb[0])]
    return RecursivelyEnumeratedSet(seeds, succ, structure='forest')

@cached_function
def wiggly_complex(P):
    r"""
    Return the wiggly complex on the point set P.

    EXAMPLES::
        sage: P = ((0,0), (1,0), (2,0))
        sage: WC = wiggly_complex(P)
        sage: WC.f_vector()
        [1, 2]
        sage: WC.is_pure()
        True
        sage: WC.is_pseudomanifold()
        True
        sage: add((-1)^i * WC.f_vector()[i] for i in range(len(WC.f_vector())))
        -1

        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: WC = wiggly_complex(P)
        sage: WC.is_pure()
        True
        sage: WC.is_pseudomanifold()
        True
        sage: WC.f_vector()
        [1, 7, 15, 10]
        sage: add((-1)^i * WC.f_vector()[i] for i in range(len(WC.f_vector())))
        -1

        sage: P = ((0,0), (1,1), (2,0), (0,2), (2,2))
        sage: WC = wiggly_complex(P)
        sage: WC.is_pure()
        True
        sage: WC.is_pseudomanifold()
        True
        sage: WC.f_vector()
        [1, 8, 18, 12]
        sage: add((-1)^i * WC.f_vector()[i] for i in range(len(WC.f_vector())))
        -1

        sage: P = ((0,0), (6,0), (2,2), (4,2), (3,3), (0,6), (6,6))
        sage: WC = wiggly_complex(P)
        sage: WC.is_pure()
        True
        sage: WC.is_pseudomanifold()
        True
        sage: WC.f_vector()
        [1, 27, 242, 1051, 2500, 3331, 2331, 666]
        sage: add((-1)^i * WC.f_vector()[i] for i in range(len(WC.f_vector())))
        -1

        sage: P = ((0,0), (6,0), (2,2), (3,3), (4,4), (0,6), (6,6))
        sage: WC = wiggly_complex(P)
        sage: WC.is_pure()
        True
        sage: WC.is_pseudomanifold()
        True
        sage: WC.f_vector()
        [1, 34, 350, 1664, 4190, 5764, 4088, 1168]
        sage: add((-1)^i * WC.f_vector()[i] for i in range(len(WC.f_vector())))
        -1
    """
    return SimplicialComplex([x[0] for x in wiggly_complex_enum(P)])

@cached_function
def rigidity_vector(P, a):
    r"""
    Return the rigidity vector of the arc a.

    EXAMPLES::
        sage: P = ((0,0), (1,0), (2,0))
        sage: for a in arcs(P):
        ....:     print(a, rigidity_vector(P, a))
        (0, 1, {}, {}) (1, 0, -1, 0, 0, 0)
        (0, 2, {}, {1}) (2, 1, 0, -2, -2, 1)
        (0, 2, {1}, {}) (2, -1, 0, 2, -2, -1)
        (1, 2, {}, {}) (0, 0, 1, 0, -1, 0)

        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for a in arcs(P):
        ....:     print(a, rigidity_vector(P, a))
        ....:
        (0, 1, {}, {}) (1, 1, -1, -1, 0, 0, 0, 0, 0, 0)
        (0, 2, {}, {}) (3, 0, 0, 0, -3, 0, 0, 0, 0, 0)
        (1, 2, {}, {}) (0, 0, 2, -1, -2, 1, 0, 0, 0, 0)
        (0, 3, {}, {}) (0, 3, 0, 0, 0, 0, 0, -3, 0, 0)
        (1, 3, {}, {}) (0, 0, -1, 2, 0, 0, 1, -2, 0, 0)
        (2, 3, {}, {}) (0, 0, 0, 0, -3, 3, 3, -3, 0, 0)
        (0, 4, {}, {1}) (1, 5, 3, -3, 0, 0, 0, 0, -4, -2)
        (0, 4, {1}, {}) (5, 1, -3, 3, 0, 0, 0, 0, -2, -4)
        (1, 4, {}, {}) (0, 0, 2, 2, 0, 0, 0, 0, -2, -2)
        (2, 4, {}, {}) (0, 0, 0, 0, 0, 3, 0, 0, 0, -3)
        (3, 4, {}, {}) (0, 0, 0, 0, 0, 0, 3, 0, -3, 0)

        sage: P = ((0,0), (1,1), (2,0), (0,2), (2,2))
        sage: for a in arcs(P):
        ....:     print(a, rigidity_vector(P, a))
        ....:
        (0, 1, {}, {}) (1, 1, -1, -1, 0, 0, 0, 0, 0, 0)
        (0, 2, {}, {}) (2, 0, 0, 0, -2, 0, 0, 0, 0, 0)
        (1, 2, {}, {}) (0, 0, 1, -1, -1, 1, 0, 0, 0, 0)
        (0, 3, {}, {}) (0, 2, 0, 0, 0, 0, 0, -2, 0, 0)
        (1, 3, {}, {}) (0, 0, -1, 1, 0, 0, 1, -1, 0, 0)
        (2, 3, {}, {1}) (0, 0, 2, 2, -3, 1, 1, -3, 0, 0)
        (2, 3, {1}, {}) (0, 0, -2, -2, -1, 3, 3, -1, 0, 0)
        (0, 4, {}, {1}) (1, 3, 2, -2, 0, 0, 0, 0, -3, -1)
        (0, 4, {1}, {}) (3, 1, -2, 2, 0, 0, 0, 0, -1, -3)
        (1, 4, {}, {}) (0, 0, 1, 1, 0, 0, 0, 0, -1, -1)
        (2, 4, {}, {}) (0, 0, 0, 0, 0, 2, 0, 0, 0, -2)
        (3, 4, {}, {}) (0, 0, 0, 0, 0, 0, 2, 0, -2, 0)
    """
    (i, j, A, B) = a
    v = vector(P[i]) - vector(P[j])
    d = 0 if v[0] != 0 else 1
    w = vector([-v[1], v[0]])
    res = []
    for k in range(len(P)):
        if k == i: c = (-v + add((P[j][d]-P[l][d])/(P[j][d]-P[i][d])*w for l in A) - add((P[j][d]-P[l][d])/(P[j][d]-P[i][d])*w for l in B))[:2]
        elif k == j: c = (v + add((P[l][d]-P[i][d])/(P[j][d]-P[i][d])*w for l in A) - add((P[l][d]-P[i][d])/(P[j][d]-P[i][d])*w for l in B))[:2]
        elif k in A: c = -w[:2]
        elif k in B: c = w[:2]
        else: c = (0,0)
        res = res + list(c)
    return vector(res)

@cached_function
def rigidity_fan(P):
    r"""
    Return the rigidity fan of P.

    EXAMPLES::
        sage: P = ((0,0), (1,0), (2,0))
        sage: F
        Rational polyhedral fan in 6-d lattice N
        sage: F.dim()
        3
        sage: [len(x) for x in F.cones()]
        [1, 4, 5, 2]

        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: F = rigidity_fan(P)
        sage: F
        Rational polyhedral fan in 10-d lattice N
        sage: F.ambient_dim()
        10
        sage: F.dim()
        7
        sage: [len(x) for x in F.cones()]
        [1, 11, 49, 116, 159, 127, 55, 10]

        sage: P = ((0,0), (1,1), (2,0), (0,2), (2,2))
        sage: F = rigidity_fan(P)
        sage: F
        Rational polyhedral fan in 10-d lattice N
        sage: F.dim()
        7
        sage: [len(x) for x in F.cones()]
        [1, 12, 56, 136, 189, 152, 66, 12]

        sage: P = ((0,0), (6,0), (2,2), (4,2), (3,3), (0,6), (6,6))
        sage: F = rigidity_fan(P)
        sage: F
        Rational polyhedral fan in 14-d lattice N
        sage: F.dim()
        11
        sage: [len(x) for x in F.cones()]
        [1, 31, 356, 2185, 8265, 20632, 35101, 41027, 32474, 16651, 4995, 666]

        sage: P = ((0,0), (6,0), (2,2), (3,3), (4,4), (0,6), (6,6))
        sage: WC = wiggly_complex(P)
        sage: F = rigidity_fan(P)
        sage: F
        Rational polyhedral fan in 14-d lattice N
        sage: F.dim()
        11
        sage: [len(x) for x in F.cones()]
        [1, 38, 492, 3272, 13083, 33942, 59290, 70528, 56446, 29124, 8760, 1168]
    """
    b = len(boundary_arcs(P))
    dict_arcs = dict((a,i) for (i,a) in enumerate(relevant_arcs(P)))
    rays = [rigidity_vector(P, a) for a in boundary_arcs(P)] + [rigidity_vector(P, a) for a in relevant_arcs(P)]
    WC = wiggly_complex(P)
    return Fan(cones=[list(range(b)) + [dict_arcs[a]+b for a in f] for f in WC.faces()[len(WC.faces())-2]], rays=rays)

@cached_function
def linear_dependence_rigidity_vectors(P, flip):
    r"""
    Return the unique linear dependence between the rigidity vectors of the arcs corresponding to two wiggly pseudotriangulations related by a flip.
    """
    aa = [list(set(flip[0]).difference(flip[1]))[0]] + [list(set(flip[1]).difference(flip[0]))[0]] + list(set(flip[0]).intersection(flip[1])) + boundary_arcs(P)
    #print(Matrix([rigidity_vector(P,a) for a in aa]))
    #print(kernel(Matrix([rigidity_vector(P,a) for a in aa])))
    return kernel(Matrix([rigidity_vector(P,a) for a in aa])).basis()[0]

@cached_function
def all_linear_dependences_rigidity_vectors(P):
    r"""
    Return all linear dependencies between rigidity vectors of the arcs corresponding to two wiggly pseudotriangulations related by a flip.

    EXAMPLES::
        sage: P = ((0,0), (1,0), (2,0),(3,0))
        sage: for x in all_linear_dependences_rigidity_vectors(P):
        ....:     print(x)
        ....:
        ((((0, 1, {}, {}), (0, 2, {}, {1}), (0, 3, {2}, {1})), ((0, 1, {}, {}), (2, 3, {}, {}), (0, 3, {2}, {1}))), (1, 2, 0, -1/2, -1/6, 0))
        ((((0, 1, {}, {}), (0, 2, {}, {1}), (0, 3, {2}, {1})), ((0, 2, {}, {1}), (0, 3, {2}, {1}), (1, 3, {2}, {}))), (1, 1/2, 0, -1/4, 0, -1/12))
        ((((0, 1, {}, {}), (0, 2, {}, {1}), (0, 2, {1}, {})), ((0, 1, {}, {}), (0, 2, {}, {1}), (0, 3, {2}, {1}))), (1, 1, -1, 0, -1/3, -2/3))
        ((((0, 1, {}, {}), (0, 2, {}, {1}), (0, 2, {1}, {})), ((0, 1, {}, {}), (0, 3, {1}, {2}), (0, 2, {1}, {}))), (1, 1, 0, -1, -2/3, -1/3))
        ((((0, 1, {}, {}), (0, 2, {}, {1}), (0, 2, {1}, {})), ((0, 2, {}, {1}), (1, 2, {}, {}), (0, 2, {1}, {}))), (1, 1, -1/4, -1/4, 0, 0))
        ((((0, 1, {}, {}), (0, 3, {1}, {2}), (0, 2, {1}, {})), ((0, 1, {}, {}), (0, 3, {1}, {2}), (2, 3, {}, {}))), (1, 2, 0, -1/2, 0, -1/6))
        ((((0, 1, {}, {}), (0, 3, {1}, {2}), (0, 2, {1}, {})), ((0, 3, {1}, {2}), (0, 2, {1}, {}), (1, 3, {}, {2}))), (1, 1/2, 0, -1/4, -1/12, 0))
        ((((0, 1, {}, {}), (0, 3, {1}, {2}), (2, 3, {}, {})), ((0, 1, {}, {}), (2, 3, {}, {}), (0, 3, {2}, {1}))), (1, 1, 0, 0, -1, -1))
        ((((0, 1, {}, {}), (0, 3, {1}, {2}), (2, 3, {}, {})), ((0, 3, {1}, {2}), (2, 3, {}, {}), (1, 3, {}, {2}))), (1, 1/2, 0, -1/4, -1/12, 0))
        ((((0, 1, {}, {}), (2, 3, {}, {}), (0, 3, {2}, {1})), ((2, 3, {}, {}), (0, 3, {2}, {1}), (1, 3, {2}, {}))), (1, 1/2, -1/4, 0, 0, -1/12))
        ((((0, 2, {}, {1}), (1, 2, {}, {}), (0, 2, {1}, {})), ((1, 2, {}, {}), (0, 2, {1}, {}), (1, 3, {}, {2}))), (1, 1, -2, 0, -1/2, -1/6))
        ((((0, 2, {}, {1}), (0, 3, {2}, {1}), (1, 3, {2}, {})), ((0, 2, {}, {1}), (1, 2, {}, {}), (1, 3, {2}, {}))), (1, 2, -1, -1, -1/6, -1/6))
        ((((0, 2, {}, {1}), (0, 3, {2}, {1}), (1, 3, {2}, {})), ((2, 3, {}, {}), (0, 3, {2}, {1}), (1, 3, {2}, {}))), (1, 2, 0, -1/2, -1/6, 0))
        ((((0, 2, {}, {1}), (1, 2, {}, {}), (1, 3, {2}, {})), ((0, 2, {}, {1}), (1, 2, {}, {}), (0, 2, {1}, {}))), (1, 1, 0, -2, -1/6, -1/2))
        ((((0, 2, {}, {1}), (1, 2, {}, {}), (1, 3, {2}, {})), ((1, 2, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2}))), (1, 1, 0, -2, -1/2, -1/6))
        ((((0, 3, {1}, {2}), (0, 2, {1}, {}), (1, 3, {}, {2})), ((1, 2, {}, {}), (0, 2, {1}, {}), (1, 3, {}, {2}))), (1, 2, -1, -1, -1/6, -1/6))
        ((((0, 3, {1}, {2}), (0, 2, {1}, {}), (1, 3, {}, {2})), ((0, 3, {1}, {2}), (2, 3, {}, {}), (1, 3, {}, {2}))), (1, 2, 0, -1/2, 0, -1/6))
        ((((1, 2, {}, {}), (0, 2, {1}, {}), (1, 3, {}, {2})), ((1, 2, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2}))), (1, 1, 0, -2, -1/6, -1/2))
        ((((0, 3, {1}, {2}), (2, 3, {}, {}), (1, 3, {}, {2})), ((2, 3, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2}))), (1, 1, -1, 0, -1/3, -2/3))
        ((((2, 3, {}, {}), (0, 3, {2}, {1}), (1, 3, {2}, {})), ((2, 3, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2}))), (1, 1, -1, 0, -2/3, -1/3))
        ((((1, 2, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2})), ((2, 3, {}, {}), (1, 3, {2}, {}), (1, 3, {}, {2}))), (1, 1, -1/4, -1/4, 0, 0))

        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: for x in all_linear_dependences_rigidity_vectors(P):
        ....:     print(x)
        ....:
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {}))), (1, 3, -1/2, 3/2, 0, -2, 0, -1))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {}, {}))), (1, 4/3, -1, -2, -4/3, 0, -4/3, -2/3))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 2/5, -1/5, 2/5, 0, -4/15, 0, -2/15))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {}))), (1, 1/3, -1/6, 1/2, -2/3, 0, -1/3, 0))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 1, {}, {}), (2, 3, {}, {}), (1, 2, {}, {}))), (1, 4/3, -1, -2, 0, -4/3, -2/3, -4/3))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 2/5, -1/5, 2/5, -4/15, 0, -2/15, 0))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {})), ((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 1/2, -1/6, -1/6, 0, 0, 0, 0))
        ((((0, 1, {}, {}), (2, 3, {}, {}), (1, 2, {}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {}, {}))), (3, 3, 3, -1, -1, -1, 0, 0))
        ((((0, 1, {}, {}), (2, 3, {}, {}), (1, 2, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {}), (1, 2, {}, {}))), (3, 3, 3, -1, -1, -1, 0, 0))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {}), (1, 2, {}, {}))), (3, 3, 3, -1, -1, -1, 0, 0))
        ((((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 12/5, -3/5, -1/5, 0, -8/5, 0, -4/5))
        ((((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {})), ((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 12/5, -3/5, -1/5, -8/5, 0, -4/5, 0))
        ((((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 5/2, -3/4, 1/2, -5/3, -1/3, -5/6, -1/6))
        ((((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 5/2, -3/4, 1/2, -1/3, -5/3, -1/6, -5/6))
        ((((1, 3, {}, {}), (1, 2, {}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {}), (1, 2, {}, {}))), (3, 4, -6, -6, 0, 0, -2, -2))

        sage: P = ((0,0), (1,1), (2,0), (0,2), (2,2))
        sage: for x in all_linear_dependences_rigidity_vectors(P):
        ....:     print(x)
        ....:
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {}))), (1, 8/3, -1/3, 0, 0, -4/3, 0, -4/3))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {1}, {}))), (1, 1, -2, -2, -1/2, 1/2, -3/2, -1/2))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (0, 4, {1}, {})), ((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 1, -1/3, 2/3, 0, -1/3, 0, -1/3))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {}))), (1, 3/8, -1/8, 0, -1/2, 0, -1/2, 0))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 1, {}, {}), (1, 2, {}, {}), (2, 3, {1}, {}))), (1, 1, -2, -2, 1/2, -1/2, -1/2, -3/2))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (1, 2, {}, {})), ((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 1, -1/3, 2/3, -1/3, 0, -1/3, 0))
        ((((0, 4, {}, {1}), (0, 1, {}, {}), (0, 4, {1}, {})), ((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 1, -1/4, -1/4, 0, 0, 0, 0))
        ((((0, 1, {}, {}), (1, 2, {}, {}), (2, 3, {1}, {})), ((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {1}, {}))), (1, 1, 2/3, -1/3, -1/3, -1/3, 0, 0))
        ((((0, 1, {}, {}), (1, 2, {}, {}), (2, 3, {1}, {})), ((2, 3, {}, {1}), (1, 2, {}, {}), (2, 3, {1}, {}))), (1, 3/8, 0, -1/8, -1/2, -1/2, 0, 0))
        ((((1, 3, {}, {}), (0, 1, {}, {}), (2, 3, {1}, {})), ((1, 3, {}, {}), (2, 3, {}, {1}), (2, 3, {1}, {}))), (1, 3/8, 0, -1/8, -1/2, -1/2, 0, 0))
        ((((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {}))), (1, 8/3, 0, -1/3, 0, -4/3, 0, -4/3))
        ((((0, 4, {}, {1}), (0, 4, {1}, {}), (1, 4, {}, {})), ((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 8/3, 0, -1/3, -4/3, 0, -4/3, 0))
        ((((1, 3, {}, {}), (0, 4, {1}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {1}), (1, 4, {}, {}))), (1, 1, -2, -2, -3/2, -1/2, -1/2, 1/2))
        ((((0, 4, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {})), ((2, 3, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {}))), (1, 1, -2, -2, -1/2, -3/2, 1/2, -1/2))
        ((((2, 3, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {})), ((2, 3, {}, {1}), (1, 2, {}, {}), (2, 3, {1}, {}))), (1, 3/8, -1/8, 0, 0, 0, -1/2, -1/2))
        ((((2, 3, {}, {1}), (1, 2, {}, {}), (1, 4, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {1}), (1, 4, {}, {}))), (1, 1, -1/3, 2/3, 0, 0, -1/3, -1/3))
        ((((2, 3, {}, {1}), (1, 2, {}, {}), (2, 3, {1}, {})), ((1, 3, {}, {}), (2, 3, {}, {1}), (2, 3, {1}, {}))), (1, 1, -1/4, -1/4, 0, 0, 0, 0))
        ((((1, 3, {}, {}), (2, 3, {}, {1}), (1, 4, {}, {})), ((1, 3, {}, {}), (2, 3, {}, {1}), (2, 3, {1}, {}))), (1, 3/8, -1/8, 0, 0, 0, -1/2, -1/2))
    """
    return [(e, linear_dependence_rigidity_vectors(P, e)) for e in wiggly_complex(P).flip_graph().edges(labels=None)]

@cached_function
def reduced_rigidity_fan(P):
    r"""
    Return the reduced rigidity fan of P.

    EXAMPLES::
        sage: P = ((0,0), (1,0), (2,0))
        sage: G = reduced_rigidity_fan(P)
        sage: G
        Rational polyhedral fan in 6-d lattice N
        sage: G.dim()
        2
        sage: [len(x) for x in G.cones()]
        [1, 2]

        sage: P = ((0,0), (1,1), (3,0), (0,3), (3,3))
        sage: F = rigidity_fan(P)
        sage: F
        Rational polyhedral fan in 10-d lattice N
        sage: F.ambient_dim()
        10
        sage: F.dim()
        6
        sage: [len(x) for x in F.cones()]
        [1, 7, 15, 10]
    """
    dict_arcs = dict((a,i) for (i,a) in enumerate(relevant_arcs(P)))
    rays = [rigidity_vector(P, a) for a in relevant_arcs(P)]
    WC = wiggly_complex(P)
    return Fan(cones=[[dict_arcs[a] for a in f] for f in WC.faces()[len(WC.faces())-2]], rays=rays)
