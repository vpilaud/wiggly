from collections import defaultdict
from functools import cmp_to_key

### WIGGLY COMPLEX ON ARCS

@cached_function
def number_arcs(n):
    r"""
    Return the number of arcs on n+2 points.
    This is add(2**j * (n+1-j) for j in range(n+1)) = 2**(n+2)-n-3

    EXAMPLES::
        sage: [number_arcs(n) for n in range(10)]
        [1, 4, 11, 26, 57, 120, 247, 502, 1013, 2036]
    """
    return 2**(n+2)-n-3

def number_relevant_arcs(n):
    r"""
    Return the number of arcs on n+2 points.

    EXAMPLES::
        sage: [number_relevant_arcs(n) for n in range(10)]
        [0, 2, 9, 24, 55, 118, 245, 500, 1011, 2034]
    """
    return max(0, number_arcs(n)-2)

@cached_function
def arcs(n):
    r"""
    Return all arcs on n+2 points.

    EXAMPLES::
        sage: arcs(2)
        [(0, 1, {}, {}),
         (0, 2, {}, {1}),
         (0, 2, {1}, {}),
         (1, 2, {}, {}),
         (0, 3, {}, {1, 2}),
         (0, 3, {1}, {2}),
         (0, 3, {2}, {1}),
         (0, 3, {1, 2}, {}),
         (1, 3, {}, {2}),
         (1, 3, {2}, {}),
         (2, 3, {}, {})]
    """
    return [(i, j, A, Set(range(i+1, j)).difference(A)) for j in range(n+2) for i in range(j) for A in Subsets(range(i+1, j))]

@cached_function
def relevant_arcs(n):
    r"""
    Return all relevant arcs on n+2 points.

    EXAMPLES::
        sage: relevant_arcs(2)
        [(0, 1, {}, {}),
         (0, 2, {}, {1}),
         (0, 2, {1}, {}),
         (1, 2, {}, {}),
         (0, 3, {1}, {2}),
         (0, 3, {2}, {1}),
         (1, 3, {}, {2}),
         (1, 3, {2}, {}),
         (2, 3, {}, {})]
    """
    return [(i, j, A, Set(range(i+1, j)).difference(A)) for j in range(n+2) for i in range(j) for A in Subsets(range(i+1, j)) if j-i < n+1 or len(A) not in [0, n]]


def tuplize(arc):
    r"""
    Tuplize an arc (to use it in graphs, posets, etc).

    EXAMPLES::
        sage: [tuplize(arc) for arc in arcs(2)]
        [(0, 1, (), ()),
         (0, 2, (), (1,)),
         (0, 2, (1,), ()),
         (1, 2, (), ()),
         (0, 3, (), (1, 2)),
         (0, 3, (1,), (2,)),
         (0, 3, (2,), (1,)),
         (0, 3, (1, 2), ()),
         (1, 3, (), (2,)),
         (1, 3, (2,), ()),
         (2, 3, (), ())]
    """
    return (arc[0], arc[1], tuple(sorted(arc[2])), tuple(sorted(arc[3])))

def are_non_pointed(arc1, arc2):
    r"""
    Check if two arcs arc1 and arc2 are non pointed.
    """
    return arc1[0] == arc2[1] or arc1[1] == arc2[0]

def are_crossing(arc1, arc2):
    r"""
    Check if two arcs arc1 and arc2 are crossing.
    """
    (i1, j1, A1, B1) = arc1
    C1 = Set([i1, j1])
    (i2, j2, A2, B2) = arc2
    C2 = Set([i2, j2])
    return (A1.intersection(B2) or C1.intersection(B2) or A1.intersection(C2)) and (B1.intersection(A2) or C1.intersection(A2) or B1.intersection(C2))

def compatible(arc1, arc2):
    r"""
    Check if two arcs arc1 and arc2 are compatible.
    """
    return arc1 != arc2 and not are_non_pointed(arc1, arc2) and not are_crossing(arc1, arc2)

def incompatibility_degree(arc1, arc2):
    r"""
    Return the incompatibility number of two arcs arc1 and arc2. This counts the number of ways arc1 and arc2 are incompatible, that is, their number of crossings if they cross, and 1 if they share an endpoint.
    """
    (i1, j1, A1, B1) = arc1
    C1 = Set([i1, j1])
    (i2, j2, A2, B2) = arc2
    C2 = Set([i2, j2])
    X = (A1.intersection(B2).union(C1.intersection(B2)).union(A1.intersection(C2)))
    Y = (B1.intersection(A2).union(C1.intersection(A2)).union(B1.intersection(C2)))
    return int(i1 == j2) + int(i2 == j1) + len(Poset([[(x,y) for x in X for y in Y if x < y], lambda p1, p2: p1 != p2 and p1[0] >= p2[0] and p1[1] <= p2[1]]).minimal_elements()) + len(Poset([[(y,x) for x in X for y in Y if y < x], lambda p1, p2: p1 != p2 and p1[0] >= p2[0] and p1[1] <= p2[1]]).minimal_elements())

@cached_function
def incompatibility_numbers_dictionary(n):
    r"""
    Return the dictionary of incompatibility numbers of the arcs. The incompatibility number of an arc is the sum of its incompatibility degrees with all arcs.

    EXAMPLES::
        sage: incompatibility_numbers_dictionary(2)
        {(0, 1, (), ()): 3,
         (0, 2, (), (1,)): 3,
         (0, 2, (1,), ()): 3,
         (1, 2, (), ()): 4,
         (0, 3, (1,), (2,)): 4,
         (0, 3, (2,), (1,)): 4,
         (1, 3, (), (2,)): 3,
         (1, 3, (2,), ()): 3,
         (2, 3, (), ()): 3}
    """
    return dict([(tuplize(arc), add(incompatibility_degree(arc, arc2) for arc2 in relevant_arcs(n))) for arc in relevant_arcs(n)])

@cached_function
def wiggly_compatibililty_graph(n):
    r"""
    Return the compatibility graph on relevant arcs on n+2 points.

    EXAMPLES::
        sage: wiggly_compatibililty_graph(2)
        Graph on 9 vertices
        sage: wiggly_compatibililty_graph(2).edges(labels=None)
        [((0, 1, {}, {}), (0, 2, {1}, {})), ((0, 1, {}, {}), (0, 2, {}, {1})), ((0, 1, {}, {}), (0, 3, {2}, {1})), ((0, 1, {}, {}), (0, 3, {1}, {2})), ((0, 1, {}, {}), (2, 3, {}, {})), ((0, 2, {1}, {}), (0, 3, {1}, {2})), ((0, 2, {1}, {}), (1, 2, {}, {})), ((0, 2, {1}, {}), (1, 3, {}, {2})), ((0, 2, {}, {1}), (0, 2, {1}, {})), ((0, 2, {}, {1}), (0, 3, {2}, {1})), ((0, 2, {}, {1}), (1, 2, {}, {})), ((0, 2, {}, {1}), (1, 3, {2}, {})), ((0, 3, {2}, {1}), (1, 3, {2}, {})), ((0, 3, {2}, {1}), (2, 3, {}, {})), ((0, 3, {1}, {2}), (1, 3, {}, {2})), ((0, 3, {1}, {2}), (2, 3, {}, {})), ((1, 2, {}, {}), (1, 3, {2}, {})), ((1, 2, {}, {}), (1, 3, {}, {2})), ((1, 3, {2}, {}), (2, 3, {}, {})), ((1, 3, {}, {2}), (1, 3, {2}, {})), ((1, 3, {}, {2}), (2, 3, {}, {}))]
    """
    return Graph([relevant_arcs(n), compatible])

@cached_function
def wiggly_complex(n):
    r"""
    Return the wiggly complex. It is the clique complex on the compatibility graph on relevant arcs.

    EXAMPLES::
        sage: [wiggly_complex(n).f_vector() for n in range(5)]
        [[1],
         [1, 2],
         [1, 9, 21, 14],
         [1, 24, 154, 396, 440, 176],
         [1, 55, 729, 4002, 10930, 15684, 11312, 3232]]
    """
    G = wiggly_compatibililty_graph(n)
    G.relabel(tuplize)
    return G.clique_complex()

@cached_function
def wiggly_polynomial(n):
    r"""
    Return the polynomial whose ith coefficient is the number of wiggly pseudotriangulations of n with i relevant arcs ending at n+1.

    EXAMPLE::
        sage: wiggly_polynomial(1)
        x + 1
        sage: wiggly_polynomial(2)
        3*x^3 + 5*x^2 + 4*x + 2
        sage: wiggly_polynomial(3)
        15*x^5 + 35*x^4 + 44*x^3 + 40*x^2 + 28*x + 14

    This can also be computed as follows:
        sage: F = 1
        ....: f(x) = x+1
        ....: for i in range(1, 10):
        ....:     F = F + f * y^i
        ....:     f = expand(factor((f(1) + x^2 * f * (2*x-3) + diff(f) * x^3 * (x-1)) / (1-x)^2) + diff(f) * x^3 + f * 2 * x^2)
        ....: F
        x |--> (34459425*x^17 + 218243025*x^16 + 745945200*x^15 + 1827625800*x^14 + 3595992400*x^13 + 6034818790*x^12 + 8959301120*x^11 + 12050255360*x^10 + 14923333360*x^9 + 17207539360*x^8 + 18611188864*x^7 + 18962789832*x^6 + 18223650336*x^5 + 16476603080*x^4 + 13897876544*x^3 + 10722777024*x^2 + 7205540800*x + 3602770400)*y^9 + (2027025*x^15 + 11486475*x^14 + 35315280*x^13 + 78228150*x^12 + 139784260*x^11 + 213823610*x^10 + 290120320*x^9 + 357161904*x^8 + 404873680*x^7 + 426519128*x^6 + 419495744*x^5 + 385156800*x^4 + 327939216*x^3 + 254236280*x^2 + 171068352*x + 85534176)*y^8 + (135135*x^13 + 675675*x^12 + 1843380*x^11 + 3643640*x^10 + 5837300*x^9 + 8033718*x^8 + 9825760*x^7 + 10901304*x^6 + 11100912*x^5 + 10423560*x^4 + 8994688*x^3 + 7020360*x^2 + 4732496*x + 2366248)*y^7 + (10395*x^11 + 45045*x^10 + 107100*x^9 + 185570*x^8 + 261800*x^7 + 318066*x^6 + 343136*x^5 + 333920*x^4 + 294144*x^3 + 231920*x^2 + 156768*x + 78384)*y^6 + (945*x^9 + 3465*x^8 + 7000*x^7 + 10360*x^6 + 12520*x^5 + 13006*x^4 + 11872*x^3 + 9520*x^2 + 6464*x + 3232)*y^5 + (105*x^7 + 315*x^6 + 520*x^5 + 630*x^4 + 620*x^3 + 514*x^2 + 352*x + 176)*y^4 + (15*x^5 + 35*x^4 + 44*x^3 + 40*x^2 + 28*x + 14)*y^3 + (3*x^3 + 5*x^2 + 4*x + 2)*y^2 + (x + 1)*y + 1

    Hence, the number of wiggly pseudotriangulations is computed as follows:
        sage: F = 1
        ....: f(x) = x+1
        ....: for i in range(1, 15):
        ....:     F = F + f(1) * y^i
        ....:     f = expand(factor((f(1) + x^2 * f * (2*x-3) + diff(f) * x^3 * (x-1)) / (1-x)^2) + diff(f) * x^3 + f * 2 * x^2)
        ....: F
        209288750022132980352*y^14 + 2681213937595142656*y^13 + 37206559614499840*y^12 + 563142033172480*y^11 + 9373542317760*y^10 + 173300710720*y^9 + 3602770400*y^8 + 85534176*y^7 + 2366248*y^6 + 78384*y^5 + 3232*y^4 + 176*y^3 + 14*y^2 + 2*y + 1

    """
    var('x')
    return add(x**len([arc for arc in pt if arc[1] == n+1]) for pt in wiggly_complex(n).facets())

@cached_function
def wiggly_series(n):
    r"""
    Return the bivariate series whose (i,n)th coefficient is the number of wiggly pseudotriangulations of n with i relevant arcs ending at n+1.

    EXAMPLES::
        sage: wiggly_series(3)
        (15*x^5 + 35*x^4 + 44*x^3 + 40*x^2 + 28*x + 14)*y^3 + (3*x^3 + 5*x^2 + 4*x + 2)*y^2 + (x + 1)*y + 1
    """
    var('y')
    return add(wiggly_polynomial(m) * y**m for m in range(n+1))

### WIGGLY PERMUTATIONS AND WIGGLY LATTICE

@cached_function
def wiggly_permutations_among_all_permutations(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    Here, they are obtained by considering all permutations and deleting those containing the patterns.
    See also wiggly_permutations_by_inversions(n) and wiggly_permutations(n).

    EXAMPLES::
        sage: wiggly_permutations_among_all_permutations(2)
        [[1, 2, 3, 4],
         [1, 2, 4, 3],
         [1, 3, 4, 2],
         [1, 4, 2, 3],
         [1, 4, 3, 2],
         [2, 1, 3, 4],
         [2, 1, 4, 3],
         [3, 4, 1, 2],
         [3, 4, 2, 1],
         [4, 1, 2, 3],
         [4, 1, 3, 2],
         [4, 2, 1, 3],
         [4, 3, 1, 2],
         [4, 3, 2, 1]]
    """
    return [p for p in Permutations(2*n) if all(p(j) > 2*i for i in range(2,n+1) for j in range(p.inverse()(2*i-1)+1, p.inverse()(2*i))) and all(p(j) < 2*i for i in range(1,n) for j in range(p.inverse()(2*i)+1, p.inverse()(2*i-1)))]

@cached_function
def wiggly_permutations_by_inversions(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    Here, they are characterised by their inversion sets.
    See also wiggly_permutations_among_all_permutations(n) and wiggly_permutations(n).

    EXAMPLES::
        sage: wiggly_permutations_by_inversions(2)
        [[1, 2, 3, 4],
         [1, 2, 4, 3],
         [1, 4, 2, 3],
         [1, 3, 4, 2],
         [1, 4, 3, 2],
         [2, 1, 3, 4],
         [2, 1, 4, 3],
         [4, 1, 2, 3],
         [4, 1, 3, 2],
         [4, 2, 1, 3],
         [3, 4, 1, 2],
         [4, 3, 1, 2],
         [3, 4, 2, 1],
         [4, 3, 2, 1]]
    """
    return [p.inverse() for p in Permutations(2*n) if all(all(p(i) < p(2*j-1) or p(i) > p(2*j) for i in range(1, 2*j-1)) and all(p(2*j-1) < p(k) or p(2*j) > p(k) for k in range(2*j+1, 2*n+1)) for j in range(1,n+1))]

def flip(wp, j):
    r"""
    Flip an ascent j in a wiggly permutation to get a new wpecial permutation.
    """
    i = min(j-1, wp.inverse()(wp(j)+1)-1) if wp(j) % 2 else j-1
    k = max(j+1, wp.inverse()(wp(j+1)+1)) if wp(j+1) % 2 else j+1
    l = list(wp)
    return Permutation(l[:i] + l[j:k] + l[i:j] + l[k:])

@cached_function
def wiggly_permutations(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    Here, they are constructed by exploring the flip graph.
    See also wiggly_permutations_among_all_permutations(n) and wiggly_permutations_by_inversions(n).

    EXAMPLE::
        sage: wiggly_permutations(2)
        [[1, 2, 3, 4],
         [2, 1, 3, 4],
         [1, 3, 4, 2],
         [1, 2, 4, 3],
         [3, 4, 2, 1],
         [2, 1, 4, 3],
         [3, 4, 1, 2],
         [1, 4, 3, 2],
         [1, 4, 2, 3],
         [4, 3, 2, 1],
         [4, 2, 1, 3],
         [4, 3, 1, 2],
         [4, 1, 3, 2],
         [4, 1, 2, 3]]
        sage: [len(wiggly_permutations(n)) for n in range(1,7)]
        [2, 14, 176, 3232, 78384, 2366248]
    """
    seeds = [Permutation(range(1,2*n+1))]
    def succ(wp):
        return [flip(wp, j) for j in wp.descents(positive=True)]
    return list(RecursivelyEnumeratedSet(seeds, succ))

@cached_function
def wiggly_permutations_digraph(n):
    r"""
    Return the graph of oriented flips on wiggly permutations.

    EXAMPLES::
        sage: wiggly_permutations_digraph(2)
        Digraph on 14 vertices
        sage: wiggly_permutations_digraph(2).in_degree_sequence()
        [3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0]
    """
    return DiGraph(dict((wp, [flip(wp, j) for j in wp.descents(positive=True)]) for wp in wiggly_permutations(n)))

@cached_function
def wiggly_permutations_graph(n):
    r"""
    Return the graph of unoriented flips on wiggly permutations.

    EXAMPLES::
        sage: wiggly_permutations_graph(2)
        Graph on 14 vertices
        sage: wiggly_permutations_graph(2).degree_histogram()
        [0, 0, 0, 14]
    """
    return Graph(wiggly_permutations_digraph(n))

@cached_function
def wiggly_permutations_lattice(n):
    r"""
    Return the lattice on wiggly permutations. This is the sublattice of the weak order induced by wiggly permutations.
    """
    return LatticePoset(wiggly_permutations_digraph(n))
    # equivalent to the following:
    # return LatticePoset([wiggly_permutations(n), lambda x, y: Set(x.inverse().inversions()).issubset(Set(y.inverse().inversions()))])

@cached_function
def c_vectors(wp, essential=True):
    r"""
    Return the c-vectors of the wiggly permutation.
    """
    res = []
    m = len(wp)
    for i in range(1,m):
        if wp(i) < wp(i+1):
            v = basis_vector(wp(i), m)/2 + basis_vector(wp(i)+1 if wp(i) % 2 and wp.inverse()(wp(i)+1) < i else wp(i), m)/2 - basis_vector(wp(i+1), m)/2 - basis_vector(wp(i+1)+1 if wp(i+1) % 2 and wp.inverse()(wp(i+1)+1) > i+1 else wp(i+1), m)/2 + add((-1)**(wp.inverse()(2*x) < i) * (basis_vector(2*x-1, m) - basis_vector(2*x, m))/2 for x in range((wp(i)+1)//2 + 1, (wp(i+1)+1)//2))
        else:
            v = basis_vector(wp(i), m)/2 + basis_vector(wp(i)-1 if not(wp(i) % 2) and wp.inverse()(wp(i)-1) < i else wp(i), m)/2 - basis_vector(wp(i+1), m)/2 - basis_vector(wp(i+1)-1 if not(wp(i+1) % 2) and wp.inverse()(wp(i+1)-1) > i+1 else wp(i+1), m)/2 + add((-1)**(wp.inverse()(2*x) < i) * (basis_vector(2*x, m) - basis_vector(2*x-1, m))/2 for x in range((wp(i+1)+1)//2 + 1, (wp(i)+1)//2))
        if essential: v = essentialize_root(v)
        res.append(v)
    return res

### The next few functions extend wiggly permutations to even numbers of elements.
### It is just a curiosity, we do not need it in the paper.

@cached_function
def wiggly_permutations_odd_even(m):
    r"""
    Return the permutations avoiding 2j-1 -- small -- 2j and 2j -- big -- 2j-1.
    
    EXAMPLES::
        sage: [len(wiggly_permutations(i)) for i in range(9)]
        [1, 1, 2, 5, 14, 51, 176, 807, 3232]
    """
    return [p for p in Permutations(m) if all(p(j) > 2*i for i in range(2,m//2+1) for j in range(p.inverse()(2*i-1)+1, p.inverse()(2*i))) and all(p(j) < 2*i for i in range(1,m//2+1) for j in range(p.inverse()(2*i)+1, p.inverse()(2*i-1)))]

@cached_function
def wiggly_permutations_lattice_odd_even(m):
    return LatticePoset([wiggly_permutations_odd_even(m), lambda x, y: Set(x.inverse().inversions()).issubset(Set(y.inverse().inversions()))])

def remove_last(p):
    l = list(p)
    l.remove(len(l))
    return Permutation(l)

### BIJECTION FROM WIGGLY PSEUDOTRIANGULATIONS TO WIGGLY PERMUTATIONS

def lower_than(arc1, arc2):
    r"""
    Is arc1 below arc2?
    """
    (i1, j1, A1, B1) = arc1
    C1 = Set([i1, j1])
    (i2, j2, A2, B2) = arc2
    C2 = Set([i2, j2])
    return bool(Set(B1).intersection(Set(A2)) or C1.intersection(Set(A2)) or Set(B1).intersection(C2))

def ordered_arcs(pt):
    r"""
    Return the total order on the arcs of a wiggly pseudotriangulation pt
    """
    return Poset([pt, lower_than]).linear_extension()

def wpt2wp(pt):
    r"""
    Return the wiggly permutation corresponding to a wiggly pseudotriangulation pt.

    EXAMPLES::
        sage: n = 3; Set([wpt2wp(pt) for pt in wiggly_complex(n).facets()]) == Set(wiggly_permutations(n))
        True
        sage: [(wpt2wp(pt), ordered_arcs(pt)) for pt in WC2.facets()]
        [([4, 3, 1, 2], [(0, 2, (), (1,)), (0, 3, (2,), (1,)), (1, 3, (2,), ())]),
         ([4, 2, 1, 3], [(0, 2, (), (1,)), (0, 1, (), ()), (0, 2, (1,), ())]),
         ([3, 4, 1, 2], [(2, 3, (), ()), (0, 3, (2,), (1,)), (1, 3, (2,), ())]),
         ([1, 4, 2, 3], [(1, 3, (), (2,)), (1, 2, (), ()), (0, 2, (1,), ())]),
         ([1, 2, 3, 4], [(1, 3, (), (2,)), (0, 3, (1,), (2,)), (2, 3, (), ())]),
         ([4, 1, 2, 3], [(0, 2, (), (1,)), (1, 2, (), ()), (0, 2, (1,), ())]),
         ([2, 1, 3, 4], [(0, 1, (), ()), (0, 3, (1,), (2,)), (2, 3, (), ())]),
         ([1, 3, 4, 2], [(1, 3, (), (2,)), (2, 3, (), ()), (1, 3, (2,), ())]),
         ([3, 4, 2, 1], [(2, 3, (), ()), (0, 3, (2,), (1,)), (0, 1, (), ())]),
         ([2, 1, 4, 3], [(0, 1, (), ()), (0, 3, (1,), (2,)), (0, 2, (1,), ())]),
         ([4, 3, 2, 1], [(0, 2, (), (1,)), (0, 3, (2,), (1,)), (0, 1, (), ())]),
         ([1, 4, 3, 2], [(1, 3, (), (2,)), (1, 2, (), ()), (1, 3, (2,), ())]),
         ([1, 2, 4, 3], [(1, 3, (), (2,)), (0, 3, (1,), (2,)), (0, 2, (1,), ())]),
         ([4, 1, 3, 2], [(0, 2, (), (1,)), (1, 2, (), ()), (1, 3, (2,), ())])]    """
    pt = ordered_arcs(pt)
    n = (len(pt)+1)//2
    opening = [tuple([]) for i in range(n+2)]
    closing = [tuple([]) for i in range(n+2)]
    for arc in pt:
        if not opening[arc[0]]: opening[arc[0]] = arc
        if not opening[arc[1]]: opening[arc[1]] = arc
    for arc in pt[::-1]:
        if not closing[arc[0]]: closing[arc[0]] = arc
        if not closing[arc[1]]: closing[arc[1]] = arc
    opening_dict = defaultdict(list)
    closing_dict = defaultdict(list)
    for i in range(1,n+1):
        opening_dict[opening[i]].append(i)
        closing_dict[closing[i]].append(i)
    res = []
    for arc in pt:
        if opening_dict[arc]:
            i = opening_dict[arc][0]
            res.append(2*i-1 if i == arc[0] else 2*i)
        if closing_dict[arc]:
            i = closing_dict[arc][0]
            res.append(2*i if i == arc[0] else 2*i-1)
    return Permutation(res)

### BIJECTION FROM WIGGLY PERMUTATIONS TO WIGGLY PSEUDOTRIANGULATIONS

def wp2wpt(wp):
    r"""
    
    """
    res = []
    wpinv = wp.inverse()
    n = len(wp)//2
    for k in range(1,len(wp)):
        i = max([0] + [i for i in range(1,n+1) if wpinv(2*i-1) <= k < wpinv(2*i)])
        j = min([n+1] + [j for j in range(1,n+1) if wpinv(2*j) <= k < wpinv(2*j-1)])
        A = Set([l for l in range(i+1,j) if wpinv(2*l) <= k])
        B = Set([l for l in range(i+1,j) if wpinv(2*l) > k])
        res.append((i, j, A, B))
    return res

### WIGGLYHEDRON

def basis_vector(i,m):
    r"""
    Return the i-th vector of the standard basis.
    """
    return vector([0]*(i-1) + [1] + [0]*(m-i))

def g_vector(arc, n, essential=True):
    r"""
    Return the g-vector of an arc.
    This is defined as 00000+-......-+00000 where the ...... have  to be replaced by ++ if passing above and -- if passing below.
    If essential, then we essentialize.

    EXAMPLES::
        sage: [(arc, g_vector(arc, 2, essential=False)) for arc in arcs(n)]
        [((0, 1, {}, {}), (-1, 1, 0, 0)),
         ((0, 2, {}, {1}), (-1, -1, -1, 1)),
         ((0, 2, {1}, {}), (1, 1, -1, 1)),
         ((1, 2, {}, {}), (1, -1, -1, 1)),
         ((0, 3, {}, {1, 2}), (-1, -1, -1, -1)),
         ((0, 3, {1}, {2}), (1, 1, -1, -1)),
         ((0, 3, {2}, {1}), (-1, -1, 1, 1)),
         ((0, 3, {1, 2}, {}), (1, 1, 1, 1)),
         ((1, 3, {}, {2}), (1, -1, -1, -1)),
         ((1, 3, {2}, {}), (1, -1, 1, 1)),
         ((2, 3, {}, {}), (0, 0, 1, -1))]
        sage: [(arc, g_vector(arc, 2)) for arc in arcs(n)]
        [((0, 1, {}, {}), (-1, 1/2, 0)),
         ((0, 2, {}, {1}), (0, 0, -1)),
         ((0, 2, {1}, {}), (0, 1, -1)),
         ((1, 2, {}, {}), (1, 0, -1)),
         ((0, 3, {}, {1, 2}), (0, 0, 0)),
         ((0, 3, {1}, {2}), (0, 1, 0)),
         ((0, 3, {2}, {1}), (0, -1, 0)),
         ((0, 3, {1, 2}, {}), (0, 0, 0)),
         ((1, 3, {}, {2}), (1, 0, 0)),
         ((1, 3, {2}, {}), (1, -1, 0)),
         ((2, 3, {}, {}), (0, -1/2, 1))]
    """
    (i, j, A, B) = arc
    positives = Set([2*i-1, 2*j] + [2*a-1 for a in A] + [2*a for a in A]).difference(Set([-1, 2*n+2]))
    negatives = Set([2*i, 2*j-1] + [2*b-1 for b in B] + [2*b for b in B]).difference(Set([0, 2*n+1]))
    vect = add(basis_vector(pos, 2*n) for pos in positives) - add(basis_vector(neg, 2*n) for neg in negatives) # + (len(positives) - len(negatives)) * vector([1]*(2*n))
    if essential: vect = essentialize_weight(vect)
    return vect

@cached_function
def wiggly_fan(n):
    r"""
    Return the wiggly fdan.
    It is constructed using the g-vectors as support, and the wiggly complex as cones.

    EXAMPLES::
        sage: wiggly_fan(2)
        Rational polyhedral fan in 3-d lattice N
    """
    arc_dict = dict([(tuplize(arc), k) for (k, arc) in enumerate(arcs(n))])
    return Fan(cones = [[arc_dict[arc] for arc in f] for f in wiggly_complex(n).faces()[2*n-2]], rays = [g_vector(arc, n) for arc in arcs(n)])

@cached_function
def wigglyhedron_from_inequalities(n, essential=True):
    r"""
    Return the wigglyhedron.
    Here, it is constructed using its inequality description.
    There is a facet for each relevant arc: the normal vector is the g-vector of the arc, and the right hand side is the incompatibility number of the arc.

    EXAMPLES::
        sage: wigglyhedron_from_inequalities(1)
        A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
        sage: wigglyhedron_from_inequalities(2)
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 14 vertices
        sage: wigglyhedron_from_inequalities(3)
        A 5-dimensional polyhedron in QQ^5 defined as the convex hull of 176 vertices
        sage: wigglyhedron_from_inequalities(4)
        A 7-dimensional polyhedron in QQ^7 defined as the convex hull of 3232 vertices
    """
    m = len(relevant_arcs(n))
    return Polyhedron(ieqs = [[incompatibility_numbers_dictionary(n)[tuplize(arc)]] + list(-g_vector(arc, n, essential=essential)) for arc in relevant_arcs(n)])

def wiggly_vertex(pt, essential=True):
    r"""
    Return the vertex corresponding to the wiggly pseudotriangulation pt.
    """
    pt = ordered_arcs(pt)
    n = (len(pt)+1)//2
    cvect = c_vectors(wpt2wp(pt), essential=essential)
    return add([incompatibility_numbers_dictionary(n)[tuplize(pt[i])] * cvect[i] for i in range(len(pt))])

@cached_function
def wigglyhedron_from_vertices(n, essential=True):
    r"""
    Return the wigglyhedron.
    Here, it is constructed using its vertex description.
    There is a vertex for each facet of the wiggly complex: it is the sum of the c-vectors multiplied by the incompatibility numbers.

    EXAMPLES::
    """
    return Polyhedron([wiggly_vertex(pt, essential=essential) for pt in wiggly_complex(n).facets()])

@cached_function
def wigglyhedron(n, method='inequalities', essential=True):
    r"""
    Return the wigglyhedron, using wigglyhedron_from_inequalities or wigglyhedron_from_vertices depending on the parameter method.
    """
    if method == 'inequalities': return wigglyhedron_from_inequalities(n, essential=essential)
    if method == 'vertices': return wigglyhedron_from_vertices(n, essential=essential)
    print('method should be inequalities or vertices')

@cached_function
def directed_skeleton_wigglyhedron(n, essential=True, direction=None):
    r"""
    Return the skeleton of the wigglyhedron directed in the given direction.

    EXAMPLES::
        sage: n = 3; directed_skeleton_wigglyhedron(n).is_isomorphic(wiggly_permutations_digraph(n))
        True
        sage: n = 3; directed_skeleton_wigglyhedron(n, essential=False).is_isomorphic(wiggly_permutations_digraph(n))
        True
    """
    if not direction:
        direction = vector(add([[4*n*i,4*n*i] for i in range(1,n+1)], [])) + vector([i for i in range(1,2*n+1)])
    if essential:
        return directed_skeleton(wigglyhedron(n), essentialize_weight(direction))
    else:
        return directed_skeleton(wigglyhedron(n, method='vertices', essential=False), direction)
True

### CAMBRIAN LATTICES

@cached_function
def cambrian_face(signature):
    r"""
    Return a collection of wiggly arcs corresponding to a given signature.
    """
    n = len(signature)
    return Set([(0, i+1, tuple([]), tuple(range(1,i+1))) for i in range(n) if signature[i] < 0] + [(0, i+1, tuple(range(1,i+1)), tuple([])) for i in range(n) if signature[i] > 0]) # left - left
    # return Set([(i+1, n+1, tuple([]), tuple(range(i+2,n+1))) for i in range(n) if signature[i] < 0] + [(0, i+1, tuple(range(1,i+1)), tuple([])) for i in range(n) if signature[i] > 0]) # right - left
    # return Set([(0, i+1, tuple([]), tuple(range(1,i+1))) for i in range(n) if signature[i] < 0] + [(i+1, n+1, tuple(range(i+2,n+1)), tuple([])) for i in range(n) if signature[i] > 0]) # left - right
    # return Set([(i+1, n+1, tuple([]), tuple(range(i+2,n+1))) for i in range(n) if signature[i] < 0] + [(i+1, n+1, tuple(range(i+2,n+1)), tuple([])) for i in range(n) if signature[i] > 0]) # right - right

@cached_function
def cambrian_complex(signature):
    r"""
    Return the subcomplex of the wiggly complex corresponding to the given Cambrian associahedron.
    """
    return wiggly_complex(len(signature)).link(cambrian_face(signature))

@cached_function
def cambrian_wiggly_permutations(signature):
    r"""
    Return the set of Cambrian wiggly permutations for the given signature.
    """
    n = len(signature)
    return [wp for wp in wiggly_permutations(n) if all(wp.inverse()(2*j+2) < wp.inverse()(i+1) for j in range(n) if signature[j] < 0 for i in range(2*j+1)) and all(wp.inverse()(i+1) <= wp.inverse()(2*j+1) for j in range(n) if signature[j] > 0 for i in range(2*j+2))] # left - left
    # return [wp for wp in wiggly_permutations(n) if all(wp.inverse()(2*j+1) < wp.inverse()(k+1) for j in range(n) if signature[j] < 0 for k in range(2*j+1, 2*n)) and all(wp.inverse()(i+1) <= wp.inverse()(2*j+1) for j in range(n) if signature[j] > 0 for i in range(2*j+2))] # right - left
    # return [wp for wp in wiggly_permutations(n) if all(wp.inverse()(2*j+2) < wp.inverse()(i+1) for j in range(n) if signature[j] < 0 for i in range(2*j+1)) and all(wp.inverse()(i+1) <= wp.inverse()(2*j+2) for j in range(n) if signature[j] > 0 for i in range(2*j, 2*n))] # left - right
    # return [wp for wp in wiggly_permutations(n) if all(wp.inverse()(2*j+1) < wp.inverse()(k+1) for j in range(n) if signature[j] < 0 for k in range(2*j+1, 2*n)) and all(wp.inverse()(i+1) <= wp.inverse()(2*j+2) for j in range(n) if signature[j] > 0 for i in range(2*j, 2*n))] # right - right

@cached_function
def cambrian_wiggly_permutations_digraph(signature):
    r"""
    Return the graph of oriented flips on Cambrian wiggly permutations for the given signature.

    EXAMPLES::
    """
    return DiGraph(dict((wp, [flip(wp, j) for j in wp.descents(positive=True)  if flip(wp, j) in cambrian_wiggly_permutations(signature)]) for wp in cambrian_wiggly_permutations(signature)))

@cached_function
def cambrian_wiggly_permutations_lattice(signature):
    r"""
    Return the lattice on Cambrian wiggly permutations for the given signature.

    EXAMPLES::
        sage: cambrian_wiggly_permutations_lattice((1,-1,1,-1,1)).is_isomorphic(permutree_congruence([1,-1,1,-1,1]).quotient())
        True
    """
    return LatticePoset(cambrian_wiggly_permutations_digraph(signature))

@cached_function
def cambrian_wiggly_permutations_sublattice(signature):
    r"""
    Return the sublattice of the wiggly lattice corresponding to the given Cambrian lattice.
    """
    n = len(signature)
    f = cambrian_face(signature)
    return wiggly_permutations_lattice(n).sublattice([wpt2wp(list(f) + list(g)) for g in cambrian_complex(signature).facets()])

### OTHER ASSOCIAHEDRA

def multigraph(d,n):
    r"""
    Return the multigraph corresponding to a wiggly dissection d.
    """
    G = Graph(n+2, multiedges=True)
    G.add_edges([[a[0], a[1]] for a in d] + [[0,n+1],[0,n+1]])
    return G

@cached_function
def wiggly_trees(n):
    r"""
    Return the wiggly dissections whose multigraph is a tree.

    EXAMPLES::
        sage: wiggly_trees(2)
        [((1, 2, (), ()), (1, 3, (2,), ())),
         ((0, 2, (1,), ()), (0, 1, (), ())),
         ((1, 2, (), ()), (0, 2, (1,), ())),
         ((2, 3, (), ()), (1, 3, (), (2,))),
         ((0, 2, (1,), ()), (1, 3, (), (2,))),
         ((1, 2, (), ()), (1, 3, (), (2,))),
         ((2, 3, (), ()), (1, 3, (2,), ())),
         ((2, 3, (), ()), (0, 1, (), ())),
         ((0, 2, (), (1,)), (1, 3, (2,), ())),
         ((0, 2, (), (1,)), (0, 1, (), ())),
         ((0, 2, (), (1,)), (1, 2, (), ()))]
    """
    return [d for d in wiggly_complex(n).faces()[n-1] if multigraph(d,n).is_connected()]

@cached_function
def wiggly_associahedron(wt, n):
    r"""
    Return the subcomplex of the wiggly complex corresponding to the given wiggly tree.

    EXAMPLES::
    """
    return wiggly_complex(n).link(wt)

@cached_function
def has_long_arc(wt, n):
    r"""
    Check whether the wiggly associahedron of wt contains a long arc.
    """
    return any([a[0] == 0 and a[1] == n+1 for a in wiggly_associahedron(wt, n).vertices()])
    
### MULTI WIGGLY COMPLEXES
### This is very much in progress. Does not work at the moment.

def lt(arc1, arc2):
    (i1, j1, A1, B1) = arc1
    (i2, j2, A2, B2) = arc2
    if i1 != i2: return i1 < i2
    if j1 != j2: return j1 < j2
    if A1 != A2: return min(A1.symmetric_difference(A2)) in A1
    if B1 != B2: return min(B1.symmetric_difference(B2)) in B1
    return False

@cached_function
def multi_wiggly_elements(n, k):
    seeds = [([arc], arc, Graph([[arc], []])) for arc in arcs(n)]
    
    def succ(X):
        res = []
        for arc in arcs(n):
            if lt(arc, X[1]):
                Y = (copy(X[0]), arc, copy(X[2]))
                Y[0].append(arc)
                if not any(are_non_pointed(arc, arcb) for arcb in X[0]):
                    for arcb in X[0]:
                        if not compatible(arc, arcb):
                            Y[2].add_edge(arc, arcb)
                if Y[2].clique_number() < k+1:
                    res.append(Y)
        return res
    
    return RecursivelyEnumeratedSet(seeds, succ, structure='forest')

@cached_function
def multi_wiggly_complex(n, k):
    return SimplicialComplex([x[0] for x in multi_wiggly_elements(n, k)])
