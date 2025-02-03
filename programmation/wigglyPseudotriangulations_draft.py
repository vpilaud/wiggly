@cached_function
def number_arcs(n):
    return add(2**j * (n+1-j) for j in range(n+1)) - 2

@cached_function
def arcs(n):
    return [(i, j, A, Set(range(i+1, j)).difference(A)) for j in range(n+2) for i in range(j) for A in Subsets(range(i+1, j))]

def tuplize(arc):
    return (arc[0], arc[1], tuple(sorted(arc[2])), tuple(sorted(arc[3])))

@cached_function
def relevant_arcs(n):
    return [(i, j, A, Set(range(i+1, j)).difference(A)) for j in range(n+2) for i in range(j) for A in Subsets(range(i+1, j)) if j-i < n+1 or len(A) not in [0, n]]

def are_crossing(arc1, arc2):
    (i1, j1, A1, B1) = arc1
    C1 = Set([i1, j1])
    (i2, j2, A2, B2) = arc2
    C2 = Set([i2, j2])
    return (A1.intersection(B2) or C1.intersection(B2) or A1.intersection(C2)) and (B1.intersection(A2) or C1.intersection(A2) or B1.intersection(C2))

@cached_function
def wiggly_crossing_graph(n):
    return Graph([arcs(n), are_crossing])

def compatible(arc1, arc2):
    (i1, j1, A1, B1) = arc1
    (i2, j2, A2, B2) = arc2
    return arc1 != arc2 and i1 != j2 and j1 != i2 and not are_crossing(arc1, arc2)

@cached_function
def wiggly_compatibililty_graph(n):
    return Graph([relevant_arcs(n), compatible])

def wiggly_incompatibility_number(arc1, arc2):
    (i1, j1, A1, B1) = arc1
    C1 = Set([i1, j1])
    (i2, j2, A2, B2) = arc2
    C2 = Set([i2, j2])
    X = (A1.intersection(B2).union(C1.intersection(B2)).union(A1.intersection(C2)))
    Y = (B1.intersection(A2).union(C1.intersection(A2)).union(B1.intersection(C2)))
    return int(i1 == j2) + int(i2 == j1) + len(Poset([[(x,y) for x in X for y in Y if x < y], lambda p1, p2: p1 != p2 and p1[0] >= p2[0] and p1[1] <= p2[1]]).minimal_elements()) + len(Poset([[(y,x) for x in X for y in Y if y < x], lambda p1, p2: p1 != p2 and p1[0] >= p2[0] and p1[1] <= p2[1]]).minimal_elements())

@cached_function
def incompatibility_numbers_dictionary(n):
    return dict([(tuplize(arc), add(wiggly_incompatibility_number(arc, arc2) for arc2 in relevant_arcs(n))) for arc in relevant_arcs(n)])

@cached_function
def wiggly_complex(n):
    G = wiggly_compatibililty_graph(n)
    G.relabel(tuplize)
    return G.clique_complex()

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
                for arcb in X[0]:
                    if not compatible(arc, arcb):
                        Y[2].add_edge(arc, arcb)
                if Y[2].clique_number() < k+1:
                    res.append(Y)
        return res
    
    return RecursivelyEnumeratedSet(seeds, succ, structure='forest')

@cached_function
def multi_wiggly_complex(n, k):
    return SimplicialComplex([x[0] for x in multiWigglyElements(n, k)])

@cached_function
def simplicial_asso(n):
    return Polyhedron([[(n+j-i+2)/(j-i+1) if i <= k <= j else -(j-i+2)/(n-j+i-1) for k in range(n)] for j in range(n) for i in range(j+1) if (i,j) != (0, n-1)])

def weird_set(arc, n):
    (i, j, A, B) = arc
    return Set([2*i, 2*j-1] + [2*a-1 for a in A] + [2*a for a in A]).difference(Set([0,2*n+1]))

def standardized_characteristic_vector(X, d):
    x = len(X)
    return [(d+x+1)/x if i in X else -(x+1)/(d-x) for i in range(1,d+1)]

@cached_function
def failed_wiggly_associahedron_1(n):
    return Polyhedron([standardized_characteristic_vector(weird_set(arc, n), 2*n) for arc in relevant_arcs(n)])

@cached_function
def special_permutations_among_all_permutations(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    """
    return [p for p in Permutations(2*n) if all(p(j) > 2*i for i in range(2,n+1) for j in range(p.inverse()(2*i-1)+1, p.inverse()(2*i))) and all(p(j) < 2*i for i in range(1,n) for j in range(p.inverse()(2*i)+1, p.inverse()(2*i-1)))]

@cached_function
def special_permutations_by_inversions(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    Here, they are characterised by their inversion sets.
    """
    return [p.inverse() for p in Permutations(2*n) if all(all(p(i) < p(2*j-1) or p(i) > p(2*j) for i in range(1, 2*j-1)) and all(p(2*j-1) < p(k) or p(2*j) > p(k) for k in range(2*j+1, 2*n+1)) for j in range(1,n+1))]

def flip(sp, j):
    r"""
    Flip an ascent j in a special permutation to get a new special permutation.
    """
    i = min(j-1, sp.inverse()(sp(j)+1)-1) if sp(j) % 2 else j-1
    k = max(j+1, sp.inverse()(sp(j+1)+1)) if sp(j+1) % 2 else j+1
    l = list(sp)
    return Permutation(l[:i] + l[j:k] + l[i:j] + l[k:])

@cached_function
def special_permutations_digraph(n):
    return DiGraph(dict((sp, [flip(sp, j) for j in sp.descents(positive=True)]) for sp in special_permutations(n)))

@cached_function
def special_permutations_graph(n):
    return Graph(special_permutations_digraph(n))

@cached_function
def special_permutations(n):
    r"""
    Return the set of permutations avoiding 2j-1 -- s -- 2j for s < 2j-1 and 2j -- b -- 2j-1 for b > 2j.
    
    EXAMPLE::
        sage: [len(special_permutations(n)) for n in range(1,7)]
        [2, 14, 176, 3232, 78384, 2366248]
    """
    seeds = [Permutation(range(1,2*n+1))]
    def succ(sp):
        return [flip(sp, j) for j in sp.descents(positive=True)]
    return list(RecursivelyEnumeratedSet(seeds, succ))

@cached_function
def special_permutations_lattice(n):
    return LatticePoset([special_permutations(n), lambda x, y: Set(x.inverse().inversions()).issubset(Set(y.inverse().inversions()))])

@cached_function
def special_permutations_singletons(n):
    return [sp for sp in special_permutations(n) if not any((sp(j) % 2 and sp.inverse()(sp(j)+1) < j) or (sp(j+1) % 2 and sp.inverse()(sp(j+1)+1) > j+1) for j in sp.descents(positive=True))]

@cached_function
def special_permutations_singletons_prefixes(n):
    return Set([Set(sp[:i+1]) for sp in special_permutations_singletons(n) for i in range(2*n-1)])

@cached_function
def failed_wiggly_associahedron_2(n):
    return Polyhedron(ieqs = [[-binomial(len(pref)+1,2)] + [int(i in pref) for i in range(1,2*n+1)] for pref in special_permutations_singletons_prefixes(n)], eqns = [[-binomial(2*n+1,2)] + [1]*(2*n)])

@cached_function
def failed_wiggly_associahedron_3(n):
    return Polyhedron(ieqs = [[-3**(len(pref))] + [int(i in pref) for i in range(1,2*n+1)] for pref in special_permutations_singletons_prefixes(n)], eqns = [[-3**(2*n)] + [1]*(2*n)])

def poset_from_special_permutation(sp):
    m = len(sp)
    invs = sp.inverse().inversions()
    decreasing_relations = [(j,i) for (i,j) in invs]
    increasing_relations = Set([(i,j) for j in range(1,m+1) for i in range(1,j)]).difference(Set(invs))
    forgotten_increasing_relations = Set([])
    for j in sp.descents(positive=True):
        i = min(j, sp.inverse()(sp(j)+1)) if sp(j) % 2 else j
        k = max(j+1, sp.inverse()(sp(j+1)+1)) if sp(j+1) % 2 else j+1
        increasing_relations = increasing_relations.difference(Set([(sp(x), sp(y)) for x in range(1, j+1) for y in range(j+1, k+1)])).difference(Set([(sp(x), sp(y)) for x in range(i, j+1) for y in range(j+1, m+1)]))
        forgotten_increasing_relations = forgotten_increasing_relations.union(Set([(sp(x), sp(k)) for x in range(1, i+1)])).union(Set([(sp(i), sp(y)) for y in range(k, m+1)]))
    increasing_relations = increasing_relations.union(forgotten_increasing_relations)
    return Poset([list(range(1,m+1)), list(increasing_relations) + list(decreasing_relations)])

def basis_vector(i,m):
    return vector([0]*(i-1) + [1] + [0]*(m-i))

def weird_vector(i,j,k,m):
    return vector([0]*(i-1) + [1/(j-i+1)]*(j-i+1) + [-1/(k-j)]*(k-j) + [0]*(m-k))

def incidence_cone(pos):
    n = len(pos)
    return Polyhedron(ieqs = [[0] + list(basis_vector(i,n) - basis_vector(j,n)) for (i,j) in pos.cover_relations()], eqns = [[0] + [1]*n])

@cached_function
def failed_wiggly_fan(n):
    return Fan([Cone(incidence_cone(poset_from_special_permutation(sp))) for sp in special_permutations(n)])

def wiggly_ineq_ascent(sp, j):
    r"""
    Return the inequality corresponding to an ascent j in a special permutation sp.
    """
    i = min(j, sp.inverse()(sp(j)+1)) if sp(j) % 2 else j
    k = max(j+1, sp.inverse()(sp(j+1)+1)) if sp(j+1) % 2 else j+1
    m = len(sp)
    #return [0] + list(basis_vector(sp(i),m) - basis_vector(sp(k),m)) # obviously not working
    #return [0] + list(add(basis_vector(sp(l),m) for l in range(i,j+1))/(j-i+1) - add(basis_vector(sp(l),m) for l in range(j+1,k+1))/(k-j))
    return [0] + list(basis_vector(sp(i),m) + basis_vector(sp(j),m) - basis_vector(sp(j+1),m) - basis_vector(sp(k),m))

def wiggly_ineq_descent(sp, j):
    r"""
    Return the inequality corresponding to a descent j in a special permutation sp.
    """
    i = j if sp(j) % 2 else min(j, sp.inverse()(sp(j)-1))
    k = j+1 if sp(j+1) % 2 else max(j+1, sp.inverse()(sp(j+1)-1))
    m = len(sp)
    #return [0] + list(basis_vector(sp(i),m) - basis_vector(sp(k),m)) # obviously not working
    #return [0] + list(add(basis_vector(sp(l),m) for l in range(i,j+1))/(j-i+1) - add(basis_vector(sp(l),m) for l in range(j+1,k+1))/(k-j))
    return [0] + list(basis_vector(sp(i),m) + basis_vector(sp(j),m) - basis_vector(sp(j+1),m) - basis_vector(sp(k),m))

def wiggly_cone(sp):
    r"""
    Return the cone corresponding to a special permutation sp.
    """
    print(sp, [wiggly_ineq_ascent(sp, j) for j in sp.descents(positive=True)] + [wiggly_ineq_descent(sp, j) for j in sp.descents(positive=False)])
    return Cone(Polyhedron(ieqs=[wiggly_ineq_ascent(sp, j) for j in sp.descents(positive=True)] + [wiggly_ineq_descent(sp, j) for j in sp.descents(positive=False)], eqns = [[0] + [1]*len(sp)]))
    
@cached_function
def failed_wiggly_fan_2(n):
    return Fan([wiggly_cone(sp) for sp in special_permutations(n)])

@cached_function
def wiggly_permutations(m):
    r"""
    Return the permutations avoiding 2j-1 -- small -- 2j and 2j -- big -- 2j-1.
    
    EXAMPLES::
        sage: [len(wiggly_permutations(i)) for i in range(9)]
        [1, 1, 2, 5, 14, 51, 176, 807, 3232]
    """
    return [p for p in Permutations(m) if all(p(j) > 2*i for i in range(2,m//2+1) for j in range(p.inverse()(2*i-1)+1, p.inverse()(2*i))) and all(p(j) < 2*i for i in range(1,m//2+1) for j in range(p.inverse()(2*i)+1, p.inverse()(2*i-1)))]

@cached_function
def wiggly_permutations_lattice(m):
    return LatticePoset([wiggly_permutations(m), lambda x, y: Set(x.inverse().inversions()).issubset(Set(y.inverse().inversions()))])

def remove_last(p):
    l = list(p)
    l.remove(len(l))
    return Permutation(l)

def arc_vector(arc, n, essential=True):
    (i, j, A, B) = arc
    positives = Set([2*i-1, 2*j] + [2*a-1 for a in A] + [2*a for a in A]).difference(Set([-1, 2*n+2]))
    negatives = Set([2*i, 2*j-1] + [2*b-1 for b in B] + [2*b for b in B]).difference(Set([0, 2*n+1]))
    vect = add(basis_vector(pos, 2*n) for pos in positives) - add(basis_vector(neg, 2*n) for neg in negatives) # + (len(positives) - len(negatives)) * vector([1]*(2*n))
    if essential: vect = essentialize_weight(vect)
    return vect

@cached_function
def wiggly_fan(n):
    arc_dict = dict([(tuplize(arc), k) for (k, arc) in enumerate(arcs(n))])
    return Fan(cones = [[arc_dict[arc] for arc in f] for f in wiggly_complex(n).faces()[2*n-2]], rays = [arc_vector(arc, n) for arc in arcs(n)])

@cached_function
def wiggly_associahedron(n):
    m = len(relevant_arcs(n))
    return Polyhedron(ieqs = [[incompatibility_numbers_dictionary(n)[tuplize(arc)]] + list(-arc_vector(arc, n)) for arc in relevant_arcs(n)])
