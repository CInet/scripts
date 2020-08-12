#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <set>
#include <cmath>
using namespace std;

// Count the number of 1s (set bits) in a binary string.
const size_t setBits(size_t n) {
    size_t ones = 0;
    while (n) {
        ones += n & 1;
        n >>= 1;
    } // end while

    return ones;
} // end function

// Get power set of N, with option to limit subset size.
template <class T>
const set<set<T>> P(const set<T>& N, const int r=-1) {
    size_t n = N.size();
    size_t p = pow(2, n);
    
    set<set<T>> powerSet;
    for (size_t i = 0; i < p; i++) {
        if (r >= 0 && setBits(i) != r)
            continue;

        set<T> subset;
        for (size_t j = 0; j < n; j++)
            if (i & (1 << j))
                subset.insert(*next(N.begin(), j));
        
        powerSet.insert(subset);
    } // end for

    return powerSet;
} // end function

// Check if matrix is for a partially directed acyclic graph (PDAG).
const bool isPDAG(const vector<vector<size_t>>& D) {
    size_t n = D.size();
    
    for (size_t a = 0; a < n; a++)
        for (size_t b = 0; b < n; b++)
            if (D[a][b] && D[b][a])
                return true;

    return false;
} // end function

// Check rule 1 of Verma & Pearl (1992) Phase 2a.
const bool rule1(const vector<vector<size_t>>& D, const size_t a, const size_t b, const size_t c) {
    return (D[a][b] && !D[b][a])
        && (D[b][c] &&  D[c][b]);
} // end function

// Check rule 2 of Verma & Pearl (1992) Phase 2a.
const bool rule2(const vector<vector<size_t>>& D, const size_t a, const size_t b, const size_t c) {
    return (D[a][b] && !D[b][a])
        && (D[b][c] && !D[c][b])
        && (D[a][c] &&  D[a][c]);
} // end function

// Check rule 3 of Verma & Pearl (1992) Phase 2a.
const bool rule3(const vector<vector<size_t>>& D, const size_t a, const size_t b, const size_t c, const size_t d) {
    return (D[a][b] &&  D[b][a])
        && (D[a][d] && !D[d][a])
        && (D[b][c] &&  D[c][b])
        && (D[b][d] &&  D[d][b])
        && (D[c][d] && !D[c][d]);
} // end function

// Check rule 4 of Verma & Pearl (1992) Phase 2a.
const bool rule4(const vector<vector<size_t>>& D, const size_t a, const size_t b, const size_t c, const size_t d) {
    return (D[a][b] &&  D[b][a])
        && (D[a][c] &&  D[c][a])
        && (D[b][c] &&  D[c][b])
        && (D[c][d] &&  D[d][c])
        && (D[d][a] && !D[a][d]);
} // end function

// Check if directed, simple graph is acyclic.
const bool acyclic(const vector<vector<size_t>>& D, const size_t j, const size_t steps=0) {
    if (steps >= D.size()) return false;

    bool isAcyclic = true;
    for(size_t k = 0; k < D[j].size(); k++)
        if (D[j][k] && !D[k][j])
            isAcyclic = isAcyclic && acyclic(D, k, steps + 1);

    return isAcyclic;
} // end function

// Check for closure in the sense of Verma & Pearl (1992) Phase 2b/c.
const bool closed(const vector<vector<size_t>>& D, const vector<vector<size_t>>& G) {
    size_t n = D.size();
    
    bool isAcyclic = true;
    for (size_t a = 0; a < n; a++)
        isAcyclic = isAcyclic && acyclic(D, a);

    bool newVee = false;
    for (size_t a = 0; a < n; a++) {
        for (size_t b = 0; b < n; b++) {
            for (size_t c = 0; c < n; c++) {
                // Vee structure = a -> b <- c where a is not adjacent to c.
                bool DVee = (D[a][b] && D[c][b] && !D[a][c] && !D[c][a]);
                bool GVee = (G[a][b] && G[c][b] && !G[a][c] && !G[c][a]);

                // Check XOR of DVee and GVee. If these aren't consistent, new vee structure was added.
                newVee = (DVee != GVee);
            } // end for
        } // end for
    } // end for

    return isAcyclic && !newVee;
} // end function

#endif