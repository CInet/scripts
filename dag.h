#ifndef DAG_H
#define DAG_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <cmath>
#include <algorithm>

#include "util.h"
#include "ci.h"

using namespace std;

class DAG {
    private:
        vector<vector<size_t>> adjacencyMatrix;

        const vector<vector<size_t>> generatePDAG(const CI& M, const size_t n);
        const vector<vector<size_t>> extendToDAG(const vector<vector<size_t>>& G);
        const bool active(const set<size_t>& L, const size_t j);
        const vector<vector<bool>> legal(const set<size_t>& J, const vector<bool>& descendent);
        const map<size_t, bool> reachable(const set<size_t>& L, const set<size_t>& J, const vector<bool>& descendent);
        const set<size_t> maximalIndependentSet(const set<size_t>& J, const set<size_t>& L);
    
    public:
        DAG();
        DAG(const CI& M, const size_t n);

        const vector<vector<size_t>> matrix();

        vector<size_t> topologicalOrder();
        const bool consistent(const CI& M);
        const CI toCI();
        const vector<vector<size_t>> essentialGraph();
}; // end class

// Utility functions.

// Phase 1: Generate a PDAG if possible.
const vector<vector<size_t>> DAG::generatePDAG(const CI& M, const size_t n) {
    vector<vector<size_t>> G(n, vector<size_t>(n, 0));
    map<set<size_t>, set<size_t>> S;

    set<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.insert(i);

    map<set<size_t>, set<set<size_t>>> CIStruct = M.structure;

    // Look at all pairs of nodes and check if there is a CI statement between them.
    // If not, draw an undirected edge.
    // If so, map the pair with the respective predicate sets.
    set<set<size_t>> pairs = P(N, 2);
    for (set<size_t> I : pairs) {
        auto query = CIStruct.find(I);
        size_t a = *(I.begin());
        size_t b = *next(I.begin());
        
        if (query == CIStruct.end()) {
            G[a][b] = G[b][a] = 1;
        } else {
            for (set<size_t> evidence : query->second) {
                S[evidence].insert(a);
                S[evidence].insert(b);
            } // end for
        } // end if-else
    } // end for

    for (set<size_t> I : pairs) {
        size_t a = *(I.begin());
        size_t b = *next(I.begin());
        
        // Skip adjacent nodes.
        if (G[a][b] || G[b][a]) continue;

        // Check over all predicate sets.
        for (set<size_t> J : CIStruct[I]) {
            // Complement set.
            set<size_t> C = N;
            set_difference(N.begin(), N.end(), S[J].begin(), S[J].end(), inserter(C, C.end()));

            for (size_t c : C) {
                if ((G[a][c] || G[c][a]) && (G[b][c] || G[c][b])) {
                    if ((G[c][a] && !G[a][c]) || (G[c][b] && !G[b][c]))
                        throw 20;
                    else
                        G[c][a] = G[c][b] = 0;
                } // end if
            } // end for
        } // end for
    } // end for

    return G;
} // end function

// Phase 2: Extend G to a DAG, if possible.
const vector<vector<size_t>> DAG::extendToDAG(const vector<vector<size_t>>& G) {
    stack<pair<pair<size_t, size_t>, vector<vector<size_t>>>> C;
    vector<vector<size_t>> D = G;
    
    // Create a ground set N.
    size_t n = D.size();
    vector<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.push_back(i);

    // Store last pair that was arbitrarily directed.
    pair<size_t, size_t> lastDirected;

    while (isPDAG(D)) {
        // Iterate all triplets in D.
        for (size_t a = 0; a < n; a++) {
            for (size_t b = 0; b < n; b++) {
                if (b == a) continue;

                for (size_t c = 0; c < n; c++) {
                    if (c == a || c == b) continue;

                    if (rule1(D, a, b, c))
                        D[c][b] = 0;
                    if (rule2(D, a, b, c))
                        D[c][a] = 0;

                    // If we're able to, iterate quadruplets.
                    if (n > 3) {
                        for (size_t d = 0; d < n; d++) {
                            if (d == a || d == b || d == c) continue;

                            if (rule3(D, a, b, c, d))
                                D[d][b] = 0;
                            if (rule4(D, a, b, c, d))
                                D[b][a] = D[b][c] = 0;
                        } //end for
                    } // end if
                } // end for
            } // end for
        } // end for

        // Check for closure.
        if (closed(D, G)) {
            bool fullyDirected = true;

            // Check for further undirected edges.
            for (size_t a = 0; a < n; a++) {
                for (size_t b = 0; b < n; b++) {
                    if (D[a][b] && D[b][a]) {
                        D[b][a] = 0;
                        lastDirected = make_pair(b, a);
                        
                        fullyDirected = false;
                        C.push(make_pair(lastDirected, D));
                        break;
                    } // end if
                } // end for

                if (!fullyDirected) break;
            } // end for

            if (fullyDirected)
                return D;
        } else {
            // If the stack is empty, then this PDAG cannot be extended.
            pair<pair<size_t, size_t>, vector<vector<size_t>>> prev;
            try {
                prev = C.top();
                C.pop();
            } catch (const exception& e) {
                throw 21;
            } // end try-catch

            lastDirected = prev.first;
            D = prev.second;

            D[lastDirected.first][lastDirected.second] = 1;
            D[lastDirected.second][lastDirected.first] = 0;
        } // end if-else
    } // end while
} // end function

// Recursively walk to descendents and check if in L.
const bool DAG::active(const set<size_t>& L, const size_t j) {
    // Check for node j in L.
    for (size_t l : L)
        if (j == l) return true;
    
    // Otherwise, check descendents.
    bool isActive = false;
    for (size_t k = 0; k < adjacencyMatrix[j].size(); k++)
        if (adjacencyMatrix[j][k])
            isActive = isActive || active(L, k);

    return isActive;
} // end function

// Get legal paths as ground-truth matrix in the sense of Algorithm 2 (iii).
const vector<vector<bool>> DAG::legal(const set<size_t>& L, const vector<bool>& descendent) {
    vector<vector<bool>> F(adjacencyMatrix.size(), vector<bool>(adjacencyMatrix.size(), false));

    for (size_t i = 0; i < adjacencyMatrix.size(); i++) {
        for (size_t j = 0; j < adjacencyMatrix[i].size(); j++) {
            if (!F[i][j] && adjacencyMatrix[i][j]) {
                bool inL = (L.find(adjacencyMatrix[i][j]) != L.end());
                F[i][j] = !inL;
            } // end if
        } // end for
    } // end for

    return F;
} // end function

// Algorithm 1
const map<size_t, bool> DAG::reachable(const set<size_t>& L, const set<size_t>& J, const vector<bool>& descendent) {
    size_t n = adjacencyMatrix.size();

    // Undirect D.
    vector<vector<size_t>> Dprime = adjacencyMatrix;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (Dprime[i][j]) Dprime[j][i] = 1;
    
    map<size_t, bool> R;
    queue<pair<size_t, size_t>> links;

    map<pair<size_t, size_t>, size_t> label;
    size_t i = 1;

    // Fill R with falses.
    for (size_t i = 0; i < n; i++)
        R[i] = false;

    // Set new node s to true.
    R[n] = true;

    // Set all nodes in J to true and push s -> j edges to queue.
    for (size_t j : J) {
        R[j] = true;
        links.push(make_pair(n, j));
        label[make_pair(n, j)] = i;
    } // end for

    // Start walking s -> j, looking for legal paths j -> k.
    // If any found, label R[k] = true, then walk to their descendents.
    while (!links.empty()) {
        pair<size_t, size_t> edge = links.front();

        size_t u = edge.first;
        size_t v = edge.second;

        links.pop();

        vector<pair<size_t, size_t>> candidates;
        for (size_t w = 0; w < n; w++) {
            // Ensure u != w.
            if (u == w) continue;

            if (label[edge] == i && Dprime[v][w]) {
                // Check if head-to-head in original graph. Special case: s -> v <- w is not head-to-head.
                bool headToHead;
                if (u != n)
                    headToHead = (adjacencyMatrix[u][v] && adjacencyMatrix[v][w]);
                else
                    headToHead = false;

                bool inL = (L.find(w) != L.end());

                if ((headToHead && descendent[v]) || (!headToHead && !inL))
                    candidates.push_back(make_pair(v, w));
            } // end if
        } // end for

        for (pair<size_t, size_t> unlabeled : candidates) {
            label[unlabeled] = i + 1;
            R[unlabeled.second] = true;
        } // end for
        i++;
    } // end while

    return R;
} // end function

// DAG Methods

DAG::DAG() {
    adjacencyMatrix = vector<vector<size_t>>();
} // end function

// Constructor given a CI structure and the number of vertices.
DAG::DAG(const CI& M, const size_t n) {
    vector<vector<size_t>> G;
    try {
        G = generatePDAG(M, n);
    } catch (int e) {
        if (e == 20)
            cout << "Error: PDAG generation failed. This CI structure is not DAG-representable." << endl;
        throw 20;
    } // end try-catch

    vector<vector<size_t>> D;
    try {
        D = extendToDAG(G);
    } catch (int e) {
        if (e == 21)
            cout << "Error: DAG extension failed. This CI structure is not DAG-representable." << endl;
        throw 21;
    } // end try-catch

    adjacencyMatrix = D;

    if (!consistent(M))
        throw 22;
} // end function

const vector<vector<size_t>> DAG::matrix() {
    return adjacencyMatrix;
} // end function

// Generate topological ordering via Kahn (1962).
vector<size_t> DAG::topologicalOrder() {
    vector<size_t> L;
    queue<size_t> S;

    vector<vector<size_t>> copy = adjacencyMatrix;
    size_t d = adjacencyMatrix.size();

    // Populate S with nodes having no parents.
    for (size_t a = 0; a < d; a++) {
        bool parentless = true;
        for (size_t b = 0; b < d; b++)
            if (copy[b][a])
                parentless = false;

        if (parentless)
            S.push(a);
    } // end for

    while (!S.empty()) {
        size_t n = S.front();
        S.pop();

        L.push_back(n);
        
        for (size_t m = 0; m < d; m++) {
            if (copy[n][m]) {
                copy[n][m] = 0;

                bool incoming = false;
                for (size_t l = 0; l < d; l++)
                    if (copy[l][m])
                        incoming = true;

                if (!incoming)
                    S.push(m);
            } // end if
        } // end for
    } // end while

    // Make sure the graph has no edges (D is 0 matrix).
    for (size_t n = 0; n < d; n++)
        for (size_t m = 0; m < d; m++)
            if (copy[n][m])
                throw 23;

    return L;
} // end function

// Phase 3
const bool DAG::consistent(const CI& M) {
    // Get CI structure of D, compare with M.
    size_t n = adjacencyMatrix.size();
    
    set<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.insert(i);
    set<set<size_t>> PN = P(N, 2);

    CI CIStruct = toCI();
    if (CIStruct != M)
        return false;

    // Run second ordering test.
    vector<size_t> U;
    try {
        U = topologicalOrder();
    } catch (int e) {
        if (e == 23)
            cout << "Error: Unable to generate topological ordering of graph." << endl;
    } // end try-catch
    
    for (size_t a = 0; a < n; a++) {
        // Create sets of nodes that (a) precede a in the ordering or (b) are parents of a.
        set<size_t> Ua;
        for (size_t i : U) {
            if (i == a) break;
            Ua.insert(i);
        } // end for

        set<size_t> pi;
        for (size_t b = 0; b < n; b++)
            if (adjacencyMatrix[b][a])
                pi.insert(b);

        // Complement Ua wrt pi.
        set<size_t> C;
        set_difference(Ua.begin(), Ua.end(), pi.begin(), pi.end(), inserter(C, C.end()));

        // Check for CI in complemented Ua, given pi.
        for (size_t b : C) {
            set<size_t> pair;
            pair.insert(a);
            pair.insert(b);

            if (M.structure.find(pair) == M.structure.end())
                return false;
        } // end for
    } // end for

    // We haven't failed either test, so D and M must be consistent.
    return true;
} // end function

const set<size_t> DAG::maximalIndependentSet(const set<size_t>& J, const set<size_t>& L) {
    size_t n = adjacencyMatrix.size();
    
    vector<bool> descendent(n, false);
    for (size_t i = 0; i < n; i++)
        descendent[i] = active(L, i);

    vector<vector<bool>> F = legal(L, descendent);
    map<size_t, bool> R = reachable(L, J, descendent);

    set<size_t> Kprime;
    for (size_t i = 0; i < n; i++)
        if (R[i]) Kprime.insert(i);

    for (size_t j : J)
        Kprime.insert(j);
    for (size_t l : L)
        Kprime.insert(l);

    // Get the ground set of vertices (nodes).
    set<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.insert(i);

    // Return complement of N wrt Kprime.
    set<size_t> K;
    set_difference(N.begin(), N.end(), Kprime.begin(), Kprime.end(), inserter(K, K.end()));

    return K;
} // end function

// Algorithm 2
const CI DAG::toCI() {
    map<set<size_t>, set<set<size_t>>> CIStruct;
    size_t n = adjacencyMatrix.size();

    set<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.insert(i);
    set<set<size_t>> PN = P(N, 2);

    for (size_t j : N) {
        set<size_t> J = {j};

        set<size_t> JComp;
        set_difference(N.begin(), N.end(), J.begin(), J.end(), inserter(JComp, JComp.end()));

        set<set<size_t>> PL = P(JComp);
        for (set<size_t> L : PL) {
            set<size_t> K = maximalIndependentSet(J, L);

            if (!K.empty()) {
                for (size_t k : K) {
                    set<size_t> I = {j, k};

                    if (CIStruct.find(I) == CIStruct.end())
                        CIStruct[I] = {L};
                    else
                        CIStruct[I].insert(L);
                } // end for
            } // end if
        } // end for
    } // end for

    return CI(CIStruct);
} // end function

const bool stronglyProtected(const vector<vector<size_t>>& G, const size_t a, const size_t b) {
    size_t n = G.size();
    
    for (size_t c = 0; c < n; c++) {
        if (c == a || c == b) continue;

        if (!G[a][c] && G[c][a] && !G[b][c] && !G[c][b])
            return true;

        if (!G[a][c] && !G[c][a] && !G[b][c] && G[c][b])
            return true;

        if (G[a][c] && !G[c][a] && !G[b][c] && G[c][b])
            return true;

        if (n > 3) {
            for (size_t d = 0; d < n; d++) {
                if (d == a || d == b || d == c) continue;

                if (!G[a][c] && !G[c][a] && !G[b][c] && G[c][b] && !G[a][d] && !G[d][a] && !G[b][d] && G[d][b])
                    return true;
            } // end for
        } // end if
    } // end for

    return false;
} // end function

const vector<vector<size_t>> DAG::essentialGraph() {
    vector<vector<size_t>> nextG = adjacencyMatrix;
    vector<vector<size_t>> currG;

    size_t n = adjacencyMatrix.size();

    while (nextG != currG) {
        currG = nextG;

        for (size_t a = 0; a < n; a++)
            for (size_t b = 0; b < n; b++)
                if (currG[a][b] && !currG[b][a] && !stronglyProtected(currG, a, b))
                    nextG[b][a] = 1;
    } // end while

    return currG;
} // end function

#endif