#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
using namespace std;

// Recursively walk to descendents and check if in L.
const bool active(const vector<vector<size_t>>& D, const vector<size_t>& L, size_t j) {
    // Check for node j in L.
    for (size_t i = 0; i < L.size(); i++)
        if (j == L[i]) return true;
    
    // Otherwise, check descendents.
    bool isActive = false;
    for (size_t k = 0; k < D[j].size(); k++)
        if (D[j][k])
            isActive = isActive || active(D, L, k);

    return isActive;
} // end function

// Get legal paths as ground-truth matrix in the sense of Algorithm 2 (iii).
const vector<vector<bool>> legal(const vector<vector<size_t>>& D, const vector<size_t>& J, const vector<bool>& descendent) {
    vector<vector<bool>> F(D.size(), vector<bool>(D.size(), false));

    for (size_t i = 0; i < D.size(); i++) {
        for (size_t j = 0; j < D[i].size(); j++) {
            if (!F[i][j] && D[i][j]) {
                bool inJ = false;
                for (size_t k = 0; k < J.size(); k++)
                    inJ = (D[i][j] == J[k]);

                F[i][j] = F[i][j] || descendent[j] || !inJ;
            } // end if
        } // end for
    } // end for

    return F;
} // end function

// Algorithm 1
const map<size_t, bool> reachable(const vector<vector<size_t>>& D, const vector<vector<bool>>& F, const vector<size_t>& J) {
    map<size_t, bool> R;
    // Fill R with falses.
    for (size_t i = 0; i < D.size(); i++)
        R[i] = false;

    // Set new node s to true.
    R[D.size()] = true;

    // Set all nodes in J to true.
    for (auto j = J.begin(); j != J.end(); j++)
        R[*j] = true;

    for (size_t i = 0; i < D.size(); i++)
        for (size_t j = 0; j < D[i].size(); j++)
            if (D[i][j] && !R[j])
                R[j] = F[i][j];

    return R;
} // end function

// Algorithm 2
const vector<size_t> ci(const vector<vector<size_t>>& D, const vector<size_t>& J, const vector<size_t>& L) {
    size_t n = D.size();
    
    vector<bool> descendent(n, false);
    for (size_t i = 0; i < n; i++)
        descendent[i] = active(D, L, i);

    // Undirect D.
    vector<vector<size_t>> Dprime = D;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            if (Dprime[i][j])
                Dprime[j][i] = 1;

    vector<vector<bool>> F = legal(D, J, descendent);
    map<size_t, bool> R = reachable(D, F, J);

    set<size_t> Kprime;
    for (size_t i = 0; i < n; i++)
        if (R[i]) Kprime.insert(i);

    for (size_t j : J)
        Kprime.insert(j);
    for (size_t l : L)
        Kprime.insert(l);

    // Get the ground set of vertices (nodes).
    vector<size_t> K;
    for (size_t i = 0; i < n; i++)
        K.push_back(i);

    // Remove Kprime.
    for (size_t k : Kprime)
        K.erase(K.begin() + k - n + K.size());

    return K;
} // end function

// Check if directed, simple graph is acyclic.
const bool acyclic(const vector<vector<size_t>>& D, size_t j, size_t steps) {
    if (steps >= D.size()) return false;

    bool isAcyclic = true;
    
    for(size_t k = 0; k < D[j].size(); k++)
        if (D[j][k])
            isAcyclic = isAcyclic && acyclic(D, k, steps + 1);

    return isAcyclic;
} // end function

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
const vector<vector<T>> P(const vector<T>& N, const int r=-1) {
    size_t n = N.size();
    size_t p = pow(2, n);
    
    vector<vector<T>> powerSet;
    for (size_t i = 0; i < p; i++) {
        if (r >= 0 && setBits(i) != r)
            continue;

        vector<T> subset;
        for (size_t j = 0; j < n; j++)
            if (i & (1 << j))
                subset.push_back(N[j]);
        
        powerSet.push_back(subset);
    } // end for

    return powerSet;
} // end function

int main() {
    size_t n = 3;

    for (size_t E = 0; E < pow(2, n * (n - 1)); E++) {
        // Skip non-simple graphs.
        if (setBits(E) > n - 1) continue;

        vector<vector<size_t>> D;
        for (size_t i = 0; i < n; i++) {
            vector<size_t> adjList;
            for (size_t j = 0; j < n - 1; j++) {
                if (E & (1 << (j * n + i)))
                    adjList.push_back(1);
                else
                    adjList.push_back(0);
            } // end for

            adjList.insert(adjList.begin() + i, 0);
            D.push_back(adjList);
        } // end for

        // Skip cyclic graphs.
        bool isAcyclic = true;
        for (size_t i = 0; i < D.size(); i++)
            isAcyclic = isAcyclic && acyclic(D, i, 0);
        if (!isAcyclic) continue;

        for (size_t i = 0; i < D.size(); i++) {
            for (size_t j = 0; j < D[i].size(); j++)
                cout << D[i][j] << " ";
            cout << endl;
        } // end for
        cout << endl;

        vector<size_t> N;
        for (size_t i = 0; i < n; i++)
            N.push_back(i);
        vector<vector<size_t>> PN = P(N);

        // Iterate over all possible separating sets L and predicate sets J.
        for (auto L = PN.begin(); L != PN.end(); L++) {
            vector<size_t> Lc = N;
            for (auto node = L->begin(); node != L->end(); node++)
                Lc.erase(Lc.begin() + (*node) - n + Lc.size());
            vector<vector<size_t>> PLc = P(Lc);

            for (auto J = PLc.begin(); J != PLc.end(); J++) {
                vector<size_t> Kmax = ci(D, *J, *L);
                vector<vector<size_t>> PKmax = P(Kmax);

                // Start outputing CI statements.
                cout << "J = {";
                for (size_t j : *J)
                    cout << j << ", ";
                if (J->size()) cout << "\b \b\b \b}" << endl;
                else cout << "}" << endl;

                cout << "L = {";
                for (size_t l : *L)
                    cout << l << ", ";
                if (L->size()) cout << "\b \b\b \b}" << endl;
                else cout << "}" << endl;

                cout << "K ∈  {";
                // Iterate all CI sets K.
                for (auto K = PKmax.begin(); K != PKmax.end(); K++) {
                    cout << "\n\t{";
                    for (size_t k : *K)
                        cout << k << ", ";
                    if (K->size()) cout << "\b \b\b \b}" << endl;
                    else cout << "},";
                } // end for

                cout << "\b \b\n}\n" << endl;
            } // end for
        } // end for

        cout << endl;
    } // end for

    return 0;
} // end main