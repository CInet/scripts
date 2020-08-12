#include <iostream>
#include <chrono>
#include <set>
#include <map>
#include <algorithm>

#include "util.h"
#include "ci.h"
#include "dag.h"

using namespace std;
using namespace chrono;

int main() {
    auto start = high_resolution_clock::now();

    size_t n = 3;
    vector<vector<vector<size_t>>> dags;
    set<vector<vector<size_t>>> essentials;

    set<size_t> N;
    for (size_t i = 0; i < n; i++)
        N.insert(i);

    // map<set<size_t>, set<set<size_t>>> structure;
    // structure[{0, 1}] = {{}, {2}};
    // structure[{0, 2}] = {{}, {1}};

    // CI M(structure);
    // try {
    //     DAG D(M, n);
    //     vector<vector<size_t>> graph = D.essentialGraph();
    // } catch (int e) {
    //     cout << "thrown: " << e << endl;
    // }

    // Enumerate all ijk triplets.
    // Doesn't necessarily enumerate all semigraphoids.
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;

            for (size_t k = 0; k < n; k++) {
                if (k == i || k == j) continue;

                set<size_t> ijk = {i, j, k};

                // Enumerate all combinations of L \subset N\ijk in structure.
                set<size_t> NCompIJK;
                set_difference(N.begin(), N.end(), ijk.begin(), ijk.end(), inserter(NCompIJK, NCompIJK.end()));

                set<set<set<size_t>>> PNCompIJK = P(P(NCompIJK));
                for (set<set<size_t>> PL : PNCompIJK) {
                    map<set<size_t>, set<set<size_t>>> structure;
                    
                    for (set<size_t> L : PL) {
                        structure[{i, j}].insert(L);
                        structure[{i, k}].insert(L);

                        set<size_t> jL = L;
                        jL.insert(j);
                        structure[{i, k}].insert(jL);

                        set<size_t> kL = L;
                        kL.insert(k);
                        structure[{i, j}].insert(kL);
                    } // end for

                    CI M(structure);
                    try {
                        DAG D(M, n);
                        dags.push_back(D.matrix());
                        essentials.insert(D.essentialGraph());
                    } catch (const int e) {}
                } // end for
            } // end for
        } // end for
    } // end for

    auto elapsed = duration_cast<microseconds>(high_resolution_clock::now() - start);

    cout << dags.size() << " DAG-representable semigraphoids found of size " << n << "." << endl;
    cout << essentials.size() << " equivalence classes found of size " << n << "." << endl;
    
    // Print out the DAGs.
    cout << endl;
    for (vector<vector<size_t>> D : essentials) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++)
                cout << D[i][j] << " ";
            cout << endl;
        } // end for
        
        cout << endl;
    } // end for

    cout << "Time Elapsed: " << elapsed.count() << "s" << endl;
    return 0;
} // end main