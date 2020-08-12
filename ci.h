#ifndef CI_H
#define CI_H

#include <set>
#include <map>
using namespace std;

class CI {
    public:
        // TODO: Make private?
        map<set<size_t>, set<set<size_t>>> structure;

        CI();
        CI(const map<set<size_t>, set<set<size_t>>>& M);

        const bool operator==(const CI& M);
        const bool operator!=(const CI& M);

        const bool isSemigraphoid();
}; // end class

CI::CI() {
    structure = map<set<size_t>, set<set<size_t>>>();
} // end function

CI::CI(const map<set<size_t>, set<set<size_t>>>& M) {
    // Verify all keys are pairs of variables (sets of size 2).
    for (pair<set<size_t>, set<set<size_t>>> statement : M)
        if (statement.first.size() != 2)
            throw 20;

    // Otherwise, succeed.
    structure = M;
} // end function

const bool CI::operator==(const CI& M) {
    return structure == M.structure;
} // end function

const bool CI::operator!=(const CI& M) {
    return structure != M.structure;
} // end function

const bool CI::isSemigraphoid() {
    // Check local semigraphoid axiom.
    bool semigraphoid = true;
    for (auto it = structure.begin(); it != structure.end(); it++) {
        return true; // TODO: finish
    } // end for
} // end function

#endif