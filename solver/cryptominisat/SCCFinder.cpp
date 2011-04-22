/*****************************************************************************
CryptoMiniSat -- Copyright (c) 2009 Mate Soos

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#include <iostream>
#include <vector>
#include "SolverTypes.h"
#include "SCCFinder.h"
#include "VarReplacer.h"
#include <iomanip>
#include "time_mem.h"
#include "Subsumer.h"
#include "XorSubsumer.h"

SCCFinder::SCCFinder(Solver& _solver) :
    solver(_solver)
    , varElimed1(_solver.subsumer->getVarElimed())
    , varElimed2(_solver.xorSubsumer->getVarElimed())
    , replaceTable(_solver.varReplacer->getReplaceTable())
    , totalTime(0.0)
{}

const bool SCCFinder::find2LongXors()
{
    double myTime = cpuTime();
    uint32_t oldNumReplace = solver.varReplacer->getNewToReplaceVars();

    globalIndex = 0;
    index.clear();
    index.resize(solver.nVars()*2, std::numeric_limits<uint32_t>::max());
    lowlink.clear();
    lowlink.resize(solver.nVars()*2, std::numeric_limits<uint32_t>::max());
    stackIndicator.clear();
    stackIndicator.growTo(solver.nVars()*2, false);
    assert(stack.empty());

    for (uint32_t vertex = 0; vertex < solver.nVars()*2; vertex++) {
        //Start a DFS at each node we haven't visited yet
        if (index[vertex] == std::numeric_limits<uint32_t>::max()) {
            recurDepth = 0;
            tarjan(vertex);
            assert(stack.empty());
        }
    }

    if (solver.conf.verbosity >= 2 || (solver.conflicts == 0 && solver.conf.verbosity  >= 1)) {
        std::cout << "c Finding binary XORs  T: "
        << std::fixed << std::setprecision(2) << std::setw(8) <<  (cpuTime() - myTime) << " s"
        << "  found: " << std::setw(7) << solver.varReplacer->getNewToReplaceVars() - oldNumReplace
        << std::endl;
    }
    totalTime += (cpuTime() - myTime);

    return solver.ok;
}

void SCCFinder::tarjan(const uint32_t vertex)
{
    recurDepth++;
    index[vertex] = globalIndex;  // Set the depth index for v
    lowlink[vertex] = globalIndex;
    globalIndex++;
    stack.push(vertex); // Push v on the stack
    stackIndicator[vertex] = true;

    Var vertexVar = Lit::toLit(vertex).var();
    if (!varElimed1[vertexVar] && !varElimed2[vertexVar]) {
        const vec2<Watched>& ws = solver.watches[vertex];
        for (vec2<Watched>::const_iterator it = ws.getData(), end = ws.getDataEnd(); it != end; it++) {
            if (!it->isBinary()) continue;
            const Lit lit = it->getOtherLit();

            doit(lit, vertex);
        }

        if (solver.conf.doExtendedSCC) {
            Lit vertLit = Lit::toLit(vertex);
            vector<Lit>& transCache = solver.transOTFCache[(~Lit::toLit(vertex)).toInt()].lits;
            vector<Lit>::iterator it = transCache.begin();
            vector<Lit>::iterator it2 = it;
            uint32_t newSize = 0;
            for (vector<Lit>::iterator end = transCache.end(); it != end; it++) {
                Lit lit = *it;
                lit = replaceTable[lit.var()] ^ lit.sign();
                if (lit == vertLit || varElimed1[lit.var()] || varElimed2[lit.var()]) continue;
                *it2++ = lit;
                newSize++;

                doit(lit, vertex);
            }
            transCache.resize(newSize);
        }
    }

    // Is v the root of an SCC?
    if (lowlink[vertex] == index[vertex]) {
        uint32_t vprime;
        tmp.clear();
        do {
            assert(!stack.empty());
            vprime = stack.top();
            stack.pop();
            stackIndicator[vprime] = false;
            tmp.push(vprime);
        } while (vprime != vertex);
        if (tmp.size() >= 2) {
            for (uint32_t i = 1; i < tmp.size(); i++) {
                if (!solver.ok) break;
                vec<Lit> lits(2);
                lits[0] = Lit::toLit(tmp[0]).unsign();
                lits[1] = Lit::toLit(tmp[i]).unsign();
                const bool xorEqualsFalse = Lit::toLit(tmp[0]).sign()
                                            ^ Lit::toLit(tmp[i]).sign()
                                            ^ true;
                if (solver.value(lits[0]) == l_Undef && solver.value(lits[1]) == l_Undef) {
                    solver.varReplacer->replace(lits, xorEqualsFalse, 0, true, false);
                }
            }
        }
    }
}
