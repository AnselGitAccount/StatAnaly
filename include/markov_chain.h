#ifndef MARKOV_CHAIN_SOLVER_H_
#define MARDOV_CHAIN_SOLVER_H_

#include "graph.h"
#include "hasher.h"
#include <unordered_map>
#include <map>


namespace statanaly {


/**
 * @brief A solver for markov chain type of problems.
 * 
 * Approach:
 *  - DFS with backtracking.
 *  - Recursive.
 *  - Memoization: Save the subgraph that has already been solved.
 * 
 * @param knockouts List of nodes that are unavailable during traversal at a particular time.
 * @param memory Record of computed probability of pathes. 
 * @param terms Terminal nodes.
 */
class Markovchain {
    using node = std::size_t;
    std::set<node> knockouts;
    std::map<std::tuple<node,node,std::set<node>>, double> memory;
    std::vector<node> terms;

    // Compute probability of a loop from i to f WITHOUT revisiting i,
    // on a subgraph.
    double find_prob_i2f_wo_revisit(node i, node f) {
        // Transition probability from i.
        auto& tran_i = graph.adjmat()[i];

        // Search in memory.
        auto it = memory.find( std::tuple<node,node,std::set<node>>{i,f,knockouts} );
        if (it != memory.end()) {
            return it->second;
        }

        // Copy i's outgoing edges.
        auto outgoings = tran_i;
        
        // Knockout i's outgoing edges.
        std::fill(tran_i.begin(), tran_i.end(), 0);

        // Recursively compute the probability.
        double pcif = 0.;
        for (node s=0; s<outgoings.size(); s++) {
            pcif += outgoings[s] * find_prob_i2f(s, f);
        }

        // Backtrack -- put i's outgoing edges back.
        tran_i = outgoings;
        knockouts.erase(i);

        // Memorize the transition probability of this path.
        std::pair kv(std::tuple<node,node,std::set<node>>{i,f,knockouts}, pcif);
        memory.insert( kv );

        return pcif;
    }

public:
    DWGraph<double> graph;

    Markovchain(const DWGraph<double>& g) : graph(g){}
    
    void resetGraph(const DWGraph<double>& g) {
        graph = g;
    }

    /**
     * @brief Entry Function -- find transit probability from i to all terminal states.
     * 
     * @param i Starting node.
     */
    auto find_transit_prob (node i) {
        std::map<node, double> prob;
        find_terminal_states();

        knockouts.clear();
        memory.clear();

        for (auto t : terms) {
            prob[t] = find_prob_i2f(i, t);
        }

        return prob;
    }


    /**
     * @brief Compute probability from node i to node f recursively.
     * 
     * Computation terminates when no outgoing edges OR a terminal node.
     * 
     * @param i 
     * @param f 
     * @return double 
     */
    double find_prob_i2f (node i, node f) {
        const auto& row = graph.adjmat()[i];
        bool terminate = std::all_of(row.begin(), row.end(), [](auto e) {return e==0;});
        terminate |= (std::find(terms.begin(), terms.end(), i) != terms.end());

        if (terminate) {
            if (i==f) {return 1;} 
            else {return 0;}
        }

        // Compute probability from i to f WITHOUT revisiting i.
        auto pcif = find_prob_i2f_wo_revisit(i, f);

        // Compute probabiilty of a loop from i to itself WITHOUT revisiting i.
        auto pcii = find_prob_i2f_wo_revisit(i, i);

        return pcif / (1. - pcii);
    }


    /**
     * @brief Find terminal states.
     * 
     * Two types of terminal states:
     * 1. Ending -- all zeros in a row.
     * 2. Absorbing -- one in the diagnal, and zeros elsewhere in a row.
     * 
     * @return std::vector<node> 
     */
    std::vector<node> find_terminal_states() {
        const auto& adjmat = graph.adjmat();
        const std::size_t len = graph.adjmat().len();
        terms.clear();
        
        for (std::size_t i=0; i<adjmat.len(); i++) {
            const auto& row = adjmat[i];
            double zeros[len] = {};
            bool ending = memcmp(row.data(), zeros, sizeof(double)*len)==0;

            zeros[i] = double(1);
            bool absorbing = memcmp(row.data(), zeros, sizeof(double)*len)==0;
            
            if (ending || absorbing)
                terms.push_back(i);
        }

        return terms;
    }

};


};

#endif  