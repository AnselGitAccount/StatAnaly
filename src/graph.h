#ifndef STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_
#define STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_

#include "adjacency.h"
#include <iostream>
#include <type_traits>

namespace statanaly {

template<class DERI>
class Graph {
public:
    // Expose interface functions here.
};

/* Directed Weighted graph */

template<class S>
class DWGraph : public Graph<DWGraph<S>> {
public:
    // Default to using Adjacentcy matrix.
    Adjmat<S> adjmat;

    // Constructor
    DWGraph (const Adjmat<S>& mat) { 
        adjmat.clone(mat);
    }

    // Copy operations -- clone the adjmat.
    DWGraph (const DWGraph<S>& other) {
        this->operator=(other);
    }
    DWGraph& operator = (const DWGraph<S>& other) {
        adjmat.clone(other.adjmat);
        return *this;
    }

    // Move operations -- reassign the ownership of the adjmat.
    DWGraph (DWGraph<S>&& other) = default;
    DWGraph& operator = (DWGraph<S>&& other) = default;

};


/* Directed Unweighted graph.
 * Undirected Unweighted graph is implemented as 
 * a Directed graph with bidirectional connectivity.
 */

template<class S>
class DUGraph : public Graph<DUGraph<S>> {
public:
    // Default to using Adjacentcy matrix.
    Adjmat<S> adjmat;

    // Constructor
    DUGraph (const Adjmat<S>& mat) {
        adjmat.clone(mat);
    }


    // Copy operations -- clone the adjmat.
    DUGraph (const DUGraph<S>& other) {
        this->operator=(other);
    };
    DUGraph& operator = (const DUGraph<S>& other) {
        adjmat.clone(other.adjmat);
        return *this;
    };

    // Move operations.
    DUGraph (DUGraph&& other) = default;
    DUGraph& operator = (DUGraph&& other) = default;

};


}  // namespace


#endif