#ifndef STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_
#define STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_

#include "adjacency.h"
#include <iostream>
#include <type_traits>

namespace statanaly {

template<class DERI, class S>
class Graph {
protected:
    // Default to using Adjacentcy matrix.
    Adjmat<S> mat;

public:
    // Expose interface functions here.
    constexpr const Adjmat<S>& adjmat() const {return mat;}
    constexpr       Adjmat<S>& adjmat()       {return mat;}

};

/* Directed Weighted graph */

template<class S>
class DWGraph : public Graph<DWGraph<S>, S> {
public:    
    DWGraph (const Adjmat<S>& m) {
        this->mat = m;
    }
};


/* Directed Unweighted graph.
 * Undirected Unweighted graph is implemented as 
 * a Directed graph with bidirectional connectivity.
 */

template<class S>
class DUGraph : public Graph<DUGraph<S>, S> {
public:
    DUGraph (const Adjmat<S>& m) {
        this->mat = m;
    }
};


template<class S>
inline bool operator == (const DWGraph<S>& a, const DWGraph<S>& b) {
    return a.adjmat() == b.adjmat();
};

template<class S>
inline bool operator == (const DUGraph<S>& a, const DUGraph<S>& b) {
    return a.adjmat() == b.adjmat();
};

}  // namespace


#endif