#ifndef STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_
#define STATANALY_DIRECTED_UNWEIGHTED_GRAPH_H_

#include "adjacency.h"
#include <iostream>
#include <type_traits>

/**
 * @file graph.h
 * 
 * Undirected Unweighted graph can be implemented as 
 * a Directed graph with bidirectional connectivity.
 * 
 */


namespace statanaly {

/**
 * @brief Base class for directed Weighted or Unweighted graph.
 * 
 * The graph is represented with Adjacency matrix.
 * 
 * @tparam DERI Type of graph -- Weighted or Unweighted.
 * @tparam S Type of element in Adjacency matrix.
 */
template<class DERI, class S>
class Graph {
protected:
    // Default to using Adjacency matrix.
    Adjmat<S> mat;

public:
    // Expose interface functions here.
    constexpr const Adjmat<S>& adjmat() const {return mat;}
    constexpr       Adjmat<S>& adjmat()       {return mat;}

};


/**
 * @brief Directed Weighted Graph
 * 
 * @tparam S Type of element in Adjacency matrix.
 */
template<class S>
class DWGraph : public Graph<DWGraph<S>, S> {
public:    
    DWGraph (const Adjmat<S>& m) {
        this->mat = m;
    }
};


/**
 * @brief Directed Unweighted graph.
 * 
 * @tparam S Type of element in Adjacency matrix.
 */
template<class S>
class DUGraph : public Graph<DUGraph<S>, S> {
public:
    DUGraph (const Adjmat<S>& m) {
        this->mat = m;
    }
};


/**
 * @brief Two Directed Weighted graph objects are equal if their adjacency matrixes are equal.
 * 
 * @param a Graph A.
 * @param b Graph B.
 * @return true 
 * @return false 
 */
template<class S>
inline bool operator == (const DWGraph<S>& a, const DWGraph<S>& b) {
    return a.adjmat() == b.adjmat();
};



/**
 * @brief Two Directed Unweighted graph objects are equal if their adjacency matrixes are equal.
 * 
 * @param a Graph A.
 * @param b Graph B.
 * @return true 
 * @return false 
 */
template<class S>
inline bool operator == (const DUGraph<S>& a, const DUGraph<S>& b) {
    return a.adjmat() == b.adjmat();
};

}  // namespace


#endif