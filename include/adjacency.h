/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef STATANALY_ADJACENCY_H_
#define STATANALY_ADJACENCY_H_

#include <array>
#include <string>
#include <vector>
#include <numeric>

namespace statanaly {


template<class T>
class Node {
    T val;
    std::vector<Node<T>*> nei;

public:
    Node() : 
        val(0), 
        nei(std::vector<Node<T>*>{}) {}

    Node(T val_) : 
        val(val_), 
        nei(std::vector<Node<T>*>{}) {}

    Node(T val_, std::vector<Node<T>*> nei_) : 
        val(val_),
        nei(nei_) {};

    // Can't copy construct.
    // Can't copy assign.
    Node(const Node<T>& other) = delete;
    Node& operator = (const Node<T>& other) = delete;

    std::vector<Node<T>*> getNei() const {
        return nei;
    };

    T getVal() const {
        return val;
    }
};


/**
 * @brief Adjacency Matrix.
 * 
 * Diagonal terms are deafult to zero.
 * 
 * @tparam SCALAR Matrix element type.
 */
template<class SCALAR = int32_t>
class Adjmat {
    std::vector<std::vector<SCALAR>> d;

public:
    Adjmat() = default;
    
    Adjmat(std::size_t n) {
        d.resize(n);
        for(auto& row : d) {
            row.resize(n);
        }
    }

    inline       std::vector<SCALAR>& operator[] (std::size_t idx)       { return d[idx]; }
    inline const std::vector<SCALAR>& operator[] (std::size_t idx) const { return d[idx]; }

    constexpr const std::size_t len() const { return d.size(); }
    constexpr const auto& data() const { return d; }
    constexpr       auto& data()       { return d; }
};

template<class S>
inline bool operator == (const Adjmat<S>& a, const Adjmat<S>& b) {
    return a.data() == b.data();
};


/**
 * @brief Normalize Each Row of an Adjacency Matrix.
 * 
 * Normalize row elements such that each row sums up to one.
 * 
 * @tparam S Template parameter of Adjmat.
 * @param mat Adjacency matrix.
 * @return Adjmat<double> 
 */
template<class S>
Adjmat<double> normalizeRow(Adjmat<S>& mat) {
    Adjmat<double> norm;
    for(const auto& row : mat.data()) {
        std::vector<double> temp(row.begin(), row.end());
        double rowsum = std::reduce(temp.begin(), temp.end());
        if (rowsum!=0)
            std::for_each(temp.begin(), temp.end(), [rowsum](double& e){e/=rowsum;});

        norm.data().emplace_back( temp.begin(), temp.end() );
    }
    return norm;
};


/**
 * @brief Convert ASCII Text to an Adjacency Matrix.
 * 
 * Convert ASCII Text to an Adjacency matrix for an undirected graph.
 * For example:
 * 
 * "AB1,AC1,CD1" is translated to 
 * 0 1 1 0 
 * 1 0 0 0 
 * 1 0 0 1 
 * 0 0 1 0 
 * 
 * @tparam T Template paramter of Adjmat.
 * @param str ASCII text that descripts connectivity.
 * @return Adjmat<T> 
 */
template<class T>
Adjmat<T> convert2Adjmat(std::string str) {
    // trim white space
    str.erase( std::remove(str.begin(), str.end(), ' '), str.end() );

    // seperate into chunck using delimiter ','.
    std::vector<std::string> vs;
    for(auto s_itr=str.begin(), d_itr=s_itr; 
            d_itr!=str.end(); 
            s_itr=d_itr+1) {
        d_itr = std::find(d_itr+1, str.end(), ',');
        vs.emplace_back(std::string(s_itr, d_itr));
    };

    // Find number of nodes
    int min=255, max=0; 
    for (const std::string& s : vs) {
        // Assume node labels are single character long.
        min = min < int(s[0]) ? min : int(s[0]);
        min = min < int(s[1]) ? min : int(s[1]);
        max = int(s[0]) < max ? max : int(s[0]);
        max = int(s[1]) < max ? max : int(s[1]);
    }
    int n_nodes = max - min + 1;

    // Build Adjacency matrix
    Adjmat<T> AdjMat(n_nodes);
    for(const std::string& s : vs) {
        // ASCII table offset between 'A' to 0 is 65.

        int idx1 = s[0] - 'A';
        int idx2 = s[1] - 'A';
        int weight = std::stoi(s.substr(2));

        AdjMat[idx1][idx2] = weight;
        AdjMat[idx2][idx1] = weight;
    }

    return AdjMat;
};



/**
 * @brief Convert ASCII Text to a Directed Adjacency Matrix.
 * 
 * Similar to convert2AdjMat, but for Directed graph.
 * For example:
 * 
 * "AB1,AC1,CD1" is translated to 
 * 0 1 1 0 
 * 0 0 0 0 
 * 0 0 0 1 
 * 0 0 0 0 
 * 
 * @tparam T 
 * @param str 
 * @return Adjmat<T> 
 * @see convert2Adjmat()
 */
template<class T>
Adjmat<T> convert2DirectedAdjmat(std::string str) {
    // trim white space
    str.erase( std::remove(str.begin(), str.end(), ' '), str.end() );

    // seperate into chunck using delimiter ','.
    std::vector<std::string> vs;
    for(auto s_itr=str.begin(), d_itr=s_itr; 
            d_itr!=str.end(); 
            s_itr=d_itr+1) {
        d_itr = std::find(d_itr+1, str.end(), ',');
        vs.emplace_back(std::string(s_itr, d_itr));
    };

    // Find number of nodes
    int min=255, max=0; 
    for (const std::string& s : vs) {
        // Assume node labels are single character long.
        min = min < int(s[0]) ? min : int(s[0]);
        min = min < int(s[1]) ? min : int(s[1]);
        max = int(s[0]) < max ? max : int(s[0]);
        max = int(s[1]) < max ? max : int(s[1]);
    }
    int n_nodes = max - min + 1;

    // Build Adjacency matrix
    Adjmat<T> AdjMat(n_nodes);
    for(const std::string& s : vs) {
        // ASCII table offset between 'A' to 0 is 65.

        int idx1 = s[0] - 'A';
        int idx2 = s[1] - 'A';
        int weight = std::stoi(s.substr(2));

        // Directed
        AdjMat[idx1][idx2] = weight;
    }

    return AdjMat;
}


template<class S>
void PrintAdjmat(Adjmat<S>& mat) {
    const auto NNodes = mat.len();

    for (int i=0; i<NNodes; i++) {
        for (int j=0; j<NNodes; j++) 
            std::cout << mat[i][j] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


}   // namespace


#endif