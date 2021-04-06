#ifndef STATANALY_ADJACENCY_H_
#define STATANALY_ADJACENCY_H_

#include <array>
#include <string>
#include <vector>

namespace statanaly {

/* Adjacency List */

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


/* Adjacency Matrix */
// Diagonal terms are deafult to zero.

template<class SCALAR>
class Adjmat {
    using S = std::size_t;
    S length = 0;

public:
    SCALAR** data = nullptr;
    
    Adjmat() = default;
    Adjmat(S len_) {allocate(len_);}
    ~Adjmat() {deallocate();}

    // Copy operations (aka shallow-copy) are not allowed.
    // because when other goes out of scope, the heap memory will be destroyed.
    Adjmat(const Adjmat& other) = delete;
    Adjmat& operator = (const Adjmat& other) = delete;

    // Move operations reassign the ownership of the data.
    Adjmat(Adjmat&& other) {
        resize(other.len());
        data = std::move(other.data);
        
        other.data = nullptr;
        other.resize();
    };
    Adjmat& operator = (Adjmat&& other) {
        resize(other.len());
        data = std::move(other.data);
        
        other.data = nullptr;
        other.resize();
        return *this;
    };

    // Cloning -- Deep-copy while retaining original data.
    void clone(const Adjmat& other) {
        if (length != other.len()) 
            allocate(other.len());
        
        deepcopy(other.data);
    }

    void allocate(const S len) {
        if (length>0 || data!=nullptr) 
            deallocate();

        // allocate heap space.
        length = len;
        data = new SCALAR* [length];
        for (S i=0; i<length; i++)
            data[i] = new SCALAR [length] {};
    }

    void deallocate() {
        if (length==0 || data==nullptr) 
            return;

        // deallocate heap space.
        for (S i=0; i<length; i++) 
            delete[] data[i];
        delete[] data;
        length = 0;
    }

    void deepcopy(const SCALAR* const * const other) {
        for (S i=0; i<length; i++) 
            memcpy(data[i], other[i], sizeof(SCALAR)*length);
    }

    inline       SCALAR* operator[] (S idx)       { return data[idx]; }
    inline const SCALAR* operator[] (S idx) const { return data[idx]; }

    const S len() const { return length; }
    void resize(const S n=0) { length=n; }
};


/* Convert ASCII Text to an Adjacency matrix for 
   Undirected graph */
// "AB1,AC1,CD1" is mapped to 
// 0 1 1 0 
// 1 0 0 0 
// 1 0 0 1 
// 0 0 1 0 
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
        max = int(s[1]) < max ? max : int(s[1]);
    }
    const int n_nodes = max - min + 1;


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



/* Similar to convert2AdjMat, but for 
   Directed graph */
// "AB1,AC1,CD1" is mapped to 
// 0 1 1 0 
// 0 0 0 0 
// 0 0 0 1 
// 0 0 0 0 
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
        max = int(s[1]) < max ? max : int(s[1]);
    }
    const int n_nodes = max - min + 1;

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