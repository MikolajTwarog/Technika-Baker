//
// Created by mikolajtwarog on 2021-04-28.
//

#ifndef TECHNIKA_BAKER_CYCLIC_VECTOR_H
#define TECHNIKA_BAKER_CYCLIC_VECTOR_H


#include <vector>

template <typename T>
class cyclic_vector : public std::vector<T> {
public:
    T& operator[] (int x) {
        x %= this->size();
        if (x < 0) {
            x += this->size();
        }
        return this->at(x);
    }
};


#endif //TECHNIKA_BAKER_CYCLIC_VECTOR_H
