//
// Created by mikolaj on 1/18/21.
//

#ifndef TECHNIKA_BAKER_PROBLEMS_HPP
#define TECHNIKA_BAKER_PROBLEMS_HPP

struct node
{
    int parent;
    std::vector<int> children;
    node(): parent(-1) {}
};

struct independent_set : node
{
    std::vector<int> val;
    independent_set() {
        val.push_back(0);
        val.push_back(1);
        val.push_back(1);
        val.push_back(-INT_MAX);
    }

    void merge(independent_set &two) {
        std::vector<int> copy(val);

        for(int i=0; i<4; i++){
            for(int j=0; j<2; j++) {
                val[i] = std::max(val[i], copy[(i&2) + j] + two.val[((j<<1) + (i&1))] - j);
            }
        }
    }

    void adjust() {
        if(parent == -1) {
            val[1] = -INT_MAX;
            val[2] = -INT_MAX;
            val[3]--;
        } else {
            val[3] = -INT_MAX;
        }
    }

    int result() {
        return std::max(val[0], val[3]);
    }
};

struct vertex_cover : node
{
    std::vector<int> val;
    vertex_cover() {
        val.push_back(INT16_MAX-1);
        val.push_back(1);
        val.push_back(1);
        val.push_back(2);
    }

    void merge(vertex_cover &two) {
        std::vector<int> copy(val);
        for(int i=0; i<4; i++){
            val[i] = INT16_MAX-1;
            for(int j=0; j<2; j++) {
                val[i] = std::min(val[i], copy[(i&2) + j] + two.val[((j<<1) + (i&1))] - j);
            }
        }
    }

    void adjust() {
        if(parent == -1) {
            val[1] = INT16_MAX-1;
            val[2] = INT16_MAX-1;
            val[3]--;
        } else {
            val[0] = INT16_MAX-1;
        }
    }

    int result() {
        return std::min(val[0], val[3]);
    }
};

struct dominating_set : node
{
    std::vector<int> val;
    dominating_set() {
        val.push_back(INT16_MAX-1);
        val.push_back(1);
        val.push_back(1);
        val.push_back(2);
    }

    void merge(dominating_set &two) {
//        if(two.children.empty())
//            two.val[0] = 0;
        std::vector<int> copy(val);
        for(int i=0; i<4; i++){
            val[i] = INT16_MAX-1;
            for(int j=0; j<2; j++) {
                val[i] = std::min(val[i], copy[(i&2) + j] + two.val[((j<<1) + (i&1))] - j);
            }
        }
    }

    void adjust() {
        if(parent == -1) {
            val[1] = INT16_MAX-1;
            val[2] = INT16_MAX-1;
            val[3]--;
        }
    }

    int result() {
        return std::min(val[0], val[3]);
    }
};

#endif //TECHNIKA_BAKER_PROBLEMS_HPP
