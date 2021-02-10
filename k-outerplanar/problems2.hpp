//
// Created by mikolaj on 1/18/21.
//

#ifndef TECHNIKA_BAKER_PROBLEMS_HPP
#define TECHNIKA_BAKER_PROBLEMS_HPP

struct node
{
    int parent;
    std::pair<int, int> label;
    std::vector<int> children;
    std::vector<int> face;
    node(): parent(-1) {}
};

struct independent_set : node
{
    std::vector<independent_set> component_tree;
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
    std::vector<vertex_cover> component_tree;
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
    std::vector<dominating_set> component_tree;
    std::vector<int> val;
    int child_num = 1;
    dominating_set() {
        val.push_back(INT16_MAX-1); //0
        val.push_back(1);           //1
        val.push_back(INT16_MAX-1); //2
        val.push_back(1);           //3
        val.push_back(2);           //4
        val.push_back(1);           //5
        val.push_back(INT16_MAX-1); //6
        val.push_back(1);           //7
        val.push_back(0);           //8
    }

    void merge(dominating_set &two) {
        std::vector<int> copy(val);
        for(int i=0; i<9; i++){
            int b1 = i / 3;
            int b2 = i % 3;
            val[i] = std::min(copy[b1*3 + 1] + two.val[b2 + 3] - 1, copy[b1*3 + 2] + two.val[b2]);
            val[i] = std::min(val[i], copy[b1*3] + two.val[b2 + 6]);
        }
    }

    void adjust() {
        if(parent == -1) {
            val[4]--;
        } else {
            val[3] = val[5];
            val[1] = val[7];
        }
    }

    int result() {
        return std::min(val[0], val[4]);
    }
};

#endif //TECHNIKA_BAKER_PROBLEMS_HPP
