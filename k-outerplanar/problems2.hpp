//
// Created by mikolaj on 1/18/21.
//

#ifndef TECHNIKA_BAKER_PROBLEMS_HPP
#define TECHNIKA_BAKER_PROBLEMS_HPP

#include "../utils/level_face_traversal.h"

template <typename Problem>
struct tree{
    tree<Problem>* enclosing_tree;
    int enclosing_face;
    std::vector<Problem> t;
    int root;
    int outer_face;
    int level;

    tree():enclosing_tree(nullptr){}
    explicit tree(int size): t(size), enclosing_tree(nullptr){}
    explicit tree(int size, int l): level(l), enclosing_tree(nullptr){
        for (int i = 0; i < size; i++) {
            t.emplace_back(l, this);
//            t.back().my_tree = this;
        }
    }
    tree(const tree<Problem>& t2) {
        t = t2.t;
        enclosing_tree = t2.enclosing_tree;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        level = t2.level;
        for (auto& p : t) {
            p.my_tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
    }

    Problem& operator[] (int x) {
        return t[x];
    }

    tree<Problem>& operator= (const tree<Problem>& t2) {
        t = t2.t;
        enclosing_tree = t2.enclosing_tree;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        level = t2.level;
        for (auto& p : t) {
            p.my_tree = this;
            if (!p.component_tree.empty()) {
                p.component_tree.enclosing_tree = this;
            }
        }
        return *this;
    }

    int size() {
        return t.size();
    }

    bool empty() {
        return t.empty();
    }

    void emplace_back() {
        t.emplace_back(level, this);
//        t.back().my_tree = this;
    }
};

struct node
{
    int parent;
    std::pair<int, int> label;
    int LB;
    int RB;
    std::vector<int> children;
    std::vector<int> face;
    node(): parent(-1) {}
};

struct independent_set : node
{
    tree<independent_set>* my_tree;
    tree<independent_set> component_tree;
    std::vector< std::vector<int> > val;
    int level;

    independent_set(int l): independent_set(l, nullptr) {}

    independent_set(int l, tree<independent_set>* mt):
    level(l), val(1 << l,std::vector<int>(1 << l)), my_tree(mt) {
        if(l == 1) {
            val[0][0] = 0;
            val[0][1] = 1;
            val[1][0] = 1;
            val[1][1] = -INT16_MAX;
        }
    }

    independent_set(const independent_set& two): my_tree(two.my_tree) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        level = two.level;
//        my_tree = two.my_tree;
    }

    independent_set& operator=(const independent_set& two) {
        parent = two.parent;
        label = two.label;
        LB = two.LB;
        RB = two.RB;
        children = two.children;
        face = two.face;
        val = two.val;
        component_tree = two.component_tree;
        my_tree = two.my_tree;
        level = two.level;
        return *this;
    }

    // val = {v_0, ..., v_(1 << (level*2)}
    // v_0 <- x
    // v_(1 << level) <- y

    void get_left_boundary(std::vector<int>& lb) {
        lb.push_back(label.first);

        if (my_tree->enclosing_tree == nullptr)
            return;

        my_tree->enclosing_tree->t[my_tree->enclosing_tree->t[my_tree->enclosing_face].children[LB]].get_left_boundary(lb);
    }

    void get_right_boundary(std::vector<int>& rb) {
        rb.push_back(label.second);

        if (my_tree->enclosing_tree == nullptr)
            return;

        my_tree->enclosing_tree->t[my_tree->enclosing_tree->t[my_tree->enclosing_face].children[RB - 1]].get_left_boundary(rb);
    }

    void merge(std::vector< std::vector<int> > one, std::vector< std::vector<int> > two) {
//        std::vector<int> copy(val);

        int count = 1 << level;

        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                val[u][v] = -INT16_MAX;
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }

                    val[u][v] = std::max(val[u][v], one[u][z] + two[z][v] - ones);
                }
            }
        }
    }

    void adjust() {
        int count = 1 << (level - 1);

        if(label.first == label.second) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u << 1) + 1][(v << 1) + 1]--;
                    val[u << 1][(v << 1) + 1] = -INT16_MAX;
                    val[(u << 1) + 1][v << 1] = -INT16_MAX;
                }
            }
        } else {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[(u << 1) + 1][(v << 1) + 1] = -INT16_MAX;
                }
            }
        }
    }

    void contract(independent_set& two) {
        int count = 1 << level;

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u][v] = std::max(two.val[(u << 1) + 1][(v << 1) + 1], two.val[u << 1][v << 1]);
            }
        }
    }

    template<typename Graph>
    independent_set extend(int z, Graph& g, std::set<std::pair<int, int> >& ae) {
        int count = 1 << level;

        independent_set res(level + 1);

        std::vector<int> lb;
        std::vector<int> rb;

        get_left_boundary(lb);
        get_right_boundary(rb);

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[u << 1][v << 1] = val[u][v];
                res.val[u << 1][(v << 1) + 1]= -INT16_MAX;
                res.val[(u << 1) + 1][v << 1] = -INT16_MAX;

                if (((u & 1) == 1 && check_for_edge(lb[0], z, g, ae))
                    || ((v & 1) == 1 && check_for_edge(rb[0], z, g, ae))){
                    res.val[(u << 1) + 1][(v << 1) + 1] = -INT16_MAX;
                } else {
                    res.val[(u << 1) + 1][(v << 1) + 1] = val[u][v] + 1;
                }
            }
        }

        return res;
    }

    template<typename Graph>
    void create(int child_num, Graph& g, std::set<std::pair<int, int> >& ae) {
        const std::vector<int>& children = my_tree->enclosing_tree->t[my_tree->enclosing_face].children;
        std::vector<int> vertices;

        if (child_num < children.size()) {
            independent_set& child = my_tree->enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            independent_set& child = my_tree->enclosing_tree->t[children[child_num - 1]];
            child.get_right_boundary(vertices);
        }

        int count = 1 << vertices.size();

        for (int u = 0; u < (count << 1); u++) {
            for (int v = 0; v < (count << 1); v++) {
                val[u][v] = -INT16_MAX;
            }
        }

        for (int i = 0; i < count; i++) {

            bool bad = false;
            for (int v = 0; v < vertices.size() - 1; v++) {
                if (((i >> v) & 1) && ((i >> (v + 1)) & 1)) {
                    bad = true;
                    break;
                }
            }

            if (bad) {
                continue;
            }

            int ones = 0;
            for (int j = 0; j < vertices.size(); j++) {
                ones += (i & (1 << j)) > 0;
            }

            val[i << 1][i << 1] = ones;

            if ((i & 1) == 0 || !check_for_edge(vertices[0], label.first, g,  ae)) {
                val[(i << 1) + 1][i << 1] = ones + 1;
            }

            if ((i & 1) == 0 || !check_for_edge(vertices[0], label.second, g, ae)) {
                val[i << 1][(i << 1) + 1] = ones + 1;
            }
        }
    }

    int result() {
        return std::max(val[0][0], val[1][1]);
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
