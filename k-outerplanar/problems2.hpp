//
// Created by mikolaj on 1/18/21.
//

#ifndef TECHNIKA_BAKER_PROBLEMS_HPP
#define TECHNIKA_BAKER_PROBLEMS_HPP

template <typename Problem>
struct tree{
    Problem* enclosing_face;
    std::vector<Problem> t;
    int root;
    int level;

    tree():enclosing_face(nullptr){}
    explicit tree(int size): t(size), enclosing_face(nullptr){}
    explicit tree(int size, int l): level(l), enclosing_face(nullptr){
        for (int i = 0; i < size; i++) {
            t.emplace_back(level);
        }
    }
    tree(const tree<Problem>& t2) {
        t = t2.t;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        level = t2.level;
    }

    Problem& operator[] (int x) {
        return t[x];
    }

    tree<Problem>& operator= (const tree<Problem>& t2) {
        t = t2.t;
        enclosing_face = t2.enclosing_face;
        root = t2.root;
        root = t2.root;
        level = t2.level;
        return *this;
    }

    int size() {
        return t.size();
    }

    bool empty() {
        return t.empty();
    }

    void emplace_back() {
        t.emplace_back(level);
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
    tree<independent_set>* enclosing_tree;
    int enclosing_face;
    tree<independent_set> component_tree;
    std::vector<int> val;
    int level;

    independent_set(){}

    independent_set(int l): level(l), val(1 << (l + 2)) {
//        val.push_back(0);
//        val.push_back(1);
//        val.push_back(1);
//        val.push_back(-INT_MAX);
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
        return *this;
    }

    // val = {v_0, ..., v_(1 << (level*2)}
    // v_0 <- x
    // v_(1 << level) <- y

    void get_left_boundary(std::vector<int>& lb) {
        lb.push_back(label.first);

        if (enclosing_tree == nullptr)
            return;

        enclosing_tree->t[enclosing_tree->t[enclosing_face].children[LB]].get_left_boundary(lb);
    }

    void get_right_boundary(std::vector<int>& rb) {
        rb.push_back(label.second);

        if (enclosing_tree == nullptr)
            return;

        enclosing_tree->t[enclosing_tree->t[enclosing_face].children[RB - 1]].get_left_boundary(rb);
    }

    void merge(independent_set &two) {
        std::vector<int> copy(val);

        int count = 1 << level;

        for (int u = 0; u < count; u++){
            for (int v = 0; v < count; v++) {
                for (int z = 0; z < count; z++) {
                    int ones = 0;
                    for (int i = 0; i < level; i++) {
                        ones += (z & (1 << i)) > 0;
                    }

                    val[u + (v << level)] = std::max(val[u + (v << level)],
                                                     copy[u + (z << level)] + two.val[z + (u << level)] - ones);
                }
            }
        }
    }

    void adjust() {
        int count = 1 << (level - 1);

        if(label.first == label.second) {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[((u << 1) + 1) + (((v << 1) + 1) << level)]--;
                    val[((u << 1)) + (((v << 1) + 1) << level)] = -INT_MAX;
                    val[((u << 1) + 1) + (((v << 1)) << level)] = -INT_MAX;
                }
            }
        } else {
            for (int u = 0; u < count; u++) {
                for (int v = 0; v < count; v++) {
                    val[((u << 1) + 1) + (((v << 1) + 1) << level)] = -INT_MAX;
                }
            }
        }
    }

    void contract(independent_set& two) {
        int count = 1 << level;

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                val[u + (v << level)] =
                        std::max(two.val[((u << 1) + 1) + (((v << 1) + 1) << level)],
                                 two.val[(u << 1) + ((v << 1) << level)]);
            }
        }
    }

    independent_set extend(int z) {
        int count = 1 << (level - 1);

        independent_set res(level + 1);

        std::vector<int> lb;
        std::vector<int> rb;

        get_left_boundary(lb);
        get_right_boundary(rb);

        int l = *std::find(lb.begin(), lb.end(), z);
        int r = *std::find(rb.begin(), rb.end(), z);

        for (int u = 0; u < count; u++) {
            for (int v = 0; v < count; v++) {
                res.val[(u + (v << level)) << 1] = val[u + (v << level)];

                if (((1 << l) & u) > 0 || ((1 << r) & v) > 0) {
                    res.val[((u + (v << level)) << 1) + 1] = -INT_MAX;
                } else {
                    res.val[((u + (v << level)) << 1) + 1] = val[u + (v << level)] + 1;
                }
            }
        }

        return res;
    }

    template<typename PlanarEmbedding>
    void create(PlanarEmbedding& embedding, int child_num) {
        const std::vector<int>& children = enclosing_tree->t[enclosing_face].children;
        std::vector<int> vertices;

        vertices.push_back(label.first);
        vertices.push_back(label.second);

        if (child_num < children.size()) {
            independent_set& child = enclosing_tree->t[children[child_num]];
            child.get_left_boundary(vertices);
        } else {
            independent_set& child = enclosing_tree->t[children[child_num - 1]];
            child.get_right_boundary(vertices);
        }

        int count = 1 << vertices.size();

        for (int i = 0; i < count; i++) {
            int ones = 0;
            for (int j = 0; j < level; j++) {
                ones += (i & (1 << j)) > 0;
            }

            val[i] = ones;

            for (int j = 0; j < vertices.size(); j++) {
                int v = vertices[j];

                if (((1 << j) & i) == 0) {
                    continue;
                }

                for (auto e : embedding[v]) {
                    int neighbour = e.m_source == v ? e.m_target : e.m_source;
                    int n_num = *std::find(vertices.begin(), vertices.end(), neighbour);

                    if (((1 << n_num) & i) > 0) {
                        val[i] = -INT_MAX;
                        break;
                    }
                }
            }
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
