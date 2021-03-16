//
// Created by mikolajtwarog on 2021-02-08.
//

#ifndef TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <vector>
#include <climits>
#include <queue>

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_maximal_planar.hpp>

#include "problems2.hpp"
#include "../utils/level_face_traversal.h"
#include "../utils/tree_constructor.hpp"

using namespace boost;

namespace boost {
    enum edge_faces_t { edge_faces };

    BOOST_INSTALL_PROPERTY(edge, faces);
}

typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

typedef adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int, EdgeProperty>
        >
        Graph;

typedef std::vector<std::vector< graph_traits<Graph>::edge_descriptor > > PlanarEmbedding;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::vertex_descriptor Vertex;



template<typename Problem>
class baker_impl {

//    void root_tree_with_root(::tree<Problem> &tree, int node, int root);

    Graph gr;
    PlanarEmbedding embedding;
    ::tree<Problem> tree;
    std::vector<int> vertex_level;
    std::set< std::pair<int, int> > added_edges;

    void find_dividing_points(Edge one, Edge two, std::set<int>& dividing_points) {
        int v;
        if (one.m_source == two.m_source || one.m_source == two.m_target) {
            v = one.m_source == two.m_source ? two.m_source : two.m_target;
        } else {
            v = one.m_target == two.m_source ? two.m_source : two.m_target;
        }

        int level = vertex_level[v];

        int one_it = get_edge_it(one, v, embedding);
        int two_it = get_edge_it(two, v, embedding);

        auto& edges = embedding[v];

        for (int i = (one_it + 1) % edges.size(); i != two_it; i = (i + 1) % edges.size()) {
            Edge e = edges[i];
            if (vertex_level[e.m_source] == level - 1 || vertex_level[e.m_target] == level -1) {
                dividing_points.insert(vertex_level[e.m_source] == level - 1 ? e.m_source : e.m_target);
            }
        }

        for (int i = (two_it + 1) % edges.size(); i != one_it; i = (i + 1) % edges.size()) {
            Edge e = edges[i];
            if (vertex_level[e.m_source] == level - 1 || vertex_level[e.m_target] == level -1) {
                dividing_points.insert(vertex_level[e.m_source] == level - 1 ? e.m_source : e.m_target);
            }
        }
    }

    void create_boudaries_rec(::tree<Problem>& t, int node) {
        if (t[node].children.empty()) {
            return;
        }

        for (int i : t[node].children) {
            create_boudaries_rec(t, i);
        }

        t[node].LB = t[t[node].children[0]].LB;
        t[node].RB = t[t[node].children.back()].RB;
    }

    void create_boudaries(::tree<Problem>& t, ::tree<Problem>& t2, int root) {
        if (t.enclosing_tree != nullptr) {
            Problem& f = t.enclosing_tree->t[t.enclosing_face];
            std::vector<int> y_table;
            y_table.push_back(t2[f.children[0]].label.first);

            for (int v : f.children) {
                y_table.push_back(t2[v].label.second);
            }

            std::vector<int> leaves;
            get_leaves(t, leaves, root);
            t[leaves[0]].LB = 0;
            t[leaves.back()].RB = f.children.size();

            for (int j = 1; j < leaves.size(); j++) {
                Problem& v = t[leaves[j]];
                Problem& w = t[leaves[j-1]];
                Edge one(v.label.first, v.label.second, nullptr);
                Edge two(w.label.first, w.label.second, nullptr);

                std::set<int> dividing_points;
                find_dividing_points(one, two, dividing_points);

                for (int i = w.LB; i < y_table.size(); i++) {
                    if (dividing_points.find(y_table[i]) != dividing_points.end()) {
                        v.LB = w.RB = i;
                        break;
                    }
                }
            }

            create_boudaries_rec(t, root);
        } else {
            for (auto &node : t.t) {
                node.LB = node.label.first;
                node.RB = node.label.second;
            }
        }

        for (int i = 0; i < t.size(); i++) {
            Problem& node = t[i];
            if (!node.component_tree.empty()) {
                create_boudaries(node.component_tree, t, node.component_tree.root);
            }
        }
    }

    void table (::tree<Problem>& t, int v) {
        int level = t[v].level;

        if (!t[v].children.empty() && t[v].component_tree.empty()) {
            table(t, t[v].children[0]);
            t[v].val = t[t[v].children[0]].val;

            for(int i=1; i < t[v].children.size(); i++) {
                table(t, t[v].children[i]);
                t[v].merge(t[v].val, t[t[v].children[i]].val);
            }

            t[v].adjust();
        } else if (!t[v].children.empty() && !t[v].component_tree.empty()) {
            auto& ct = t[v].component_tree;
            table(ct, ct.root);
            t[v].contract(ct[ct.root]);
            t[v].adjust();
        } else if (level > 1) {
            std::vector<int> z_table;
            z_table.push_back(t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[0]].label.first);
            for (int x : t.enclosing_tree->t[t.enclosing_face].children) {
                z_table.push_back(t.enclosing_tree->t[x].label.second);
            }

            int p = t[v].RB;
            for (int i = t[v].LB; i < t[v].RB; i++) {
                if (check_for_edge(t[v].label.second, z_table[i], gr, added_edges)) {
                    p = i;
                    break;
                }
            }

            t[v].create(p, gr, added_edges);

            int j = p - 1;

            while (j >= t[v].LB) {
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.first, gr, added_edges);
                t[v].merge(t[v].val, second.val);
                j--;
            }

            j = p;

            while (j < t[v].RB) {
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.second, gr, added_edges);
                t[v].merge(second.val, t[v].val);
                j++;
            }
        }
    }

public:
    baker_impl(const Graph& arg_g): gr(arg_g), embedding(num_vertices(arg_g)), vertex_level(num_vertices(arg_g)),
                                    visited(num_vertices(arg_g)) {
        property_map<Graph, edge_index_t>::type e_index = get(edge_index, gr);
        graph_traits<Graph>::edges_size_type edge_count = 0;
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(gr); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);


        if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = gr,
                                         boyer_myrvold_params::embedding =
                                                 &embedding[0]
        )
                )
            std::cout << "Input graph is planar" << std::endl;
        else
            std::cout << "Input graph is not planar" << std::endl;

        auto component = get(edge_index, gr);
        int bi_num = biconnected_components(gr, component);
        std::vector<Graph> bicomps(bi_num);
        std::vector<PlanarEmbedding> bicomps_emb(bi_num, PlanarEmbedding(num_vertices(gr)));
        std::vector<::tree<Problem> > bicomp_tree(bi_num);
        std::vector< std::vector<int> > bicomps_levels(bi_num, std::vector<int>(num_vertices(gr)));

        for(boost::tie(ei, ei_end) = edges(gr); ei != ei_end; ++ei) {
            add_edge(ei->m_source, ei->m_target, bicomps[component[*ei]]);
        }

        for (int i = 0; i < embedding.size(); i++) {
            for (Edge e : embedding[i]) {
                bicomps_emb[component[e]][i].push_back(e);
            }
        }

//        property_map<Graph, vertex_index_t> index;

        for (int i = 0; i < bi_num; i++) {

            bicomp_tree[i] = build_tree(v, bicomps[i], bicomps_emb[i], bicomps_levels[i]);
        }

//        int k = name_levels();

//        int v = 0;

//        for (int i = 0; i < vertex_level.size(); i++) {
//            if (vertex_level[i] == 1) {
//                v = i;
//                break;
//            }
//        }
//
//        tree = build_tree(v);

//        root_tree_with_root(tree, 0);

//        create_boudaries(tree, tree, 1);

        table(tree, 1);
    }

    int result() {
        return tree[1].result();
    }

};

template <typename Problem>
int baker2(const Graph& g) {
    baker_impl < Problem > b(g);
    return b.result();
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
