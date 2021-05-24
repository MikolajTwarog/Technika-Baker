//
// Created by mikolajtwarog on 2021-05-24.
//

#ifndef TECHNIKA_BAKER_BODLAENDER_HPP
#define TECHNIKA_BAKER_BODLAENDER_HPP

#include <vector>
#include <list>
#include <queue>
#include <set>
#include <map>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include "create_tree_decomposition.hpp"

using namespace boost;

class bodlaender {

    tree_decomposition tr;
    Graph graph;
    PlanarEmbedding embedding;
    std::vector< std::vector<node> > tables;
    int n, m;
    std::vector< std::pair<int, int> > edges;

    struct node{
        dynamic_bitset<> f;
        std::vector<int> r_values;

        node(): f(n) {}
        node(int f_num): f(n, f_num) {}
    };

    void calculate_table(int v, int p) {
        std::vector<int>& children = tr.tree[v];
        std::vector<int> vertices(tr.nodes[v]);

        for (int child : children) {
            if (child != p) {
                calculate_table(child, v);
            }
        }

        std::vector<node> temp;

        int num = 2 << vertices.size();
        for (int x = 0; x < num; x++) {
            temp.emplace_back(x);
            for (auto & e : edges) {
                temp.back().r_values.push_back(x[e.first] + x[e.second]);
            }
            temp.back().push_back(x.count());
        }

        for (int child : children) {
            if (child != p) {
                temp = calculate_temp(v, child, temp);
            }
        }

        tables[v] = temp;
    }

    std::vector<node> calculate_temp (int v, int child, std::vector<node>& temp) {
        std::vector<node> new_temp;
        std::vector<int> intersection;
        std::set_intersection(tr.nodes[v].begin(), tr.nodes[v].end(), tr.nodes[child].begin(), r.nodes[child].end()
        intersection.begin());

        for (auto& f1 : temp) {
            for (auto& f2 : tables[child]) {
                bool good = true;
                for (int ver : intersection) {
                    if (f1.f[ver] != f2.f[ver]) {
                        good = false;
                        break;
                    }
                }

                if (!good) {
                    continue;
                }

                std::vector<int> s_values;
                for (int i = 0; i < f1.r_values.size(); i++) {
                    s_values.push_back(f1.r_values[i] + f2.r_values[i]);
                }

                int f_it = -1;
                for (int i = 0; i < new_temp.size(); i++) {
                    good = true;
                    for (int r = 0; r < s_values.size() - 1; r++) {
                        if (s_values[r] != new_temp[i].r_values[r]) {
                            good = false;
                            break;
                        }
                    }

                    if (good && f1.f == new_temp[i].f) {
                        f_it = i;
                        break;
                    }
                }

                if (f_it == -1) {
                    new_temp.emplace_back();
                    new_temp.back().f = f1.f;
                    new_temp.back().r_values = s_values;
                } else {
                    new_temp[f_it].r_values.back() = std::min(new_temp[f_it].r_values.back(), s_values.back());
                }
            }
        }

        return new_temp;
    }

    bodlaender(Graph g, PlanarEmbedding emb): graph(g), embedding(emb), n(num_vertices(g)), m(num_edges(g)) {
        get_tree_decomposition(g, emb, tr);
        tables.resize(tr.tree.size());
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
            edges.emplace_back(ei->m_source, ei->m_target);
        }

        calculate_table(0, -1);
    }

    int get_result() {
        int result = INT16_MAX;

        for (auto& node : tables[0]) {
            int good = true;
            for (int r = 0; r < node.r_values.size() - 1; r++) {
                if (node.r_values[r] == 0) {
                    good = false;
                    break;
                }
            }

            if (good) {
                result = std::min(result, node.r_values.back());
            }
        }

        return result;
    }
};


#endif //TECHNIKA_BAKER_BODLAENDER_HPP
