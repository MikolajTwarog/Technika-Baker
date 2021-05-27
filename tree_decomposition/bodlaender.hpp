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

enum Problem {vc, is, ds_ecc, ds_lcc};

class bodlaender {

    struct tr_node{
        dynamic_bitset<> f;
        std::vector<int> r_values;
        int x;

        tr_node(int n): f(n, 0) {}
        tr_node(int n, int f_num): f(n, f_num) {}
    };

    tree_decomposition tr;
    Graph graph;
    PlanarEmbedding embedding;
    std::vector< std::vector<tr_node> > tables;
    int n;
    Problem problem;
    bool minimum;


    std::vector< std::vector<int> > constraints;

    void calculate_table(int v, int p) {
        std::vector<int>& children = tr.tree[v];
        std::vector<int> vertices(tr.nodes[v].begin(), tr.nodes[v].end());

        for (int child : children) {
            if (child != p) {
                calculate_table(child, v);
            }
        }

        std::vector<tr_node> temp;

        int num = 1 << vertices.size();
        for (int x = 0; x < num; x++) {
            temp.emplace_back(n);
            auto& f = temp.back().f;
            for (int i = 0; i < vertices.size(); i++) {
                if ((x >> i) & 1) {
                    f[vertices[i]] = 1;
                }
            }

            if (problem == ds_ecc) {
                for (auto &con : constraints) {
                    int r = 0;
                    for (int ver : con) {
                        r += f[ver];
                    }
                    temp.back().r_values.push_back(r);
                }
            } else {
                int r_val = 1;
                for (auto &con : constraints) {
                    int r = 0;
                    for (int ver : con) {
                        if (tr.nodes[v].find(ver) != tr.nodes[v].end()) {
                            r += f[ver];
                        } else {
                            r = 1;
                            break;
                        }
                    }
                    if (minimum) {
                        r_val &= (r > 0);
                    } else {
                        r_val &= (r < 2);
                    }
                }
                temp.back().r_values.push_back(r_val);
            }

            temp.back().r_values.push_back(f.count());
            temp.back().x = x;
        }

        for (int child : children) {
            if (child != p) {
                temp = calculate_temp(v, child, temp);
            }
        }

        tables[v] = temp;
    }

    std::vector<tr_node> calculate_temp (int v, int child, std::vector<tr_node>& temp) {
        std::vector<tr_node> new_temp;
        std::vector<int> intersection;
        std::set_intersection(tr.nodes[v].begin(), tr.nodes[v].end(), tr.nodes[child].begin(),
                              tr.nodes[child].end(), std::back_inserter(intersection));

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

                auto f3 = f1.f & f2.f;
                std::vector<int> s_values;
                if (problem == ds_ecc) {
                    for (int i = 0; i < constraints.size(); i++) {
                        auto &con = constraints[i];
                        int r = f1.r_values[i] + f2.r_values[i];
                        for (int ver : con) {
                            r -= f3[ver];
                        }
                        s_values.push_back(r);
                    }
                } else {
                    s_values.push_back(f1.r_values[0] & f2.r_values[0]);
                }

                s_values.push_back(f1.r_values.back() + f2.r_values.back() - f3.count());

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
                    new_temp.emplace_back(n);
                    new_temp.back().f = f1.f;
                    new_temp.back().r_values = s_values;
                    new_temp.back().x = f1.x | f2.x;
                } else {
                    if ((minimum && new_temp[f_it].r_values.back() > s_values.back())
                    || (!minimum && new_temp[f_it].r_values.back() < s_values.back())) {
                        new_temp[f_it].r_values.back() = s_values.back();
                        new_temp[f_it].f = f1.f;
                        new_temp[f_it].x = f1.x | f2.x;
                    }
                }
            }
        }

        return new_temp;
    }

public:
    bodlaender(Graph g, PlanarEmbedding emb, std::vector<int> outer_face, Problem pr): graph(g), embedding(emb),
    n(num_vertices(g)), problem(pr) {
        minimum = pr == vc || pr == ds_ecc || pr == ds_lcc;
        get_tree_decomposition(g, emb, outer_face, tr);
        tables.resize(tr.tree.size());

        if (problem == vc || problem == is) {
            graph_traits<Graph>::edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
                constraints.emplace_back();
                constraints.back().push_back(ei->m_source);
                constraints.back().push_back(ei->m_target);
            }
        }
        if (problem == ds_ecc || problem == ds_lcc) {
            for (int v = 0; v < embedding.size(); v++) {
                constraints.emplace_back();
                constraints.back().push_back(v);

                for (auto& e : embedding[v]) {
                    int w = v == e.m_source ? e.m_target : e.m_source;
                    constraints.back().push_back(w);
                }
            }

            if (problem == ds_lcc) {
                for (auto &node : tr.nodes) {
                    std::vector<int> to_add;
                    for (int v : node) {
                        for (auto &e : embedding[v]) {
                            int w = v == e.m_source ? e.m_target : e.m_source;
                            to_add.push_back(w);
                        }
                    }

                    for (int v : to_add) {
                        node.insert(v);
                    }
                }
            }
        }

        calculate_table(0, -1);
    }

    int get_result() {
        int result;
        if (minimum) {
            result = INT16_MAX;
        } else {
            result = -1;
        }

        for (auto& node : tables[0]) {
            int good = true;
            for (int r = 0; r < node.r_values.size() - 1; r++) {
                if (node.r_values[r] == 0) {
                    good = false;
                    break;
                }
            }

            if (good) {
                if (minimum) {
                    result = std::min(result, node.r_values.back());
                } else {
                    result = std::max(result, node.r_values.back());
                }
            }
        }

        return result;
    }
};

int bodlaender_vertex_cover(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender bd(g, emb, outer_face, vc);
    return bd.get_result();
}

int bodlaender_independent_set(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender bd(g, emb, outer_face, is);
    return bd.get_result();
}

int bodlaender_dominating_set_lcc(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender bd(g, emb, outer_face, ds_lcc);
    return bd.get_result();
}

int bodlaender_dominating_set_ecc(Graph& g, PlanarEmbedding& emb, std::vector<int>& outer_face) {
    bodlaender bd(g, emb, outer_face, ds_ecc);
    return bd.get_result();
}


#endif //TECHNIKA_BAKER_BODLAENDER_HPP
