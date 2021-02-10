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

#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "problems2.hpp"
#include "../utils/level_face_traversal.h"

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

template <typename Edge>
struct my_visitor : public planar_face_traversal_visitor
{
    std::map<Edge, std::vector<int> > &faces;
    std::vector<std::vector<int> >& vertices_in_face;
    int current_face = 0;

    my_visitor(std::map<Edge, std::vector<int> >& f, std::vector<std::vector<int> >& o)
            : faces(f), vertices_in_face(o){ }

    void begin_face() {
        vertices_in_face.emplace_back();
    }

    void end_face() {
//        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
//        std::cout << e << " " << current_face << "\n";
        faces[e].push_back(current_face);
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        vertices_in_face[current_face].push_back(v);
    }
};


template <typename Edge, typename Problem, typename PlanarEmbedding>
struct tree_builder : public planar_face_traversal_visitor
{
    std::map<Edge, std::vector<int> > &faces;
    std::vector<Problem> &tree;
    int current_face = 0;
    int last_vertex;
    Graph graph;
    PlanarEmbedding embedding;

    tree_builder(std::map<Edge, std::vector<int> >& f, std::vector<Problem>& t, Graph& g, PlanarEmbedding& emb)
            : faces(f), tree(t), graph(g), embedding(emb){ }

    void end_face() {
        current_face++;
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        tree[current_face].face.push_back(v);
        last_vertex = v;
    }

    void next_edge(Edge e)
    {
        if(current_face > 0) {
            int neighbor = faces[e][0] == current_face ? faces[e][1] : faces[e][0];

//            std::cout << "\n" << e << "\n";
//            for (int i : faces[e]) {
//                std::cout << i << "\n";
//            }

            if (neighbor == 0) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
                tree[last].label.first = last_vertex;
                tree[last].label.second = last_vertex == e.m_source ? e.m_target : e.m_source;

//                if (e != *embedding[target(e, graph)].begin())
//                    std::swap(tree[last].label.first, tree[last].label.second);
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};

template<typename Problem>
class baker_impl {

//    void root_tree_with_root(std::vector<Problem> &tree, int node, int root);

    Graph g;
    PlanarEmbedding embedding;
    std::vector<Problem> tree;
    std::vector<int> vertex_level;

    int name_levels() {
        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        planar_face_traversal(g, &embedding.front(), my_vis);

        for (int &v : vertex_level) {
            v = -1;
        }

        Edge next_level_edge;
        Edge current_edge;
        int starting_v;
        int current_v;
        int current_egde_it;

        typename boost::graph_traits<Graph>::edge_iterator fi, fi_end;
        for (boost::tie(fi, fi_end) = edges(g); fi != fi_end; ++fi) {
            if (faces[*fi][0] == 0 || faces[*fi][1] == 0) {
                next_level_edge = *fi;
                current_edge = next_level_edge;
                starting_v = next_level_edge.m_source;
                current_v = next_level_edge.m_target;
                vertex_level[current_v] = 0;
                break;
            }
        }

        current_egde_it = get_edge_it(current_edge, current_v);

        std::vector<int> current_level;
        current_level.push_back(starting_v);
        current_level.push_back(current_v);

        for (int level = 0; ; level++) {
            while (current_v != starting_v) {
                Edge e = embedding[current_v][current_egde_it];

                for (int j = (current_egde_it + 1) % embedding[current_v].size(); j != current_egde_it;
                     j = (j + 1) % embedding[current_v].size()) {
                    Edge e_j = embedding[current_v][j];
                    if (vertex_level[source(e_j, g)] == -1
                        || vertex_level[target(e_j, g)] == -1) {
                        current_edge = e_j;
                        break;
                    }
                }

                current_v = current_edge.m_source == current_v ? current_edge.m_target : current_edge.m_source;
                current_egde_it = get_edge_it(current_edge, current_v);

                current_level.push_back(current_v);
                vertex_level[current_v] = level;
            }

            vertex_level[starting_v] = level;

            for (int v : current_level) {
                for (Edge e : embedding[v]) {
                    if (vertex_level[e.m_source] == -1
                        || vertex_level[e.m_target] == -1) {
                        next_level_edge = e;
                        starting_v = e.m_source == v ? e.m_target : e.m_source;
                    }
                }
            }

            if (vertex_level[next_level_edge.m_source] > -1
                && vertex_level[next_level_edge.m_target] > -1) {
                return level + 1;
            }

            current_egde_it =
                    (get_edge_it(next_level_edge, starting_v) + 1) % embedding[starting_v].size();
            for (int j = (current_egde_it + 1) % embedding[starting_v].size(); j != current_egde_it;
                 j = (j + 1) % embedding[starting_v].size()) {
                Edge e_j = embedding[starting_v][j];
                if (vertex_level[source(e_j, g)] == -1
                    || vertex_level[target(e_j, g)] == -1) {
                    current_edge = e_j;
                    break;
                }
            }
            current_v = current_edge.m_source == starting_v ? current_edge.m_target : current_edge.m_source;

            current_level.clear();
            current_level.push_back(starting_v);
            current_level.push_back(current_v);

            vertex_level[current_v] = level + 1;
        }


    }

    void tringulate() {
        //TODO
    }

    int find_third(int one, int two) {
        auto& one_edges = embedding[one];
        auto& two_edges = embedding[two];

        int level = vertex_level[one];

        if (one == two) {
            for (auto e_i : one_edges) {
                int target_i = e_i.m_source == one ? e_i.m_target : e_i.m_source;

                if (vertex_level[target_i] == level + 1) {
                    return target_i;
                }
            }
        }

        for (auto e_i : one_edges) {
            int target_i = e_i.m_source == one ? e_i.m_target : e_i.m_source;

            if (vertex_level[target_i] == level -1) {
                auto& target_i_edges = embedding[target_i];

                for (auto e_j : target_i_edges) {
                    int target_j = e_j.m_source == target_i ? e_j.m_target : e_j.m_source;
                    if (target_j == two) {
                        return target_i;
                    }
                }
            }
        }

        return -1;
    }

    void root_tree_2(std::vector<Problem> &t, int node) {
        if (t[node].children.empty())
            return;

        int parent_it = 0;
        auto child = t[node].children.begin();
        for (int i = 0; child != t[node].children.end(); child++, i++) {
            if (*child == t[node].parent) {
                t[node].children.erase(child--);
                parent_it = i;
            } else {
                t[*child].parent = node;
                root_tree_2(t, *child);
            }
        }

        int last = t[node].children.size() - 1;

        std::rotate(t[node].children.begin(),
                    t[node].children.begin() + parent_it,
                    t[node].children.end());

        t[node].label.first = t[t[node].children[0]].label.first;
        t[node].label.second = t[t[node].children[last]].label.second;

        if (!t[node].component_tree.empty()) {
            tringulate();
            int v = find_third(t[node].label.first, t[node].label.second);
            root_tree_with_root(t[node].component_tree, v);
        }
    }

    void root_tree_with_root(std::vector<Problem> &t, int root) {
        int node;

        for (int i = 0; i < t.size(); i++) {
            for (int v : t[i].face) {
                if (v == root) {
                    node = i;
                    break;
                }
            }
        }

        root_tree_2(t, node);

        auto &children = t[node].children;

        for (int i = 0; i < children.size(); i++) {
            int child = children[i];

            if (t[child].label.first == root) {
                std::rotate(t[node].children.begin(),
                            t[node].children.begin() + i,
                            t[node].children.end());
                break;
            }
        }

        t[node].label.first = root;
        t[node].label.second = root;
    }

    template<typename Edge>
    int get_edge_it(Edge e, int v) {
        for (int i = 0; i < embedding[v].size(); i++) {
            if (embedding[v][i].m_source == e.m_source
                && embedding[v][i].m_target == e.m_target) {
                return i;
            }
        }

        return -1;
    }


    int check_for_component(const std::vector<int> &face) {
        int level = vertex_level[face[0]];

        for (int v : face) {
            for (auto e : embedding[v]) {
                int target = e.m_source == v ? e.m_target : e.m_source;
                if (vertex_level[target] == level + 1) {
                    return target;
                }
            }
        }

        return -1;
    }

    std::vector<Problem> build_tree(int v) {
        int level = vertex_level[v];

        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        level_face_traversal(g, embedding, my_vis, level, vertex_level);
        std::vector<Problem> t(my_vis.current_face);

        tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, t, g, embedding);
        level_face_traversal(g, embedding, tree_b, level, vertex_level);

        for (int i = 0; i < vertices_in_face.size(); i++) {
            auto &face = vertices_in_face[i];

            int v_in_c = check_for_component(face);
            if (v_in_c > -1) {
                t[i].component_tree = build_tree(v_in_c);
            }
        }

        return t;
    }

public:
    baker_impl(const Graph& arg_g): g(arg_g), embedding(num_vertices(arg_g)), vertex_level(num_vertices(arg_g)) {
        property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
        graph_traits<Graph>::edges_size_type edge_count = 0;
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);


        if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                         boyer_myrvold_params::embedding =
                                                 &embedding[0]
        )
                )
            std::cout << "Input graph is planar" << std::endl;
        else
            std::cout << "Input graph is not planar" << std::endl;

        int k = name_levels();

        int v = 0;

        for (int i = 0; i < vertex_level.size(); i++) {
            if (vertex_level[i] == 0) {
                v = i;
                break;
            }
        }

        tree = build_tree(v);

        root_tree_with_root(tree, 0);
    }

    int result() {
        return 0;
    }

};

template <typename Problem>
int baker2(const Graph& g) {
    baker_impl < Problem > b(g);
    return b.result();
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
