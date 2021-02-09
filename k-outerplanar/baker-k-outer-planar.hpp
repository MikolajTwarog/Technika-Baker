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
    int parent_it = -1;
    int edge_it = 0;
    Graph graph;
    PlanarEmbedding embedding;

    tree_builder(std::map<Edge, std::vector<int> >& f, std::vector<Problem>& t, Graph& g, PlanarEmbedding& emb)
            : faces(f), tree(t), graph(g), embedding(emb){ }

    void end_face() {
        current_face++;
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
                tree[last].label.first = e.m_source;
                tree[last].label.second = e.m_target;

                if (e != *embedding[target(e, graph)].begin())
                    std::swap(tree[last].label.first, tree[last].label.second);
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};



template <typename Problem>
void root_tree_2(std::vector<Problem> &tree, int node) {
    if (tree[node].children.empty())
        return;

    int parent_it = 0;
    auto child = tree[node].children.begin();
    for(int i = 0; child != tree[node].children.end(); child++, i++) {
        if(*child == tree[node].parent) {
            tree[node].children.erase(child--);
            parent_it = i;
        } else {
            tree[*child].parent = node;
            root_tree(tree, *child);
        }
    }

    int last = tree[node].children.size() - 1;

    std::rotate(tree[node].children.begin(),
                tree[node].children.begin()+parent_it,
                tree[node].children.end());

    tree[node].label.first = tree[tree[node].children[0]].label.first;
    tree[node].label.second = tree[tree[node].children[last]].label.second;
}

template <typename PlanarEmbedding, typename Edge>
int get_edge_it(const PlanarEmbedding& embedding, Edge e, int v) {
    for (int i = 0; i < embedding[v].size(); i++) {
        if (embedding[v][i].m_source == e.m_source
            && embedding[v][i].m_target == e.m_target) {
            return i;
        }
    }

    return -1;
}

template <typename PlanarEmbedding>
void name_levels(Graph& g, const PlanarEmbedding& embedding, std::vector<int>& vertex_level, int k) {
    typedef typename graph_traits<Graph>::edge_descriptor Edge;

    std::map<Edge, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    my_visitor<Edge> my_vis(faces, vertices_in_face);
    planar_face_traversal(g, &embedding.front(), my_vis);

    for (int & v : vertex_level) {
        v = -1;
    }

    Edge next_level_edge;
    Edge current_edge;
    int starting_v;
    int current_v;
    int current_egde_it;

    typename boost::graph_traits< Graph >::edge_iterator fi, fi_end;
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

    current_egde_it = get_edge_it(embedding, current_edge, current_v);

    std::vector<int> current_level;
    current_level.push_back(starting_v);
    current_level.push_back(current_v);

    for (int level = 0; level < k; level++) {
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
            current_egde_it = get_edge_it(embedding, current_edge, current_v);

            current_level.push_back(current_v);
            vertex_level[current_v] = level;
        }

        if (level < k - 1) {
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

            current_egde_it = (get_edge_it(embedding, next_level_edge, starting_v) + 1) % embedding[starting_v].size();
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



}

template <typename PlanarEmbedding>
int check_for_component(const PlanarEmbedding& embedding, const std::vector<int>& face,
                         const std::vector<int>& vertex_level) {
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

template <typename Problem, typename PlanarEmbedding>
std::vector<Problem> build_tree(Graph &g, PlanarEmbedding &embedding, int v, const std::vector<int>& vertex_level) {
    int level = vertex_level[v];

    std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    my_visitor<graph_traits<Graph>::edge_descriptor> my_vis(faces, vertices_in_face);
    level_face_traversal(g, embedding, my_vis, level, vertex_level);
    std::vector<Problem> tree(my_vis.current_face);

    tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, tree, g, embedding);
    level_face_traversal(g, embedding, tree_b, level, vertex_level);

    for (int i = 0; i < vertices_in_face.size(); i++) {
        auto& face = vertices_in_face[i];

        int v_in_c = check_for_component(embedding, face, vertex_level);
        if (v_in_c > -1) {
            tree[i].component_tree = build_tree<Problem>(g, embedding, v_in_c, vertex_level);
        }
    }

    return tree;
}

template <typename Problem, typename PlanarEmbedding>
int baker2(Graph &g, PlanarEmbedding &embedding, int k) {
    std::vector<int> vertex_level(embedding.size());
    name_levels(g, embedding, vertex_level, k);

    int v = 0;

    for (int i = 0; i < vertex_level.size(); i++) {
        if (vertex_level[i] == 0){
            v = i;
            break;
        }
    }

    std::vector<Problem> tree = build_tree<Problem>(g, embedding, v, vertex_level);

//    std::vector<Problem> tree(my_vis.current_face);




//    for(auto v : embedding) {
//        copy.emplace_back();
//        for(auto e : v) {
//            copy.back().push_back(e);
//        }
//    }
//
//    for(int level = 1; level <= k; level++) {
////        property_map<Graph, edge_faces_t>::type faces = get(edge_faces_t(), g);
////        std::vector<Edge> outer_face;
////        my_visitor<Edge> my_vis(faces, outer_face);
////        planar_face_traversal(g, &copy.front(), my_vis);
////        std::vector<Problem> tree(my_vis.current_face);
////
////        for (Edge e : outer_face) {
////            vertex_level[e.m_source] = level;
////        }
//
//        int start = -1;
//
//        for (int i = 0; i < vertex_level.size(); i++) {
//            if (vertex_level[i] == 0) {
//                start = i;
//                break;
//            }
//        }
//
//        if (start == -1)
//            break;
//
//        int before = start;
//
//        int next = -1;
//
//        for (int i = 0; i < copy[start].size(); i++) {
//            Edge nextnext = copy[start][i];
//            if ((nextnext.m_source != start && vertex_level[nextnext.m_source] == 0)
//                || (nextnext.m_target != start && vertex_level[nextnext.m_target] == 0)) {
//                if (start == nextnext.m_target)
//                    next = nextnext.m_source;
//                else
//                    next = nextnext.m_target;
//                break;
//            }
//        }
//
//        if (next == -1) {
//            vertex_level[start] = level;
//            break;
//        }
//
//        while (start != next) {
//            vertex_level[next] = level;
//
//            for (int i = 0; i < copy[next].size(); i++) {
//                if (copy[next][i].m_source == before || copy[next][i].m_target == before) {
//                    before = next;
//                    for (int j = (i + 1)%copy[next].size(); j != i; (j++)%copy[next].size()) {
//                        Edge nextnext = copy[next][j];
//                        if (vertex_level[nextnext.m_source] == 0 || vertex_level[nextnext.m_target] == 0) {
//                            if (next == nextnext.m_target)
//                                next = nextnext.m_source;
//                            else
//                                next = nextnext.m_target;
//                        }
////                        break;
//                    }
//
////                    Edge nextnext = copy[next][(i+1)%copy[next].size()];
////                    if (next == nextnext.m_target)
////                        next = nextnext.m_source;
////                    else
////                        next = nextnext.m_target;
//                    break;
//                }
//            }
//
////            std::cout << next << "\n";
//        }
//        vertex_level[start] = level;
//
////        int i = 0;
////        for (auto v_it = copy.begin(); v_it != copy.end();) {
////            if (vertex_level[i] == level) {
////                v_it = copy.erase(v_it);
////            } else {
////                ++v_it;
////
////                for (auto e_it = v_it -> begin(); e_it != v_it -> end();) {
////                    if (vertex_level[e_it -> m_source ]== level
////                        || vertex_level[e_it -> m_target] == level) {
////                        e_it = v_it -> erase(e_it);
////                    } else {
////                        ++e_it;
////                    }
////                }
////            }
////            i++;
////        }
//    }

//    for (auto i : vertex_level) {
//        std::cout << i << "\n";
//    }

    return 0;
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
