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
    ::tree<Problem> &tree;
    int current_face = 0;
    int last_vertex;
    Graph graph;
    PlanarEmbedding embedding;

    tree_builder(std::map<Edge, std::vector<int> >& f, ::tree<Problem>& t, Graph& g, PlanarEmbedding& emb)
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
        if(current_face != tree.outer_face) {
            int neighbor = faces[e][0] == current_face ? faces[e][1] : faces[e][0];

//            std::cout << "\n" << e << "\n";
//            for (int i : faces[e]) {
//                std::cout << i << "\n";
//            }

            if (neighbor == tree.outer_face) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
                tree[last].label.first = last_vertex;
                tree[last].label.second = last_vertex == e.m_source ? e.m_target : e.m_source;
//                tree[last].my_tree = &tree;

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

//    void root_tree_with_root(::tree<Problem> &tree, int node, int root);

    Graph g;
    PlanarEmbedding embedding;
    ::tree<Problem> tree;
    std::vector<int> vertex_level;
    std::set< std::pair<int, int> > added_edges;

    void find_outer_face(std::vector<int>& outer_face) {
        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        planar_face_traversal(g, &embedding.front(), my_vis);

        for (const auto& face : vertices_in_face) {
            Edge current_e;
            Edge next_e(face[0], face[1], nullptr);

            bool res = true;
            for (int i = 1; i < face.size() - 1; i++) {
                current_e = next_e;
                next_e.m_source = face[i];
                next_e.m_target = face[i + 1];

                int dis = (get_edge_it(next_e, face[i]) - get_edge_it(current_e, face[i])
                        + embedding[face[i]].size()) % embedding[face[i]].size();
                if (dis != 1 && dis != embedding[face[i]].size() - 1) {
                    res = false;
                    break;
                }

            }

            if (res) {
                for (int v : face) {
                    outer_face.push_back(v);
                }
                return;
            }
        }
    }

    int name_levels() {
//        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
//        std::vector<std::vector<int> > vertices_in_face;
//        my_visitor<Edge> my_vis(faces, vertices_in_face);
//        planar_face_traversal(g, &embedding.front(), my_vis);

        std::vector<int> outer_face;
        find_outer_face(outer_face);

        for (int &v : vertex_level) {
            v = -1;
        }

        for (int v : outer_face) {
            vertex_level[v] = 1;
        }

        std::queue<Edge> next_level_edges;

        for (int i = 0; i < outer_face.size(); i++) {
            int v = outer_face[i];
            int w = outer_face[(i + 1) % outer_face.size()];
            for (Edge e : embedding[v]) {
                if (vertex_level[e.m_source] == -1 || vertex_level[e.m_target] == -1) {
                    next_level_edges.push(e);
                }

                if (e.m_source == w) {
                    std::swap(e.m_source, e.m_target);
                    int e_it = get_edge_it(e, w);
                    std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
                }
            }
        }

        Edge next_level_edge;
        Edge current_edge;
        int starting_v;
        int current_v;
        int current_egde_it;


//        typename boost::graph_traits<Graph>::edge_iterator fi, fi_end;
//        for (boost::tie(fi, fi_end) = edges(g); fi != fi_end; ++fi) {
//            if (faces[*fi][0] == 0 || faces[*fi][1] == 0) {
//                next_level_edge = *fi;
//                current_edge = next_level_edge;
//                starting_v = next_level_edge.m_source;
//                current_v = next_level_edge.m_target;
//                vertex_level[current_v] = 0;
//                next_level_edges.push(next_level_edge);
//                break;
//            }
//        }

        int level = 1;
        std::vector<int> current_level;
//        current_level.push_back(starting_v);
//        current_level.push_back(current_v);

        while (!next_level_edges.empty()) {
            next_level_edge = next_level_edges.front();
            next_level_edges.pop();

            if (vertex_level[next_level_edge.m_source] > -1 && vertex_level[next_level_edge.m_target] > -1) {
                continue;
            }

            level = vertex_level[next_level_edge.m_source] == -1 ?
                    vertex_level[next_level_edge.m_target] : vertex_level[next_level_edge.m_source];
            level++;

            starting_v = vertex_level[next_level_edge.m_source] == -1 ?
                         next_level_edge.m_source : next_level_edge.m_target;

            current_level.push_back(starting_v);

            current_egde_it = get_edge_it(next_level_edge, starting_v);

            vertex_level[starting_v] = level;
            bool res = false;
            for (int j = (current_egde_it + 1) % embedding[starting_v].size(); j != current_egde_it;
                 j = (j + 1) % embedding[starting_v].size()) {
                Edge e_j = embedding[starting_v][j];
                if (vertex_level[e_j.m_source] == -1
                    || vertex_level[e_j.m_target] == -1) {
                    current_edge = e_j;
                    res = true;
                    break;
                }
            }
            if (!res) {
                continue;
            }

            vertex_level[starting_v] = -1;

            current_v = current_edge.m_source == starting_v ? current_edge.m_target : current_edge.m_source;
            current_egde_it = get_edge_it(current_edge, current_v);

            current_level.push_back(current_v);
            vertex_level[current_v] = level;

            while (current_v != starting_v) {
//                Edge e = embedding[current_v][current_egde_it];

                for (int j = (current_egde_it + 1) % embedding[current_v].size(); j != current_egde_it;
                     j = (j + 1) % embedding[current_v].size()) {
                    Edge e_j = embedding[current_v][j];
                    if (vertex_level[e_j.m_source] == -1
                        || vertex_level[e_j.m_target] == -1) {
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

            current_level.pop_back();

            for (int i = 0; i < current_level.size(); i++) {
                int v = current_level[i];
                int w = current_level[(i+1) % current_level.size()];
                for (Edge& e : embedding[v]) {
                    if (vertex_level[e.m_source] == -1
                        || vertex_level[e.m_target] == -1) {
                        next_level_edges.push(e);
                    }

                    if (e.m_source == w) {
                        std::swap(e.m_source, e.m_target);
                        int e_it = get_edge_it(e, w);
                        std::swap(embedding[w][e_it].m_source, embedding[w][e_it].m_target);
                    }
                }
            }
        }
        return level;
    }

    void triangulate(std::vector<int>& face, std::vector<int>& component, int turn) {
        int level = vertex_level[face[0]];
        Edge connecting_e;
        int starting_v;
        int starting_v_it = -1;
        int c_v;
        
        for (int i = 0; i < face.size(); i++) {
            int v = face[i];
            for (Edge e : embedding[v]) {
                int neighbour = e.m_source == v ? e.m_target : e.m_source;
                if (std::find(component.begin(), component.end(), neighbour) != component.end()) {
                    connecting_e = e;
                    starting_v = v;
                    starting_v_it = i;
                    c_v = neighbour;
                    break;
                }
            }
            if (starting_v_it > -1)
                break;
        }

        std::rotate(face.begin(), face.begin() + starting_v_it, face.end());
//        std::rotate(component.begin(), std::find(component.begin(), component.end(), c_v), component.end());
//
//        int face_it = 0;
//        int comp_it = 1;
//        int curr_e_it_face = get_edge_it(connecting_e, starting_v);
//        int curr_e_it_comp = get_edge_it(connecting_e, c_v);
//
//        while (face_it < face.size() && comp_it < component.size()) {
//            bool res = false;
//            int comp_curr = component[comp_it];
//            int face_curr = face[face_it];
//            int last = -1;
//
//            for (int i = curr_e_it_comp; i != (curr_e_it_comp + 1) % embedding[comp_curr].size();
//            i = (i - 1 + embedding[comp_curr].size()) % embedding[comp_curr].size()) {
//                Edge& e = embedding[comp_curr][i];
//                int neighbour = e.m_source == comp_curr ? e.m_target : e.m_source;
//                if (vertex_level[neighbour] == level && neighbour != face_curr) {
//                    res = true;
//                    last = neighbour;
//                    break;
//                }
//            }
//
//            if (res)  {
//                for (int i = face_it + 1; ; i++) {
//                    if (face[i] == last) {
//                        face_it = i;
//                        curr_e_it_face = get_edge_it(Edge(last, comp_curr, nullptr), last);
//                        comp_it++;
//                        curr_e_it_comp = get_edge_it(Edge(component[comp_it], comp_curr, nullptr), component[comp_it]);
//                        break;
//                    }
//
//                    Edge e(face[i - 1], face[i], nullptr);
//                    int e_it = get_edge_it(e, face[i]);
//                    embedding[face[i]].insert(embedding[face[i]].begin() + e_it + 1,
//                                              Edge(face[i], comp_curr, nullptr));
//                    embedding[comp_curr].insert(embedding[comp_curr].begin() +
//                    ((curr_e_it_comp - 1 + embedding[comp_curr].size()) % embedding[comp_curr].size()),
//                    Edge(face[i], comp_curr, nullptr));
//
//                    added_edges.emplace(face[i], comp_curr);
//
//                    add_edge(face[i], comp_curr, g);
//                }
//            } else {
//                embedding[face_curr].insert(embedding[face_curr].begin() + curr_e_it_face,
//                                            Edge(face_curr, comp_curr,nullptr));
//                embedding[comp_curr].insert(embedding[comp_curr].begin() +
//                                            ((curr_e_it_comp - 1 + embedding[comp_curr].size()) % embedding[comp_curr].size()),
//                                            Edge(face_curr, comp_curr, nullptr));
//
//                added_edges.emplace(face_curr, comp_curr);
//
//                add_edge(face_curr, comp_curr, g);
//
//                comp_it++;
//                curr_e_it_comp = get_edge_it(Edge(component[comp_it], comp_curr, nullptr), component[comp_it]);
//            }
//        }

        for (int i = 0; i < face.size(); i++) {
            int first = face[i];
            int second = face[(i + 1) % face.size()];

            Edge temp1(first, second, nullptr);

            int edge_it = get_edge_it(temp1, first);
            int edge_it2 = get_edge_it(temp1, second);


            edge_it = (edge_it + (1 * turn) + embedding[first].size()) % embedding[first].size();
            connecting_e = embedding[first][edge_it];

            int target = connecting_e.m_source == first ? connecting_e.m_target : connecting_e.m_source;
            Edge temp2(first, target, nullptr);

            if (!boost::edge(second, target, g).second) {
                add_edge(second, target, g);
                int target_e_it = get_edge_it(temp2, target);
                if (turn == -1) {
                    edge_it2++;
                } else {
                    target_e_it++;
                }
                embedding[second].emplace(embedding[second].begin() + edge_it2, second, target, nullptr);
                embedding[target].emplace(embedding[target].begin() + target_e_it, second, target, nullptr);
                added_edges.emplace(second, target);
            }
        }
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

            if (vertex_level[target_i] == level + 1) {
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

    void find_dividing_points(Edge one, Edge two, std::set<int>& dividing_points) {
        int v;
        if (one.m_source == two.m_source || one.m_source == two.m_target) {
            v = one.m_source == two.m_source ? two.m_source : two.m_target;
        } else {
            v = one.m_target == two.m_source ? two.m_source : two.m_target;
        }

        int level = vertex_level[v];

        int one_it = get_edge_it(one, v);
        int two_it = get_edge_it(two, v);

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

    void get_leaves(::tree<Problem>& t, std::vector<int>& leaves, int v) {
        if (t[v].children.empty()) {
            leaves.push_back(v);
            return;
        }

        for (auto child : t[v].children) {
            get_leaves(t, leaves, child);
        }
    }

    std::vector<bool> visited;

    void get_component(std::vector< std::vector<int> >& components, std::map<int, int>& vis,
                       std::vector< std::pair<int, int> >& v_in_c) {
        if (vertex_level[v_in_c[0].first] == 1) {
            find_outer_face(components[0]);

            return;
        }

        int comp_num = 0;

        for (int c = 0; c < v_in_c.size(); c++) {
            if (vis.find(v_in_c[c].first) != vis.end()) {
                continue;
            }
            components.emplace_back();

            int v = v_in_c[c].first;

            int level = vertex_level[v];
            int starting_v = v;
            int connecting_e_it = -1;

            for (int i = 0; i < embedding[v].size(); i++) {
                Edge e = embedding[v][i];
                int neighbour = e.m_source == v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level - 1) {
                    connecting_e_it = i;
                    break;
                }
            }

            vis[starting_v] = comp_num;
            components.back().push_back(starting_v);

            Edge current_e;
            int current_v;
            bool res = false;

            auto &edges = embedding[starting_v];
            for (int i = (connecting_e_it + 1) % edges.size(); i != connecting_e_it;
                 i = (i + 1) % edges.size()) {
                Edge e = edges[i];
                int neighbour = e.m_source == starting_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    res = true;
                    break;
                }
            }
            if (!res) {
                comp_num++;
                continue;
            }

            while (current_v != starting_v) {
                vis[current_v] = comp_num;
                components.back().push_back(current_v);

                int current_e_it = get_edge_it(current_e, current_v);
                auto &edges2 = embedding[current_v];
                int i = (current_e_it - 1 + edges2.size()) % edges2.size();
                res = false;
                for (; i != current_e_it; i = (i - 1 + edges2.size()) % edges2.size()) {
                    Edge e = edges2[i];
                    int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                    if (vertex_level[neighbour] == level) {
                        current_e = e;
                        current_v = neighbour;
                        res = true;
                        break;
                    }
                }

                if (res == false) {
                    Edge e = edges2[current_e_it];
                    int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                    if (vertex_level[neighbour] == level) {
                        current_e = e;
                        current_v = neighbour;
                    }
                }
            }

            comp_num++;
        }
    }

    void make_connected(std::vector< std::vector<int> >& components, std::map<int, int>& vis,
                        std::vector< std::pair<int, int> >& v_in_c) {
        for (int i = 0; i < v_in_c.size(); i++) {
            int curr = v_in_c[i].first;
            int curr_con = v_in_c[i].second;
            int next = v_in_c[(i + 1) % v_in_c.size()].first;
            int next_con = v_in_c[(i + 1) % v_in_c.size()].second;

            if (vis[curr] == vis[next] || edge(curr, next, g).second) {
                continue;
            }

            added_edges.insert(std::pair<int, int>(curr, next));
            add_edge(curr, next, g);

            embedding[curr].insert(embedding[curr].begin() + get_edge_it(Edge(curr, curr_con, nullptr), curr),
                                   Edge(curr, next, &i));

            embedding[next].insert(embedding[next].begin() + get_edge_it(Edge(next, next_con, nullptr), next) + 1,
                                   Edge(curr, next,&i));
        }
    }

    void check_for_components(const std::vector<int> &face, std::vector< std::pair<int, int> >& out) {
        if (face.size() == 2) {
            return;
        }

        int level = vertex_level[face[0]];

        int prev = face.size() - 1, curr = 0, next = 1;

        for (; curr < face.size(); prev = (prev + 1) % face.size(), curr++, next = (next + 1) % face.size()) {
            int one = get_edge_it(Edge(face[prev], face[curr], nullptr), face[curr]);
            int two = get_edge_it(Edge(face[curr], face[next], nullptr), face[curr]);
            for (int i = (two + 1) % embedding[face[curr]].size(); i != one;
            i = (i + 1) % embedding[face[curr]].size()) {
                Edge e = embedding[face[curr]][i];
                int target = e.m_source == face[curr] ? e.m_target : e.m_source;
                if (vertex_level[target] == level + 1 &&
                (out.empty() || (out.back().first != target && out.front().first != target))) {
                    out.emplace_back(target, face[curr]);
                }
            }
        }

        return;
    }

    template<typename Edge>
    int get_edge_it(Edge e, int v) {
        for (int i = 0; i < embedding[v].size(); i++) {
            if ((embedding[v][i].m_source == e.m_source
                 && embedding[v][i].m_target == e.m_target)
                || (embedding[v][i].m_source == e.m_target
                    && embedding[v][i].m_target == e.m_source)) {
                return i;
            }
        }

        return -1;
    }

    void root_tree_2(::tree<Problem> &t, int node) {
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
    }

    void root_tree_with_root(::tree<Problem> &t, int root) {
        int node;

        for (int i = 0; i < t.size(); i++) {
            for (int v : t[i].face) {
                if (v == root) {
                    node = i;
                    break;
                }
            }
        }

        t.root = node;

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

    ::tree<Problem> build_tree(std::vector<int> component, std::map<int, std::vector<Edge> > emb) {
        int level = vertex_level[component[0]];

        std::reverse(component.begin(), component.end());

        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        level_face_traversal<Graph>(emb, my_vis, level, vertex_level, component);
        ::tree<Problem> t(my_vis.current_face, level);

        int outer = 0;
        for (auto& face : vertices_in_face) {
            int start = -1;
            for (int i = 0; i < face.size(); i++) {
                if (face[i] == component[0]) {
                    start = i;
                    break;
                }
            }
            if (start == -1) {
                continue;
            }

            std::rotate(face.begin(), face.begin() + start, face.end());

            if (face == component) {
                break;
            }
            outer++;
        }

        t.outer_face = outer;
        tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, t, g, embedding);
        level_face_traversal<Graph>(emb, tree_b, level, vertex_level, component);

        t.remove_outer_face();

        return t;
    }

    ::tree<Problem> build_tree_with_dividing_points(std::vector<int> component, int root) {
        int level = vertex_level[component[0]];

        if (component.size() == 1) {
            ::tree<Problem> t(1, level);
            t[0].label.first = component[0];
            t[0].label.second = component[0];
            t.root = 0;
            return t;
        }

        std::set<int> vis;
        Graph g_temp;
        for (int v : component) {
            if (vis.find(v) != vis.end()) {
                continue;
            }

            for (Edge e : embedding[v]) {
                int target = e.m_source == v ? e.m_target : e.m_source;
                if (vertex_level[target] == level && v < target) {
                    add_edge(e.m_source, e.m_target, g_temp);
                }
            }

            vis.insert(v);
        }

        auto bicomp = get(edge_index, g_temp);
        int bi_num = biconnected_components(g_temp, bicomp);
        std::vector< ::tree<Problem> > trees(bi_num);
        std::vector< std::set<int> > v_in_bicomps(bi_num);

        std::map<std::pair<int, int>, int> bi_map;
        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(g_temp); ei != ei_end; ++ei) {
            v_in_bicomps[bicomp[*ei]].insert(ei->m_source);
            v_in_bicomps[bicomp[*ei]].insert(ei->m_target);
            bi_map[std::pair<int, int>(ei->m_source, ei->m_target)] = bicomp[*ei];
            bi_map[std::pair<int, int>(ei->m_target, ei->m_source)] = bicomp[*ei];
        }

        std::vector< std::vector<int> > bicomps(bi_num);
        for (int i = 0; i < component.size(); i++) {
            int curr = component[i];
            int next = component[(i + 1) % component.size()];
            bicomps[bi_map[std::pair<int, int>(curr, next)]].push_back(curr);
        }

        for (auto& bic : bicomps) {
            if (bic.size() == 2) {
                add_edge(bic[0], bic[1], g);
                Edge e(bic[0], bic[1], nullptr);
                int it = get_edge_it(e, bic[0]);
                embedding[bic[0]].insert(embedding[bic[0]].begin() + it, e);
                it = get_edge_it(e, bic[1]);
                embedding[bic[1]].insert(embedding[bic[1]].begin() + it, e);
            }
        }

        for (int i = 0; i < bi_num; i++) {
            std::map<int, std::vector<Edge> > emb;
            for (int v : v_in_bicomps[i]) {
                for (Edge e : embedding[v]) {
                    if (bi_map[std::pair<int, int>(e.m_source, e.m_target)] == i) {
                        emb[v].push_back(e);
                    }
                }
            }

            trees[i] = build_tree(bicomps[i], emb);
        }

        std::vector<int> art_points;
        articulation_points(g, std::back_inserter(art_points));

        std::map<int, std::vector<int> > art_to_bicomp;
        std::map<int, std::vector<int> > bicomp_to_art;

        for (int a : art_points) {
            for (int i = 0; i < v_in_bicomps.size(); i++) {
                if (v_in_bicomps[i].find(a) != v_in_bicomps[i].end()) {
                    art_to_bicomp[a].push_back(i);
                    bicomp_to_art[i].push_back(a);
                }
            }
        }

        int main_bicomp = 0;
        for (int i = 0; i < v_in_bicomps.size(); i++) {
            if (v_in_bicomps[i].find(root) != v_in_bicomps[i].end()) {
                main_bicomp = i;
                break;
            }
        }

        root_tree_with_root(trees[main_bicomp], root);

        std::vector<int> merged(bi_num);
        for (int i = 0; i < bi_num; i++) {
            merged[i] = i;
        }

        std::queue<int> q;
        q.push(main_bicomp);

        while(!q.empty()) {
            int cur_bicomp = q.front();
            q.pop();

            for (int a : bicomp_to_art[cur_bicomp]) {
                for (int b : art_to_bicomp[a]) {
                    if (merged[b] == b && b != main_bicomp) {
                        q.push(b);
                        root_tree_with_root(trees[b], a);

                        int target = cur_bicomp;
                        while (target != merged[target]) {
                            merged[target] = merged[merged[target]];
                            target = merged[target];
                        }

                        merged[b] = target;
                        trees[target].merge(trees[b]);
                    }
                }
            }
        }


        for (int i = 0; i < trees[main_bicomp].size(); i++) {
            auto& t = trees[main_bicomp][i];

            if (t.children.empty()) {
                continue;
            }

            auto& face = t.face;

            std::vector< std::pair<int, int> > v_in_c;
            check_for_components(face, v_in_c);
            if (!v_in_c.empty()) {
                std::vector< std::vector<int> > components;
                std::map<int, int> v_to_c;
                get_component(components, v_to_c, v_in_c);

                if (components.size() > 1) {
                    make_connected(components, v_to_c, v_in_c);
                    components.clear();
                    v_to_c.clear();
                    get_component(components, v_to_c, v_in_c);
                }

                triangulate(face, components[0], 1);
                triangulate(components[0], face, -1);
                int v = find_third(t.label.first, t.label.second);
                t.component_tree = build_tree_with_dividing_points(components[0], v);
                t.component_tree.enclosing_tree = &trees[main_bicomp];
                t.component_tree.enclosing_face = i;
            }
        }

        return trees[main_bicomp];
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

            t[v].adjust(g, added_edges);
        } else if (!t[v].children.empty() && !t[v].component_tree.empty()) {
            auto& ct = t[v].component_tree;
            table(ct, ct.root);
            t[v].contract(ct[ct.root]);
            t[v].adjust(g, added_edges);
        } else if (level > 1) {
            std::vector<int> z_table;
            z_table.push_back(t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[0]].label.first);
            for (int x : t.enclosing_tree->t[t.enclosing_face].children) {
                z_table.push_back(t.enclosing_tree->t[x].label.second);
            }

            int p = t[v].LB;
            for (int i = t[v].RB; i > t[v].LB; i--) {
                if (check_for_edge(t[v].label.first, z_table[i], g, added_edges)) {
                    p = i;
                    break;
                }
            }

            t[v].create(p, g, added_edges);
            t[v].adjust(g, added_edges);

            int j = p - 1;

            while (j >= t[v].LB) {
                table(*t.enclosing_tree, t.enclosing_tree->t[t.enclosing_face].children[j]);
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.first, g, added_edges);
                t[v].merge(second.val, t[v].val);
                j--;
            }

            j = p;

            while (j < t[v].RB) {
                table(*t.enclosing_tree, t.enclosing_tree->t[t.enclosing_face].children[j]);
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.second, g, added_edges);
                t[v].merge(t[v].val, second.val);
                j++;
            }
        }
    }

public:
    baker_impl(const Graph& arg_g): g(arg_g), embedding(num_vertices(arg_g)), vertex_level(num_vertices(arg_g)),
    visited(num_vertices(arg_g)) {
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
            if (vertex_level[i] == 1) {
                v = i;
                break;
            }
        }

        std::vector<int> face;
        find_outer_face(face);
        tree = build_tree_with_dividing_points(face, v);

        create_boudaries(tree, tree, tree.root);

        table(tree, tree.root);
    }

    int result() {
        return tree[tree.root].result();
    }

};

template <typename Problem>
int baker2(const Graph& g) {
    baker_impl < Problem > b(g);
    return b.result();
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
