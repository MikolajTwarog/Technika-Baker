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
    std::vector<Edge> added_edges;

    struct coord_t
    {
        std::size_t x;
        std::size_t y;
    };


    void find_outer_face(std::vector<int>& face) {
        typedef std::vector< std::vector< graph_traits< Graph >::edge_descriptor > >
                embedding_storage_t;
        typedef boost::iterator_property_map< embedding_storage_t::iterator,
                property_map< Graph, vertex_index_t >::type >
                embedding_t;

        Graph g_copy(g);

        make_maximal_planar(g_copy, &embedding[0]);

        embedding_storage_t embedding_storage(num_vertices(g));
        embedding_t embedding_copy(embedding_storage.begin(), get(vertex_index, g));

        boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g_copy,
                                     boyer_myrvold_params::embedding = embedding_copy);

        std::vector< graph_traits< Graph >::vertex_descriptor > ordering;
        planar_canonical_ordering(g_copy, embedding_copy, std::back_inserter(ordering));
        typedef std::vector< coord_t > straight_line_drawing_storage_t;
        typedef boost::iterator_property_map<
                typename straight_line_drawing_storage_t::iterator,
                property_map< Graph, vertex_index_t >::type >
                straight_line_drawing_t;

        straight_line_drawing_storage_t straight_line_drawing_storage(
                num_vertices(g));
        straight_line_drawing_t straight_line_drawing(
                straight_line_drawing_storage.begin(), get(vertex_index, g_copy));

        chrobak_payne_straight_line_drawing(
                g_copy, embedding_copy, ordering.begin(), ordering.end(), straight_line_drawing);

        auto vertexI = get(vertex_index, g_copy);
        int left_v;
        int max_left = INT_MAX;
        graph_traits< Graph >::vertex_iterator vi, vi_end;
        std::vector<coord_t> coords(num_vertices(g));

        for (boost::tie(vi, vi_end) = vertices(g_copy); vi != vi_end; ++vi) {
            coord_t coord(get(straight_line_drawing, *vi));
            coords[vertexI[*vi]] = coord;
            if (coord.x < max_left) {
                max_left = coord.x;
                left_v = vertexI[*vi];
            }
        }

        coord_t left_v_coord = coords[left_v];
        Edge outer_edge;
        double biggest_cot = INT_MAX;
        int current_v;

        for (Edge e : embedding[left_v]) {
            int neighbour = e.m_source == left_v ? e.m_target : e.m_source;
            coord_t n_coord = coords[neighbour];
            coord_t vec;
            vec.x = n_coord.x - left_v_coord.x;
            vec.y = n_coord.y - left_v_coord.y;

            if (vec.x == 0) {
                if (vec.y > 0) {
                    outer_edge = e;
                    break;
                }
                continue;
            }

            double cot = vec.y / vec.x;

            if (cot < biggest_cot) {
                outer_edge = e;
                biggest_cot = cot;
                current_v = neighbour;
            }
        }

        int current_edge_it = get_edge_it(outer_edge, current_v);
        Edge current_edge = outer_edge;

        face.push_back(current_v);

        std::vector<bool> visited(num_vertices(g));

        while (current_v != left_v) {
            for (int j = (current_edge_it + 1) % embedding[current_v].size(); j != current_edge_it;
                j = (j + 1) % embedding[current_v].size()) {
                Edge e_j = embedding[current_v][j];
                if (!visited[e_j.m_source] || !visited[e_j.m_target]) {
                    current_edge = e_j;
                    break;
                }
            }

            current_v = current_edge.m_source == current_v ? current_edge.m_target : current_edge.m_source;
            current_edge_it = get_edge_it(current_edge, current_v);

            face.push_back(current_v);
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
            for (int j = (current_egde_it + 1) % embedding[starting_v].size(); j != current_egde_it;
                 j = (j + 1) % embedding[starting_v].size()) {
                Edge e_j = embedding[starting_v][j];
                if (vertex_level[e_j.m_source] == -1
                    || vertex_level[e_j.m_target] == -1) {
                    current_edge = e_j;
                    break;
                }
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
        
        for (int i = 0; i < face.size(); i++) {
            int v = face[i];
            for (Edge e : embedding[v]) {
                int neighbour = e.m_source == v ? e.m_target : e.m_source;
                if (std::find(component.begin(), component.end(), neighbour) != component.end()) {
                    connecting_e = e;
                    starting_v = v;
                    starting_v_it = i;
                    break;
                }
            }
            if (starting_v_it > -1)
                break;
        }

        std::rotate(face.begin(), face.begin() + starting_v_it, face.end());

        for (int i = 0; i < face.size(); i++) {
            int first = face[i];
            int second = face[(i + 1) % face.size()];

            Edge temp1(first, second, nullptr);

            int edge_it = get_edge_it(temp1, first);
            int edge_it2 = get_edge_it(temp1, second);

            
            edge_it = (edge_it + (1 * turn)) % embedding[first].size();
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
                added_edges.emplace_back(second, target, nullptr);
            }
        }

//        for (int i = 0; i < component.size(); i++) {
//            int v = component[i];
//            int prior = component[(i - 1) % component.size()];
//            int next = component[(i + 1) % component.size()];
//            
//            
//        }
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

        if (!t[node].component_tree.empty()) {
            std::vector<int> face;
            for (int x : t[node].children) {
                face.push_back(t[x].label.first);
            }
            if (t[t[node].children.back()].label.second != face[0]) {
                face.push_back(t[t[node].children.back()].label.second);
            }
            std::vector<int> component;
            get_component(component, check_for_component(face));
//            std::reverse(component.begin(), component.end());
            triangulate(face, component, -1);
            triangulate(component, face, 1);
            int v = find_third(t[node].label.first, t[node].label.second);
            root_tree_with_root(t[node].component_tree, v);
        }
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

    std::vector<bool> visited;

    void get_component(std::vector<int>& component, int v) {
        if (vertex_level[v] == 1) {
            std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
            std::vector<std::vector<int> > vertices_in_face;
            my_visitor<Edge> my_vis(faces, vertices_in_face);
            planar_face_traversal(g, &embedding.front(), my_vis);

            for (int v : vertices_in_face[0]) {
                component.push_back(v);
            }

            return;
        }

        std::queue<int> q;
        q.push(v);
        int level = vertex_level[v];
        int starting_v = -1;
        int connecting_e_it = -1;

        while (!q.empty()) {
            int current_v = q.front();
            q.pop();
//            component.push_back(current_v);

            for (int i = 0; i < embedding[current_v].size(); i++) {
                Edge e = embedding[current_v][i];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (!visited[neighbour] && vertex_level[neighbour] == level){
                    q.push(neighbour);
                    visited[neighbour] = true;
                }
                
                if (vertex_level[neighbour] == level - 1) {
                    starting_v = current_v;
                    connecting_e_it = i;
                    break;
                }
            }
            
            if (starting_v > -1) {
                break;
            }
        }

        component.push_back(starting_v);

        Edge current_e;
        int current_v;

        auto &edges = embedding[starting_v];
        for (int i = (connecting_e_it - 1 + + edges.size()) % edges.size(); i != connecting_e_it;
            i = (i - 1 + + edges.size()) % edges.size()) {
            Edge e = edges[i];
            int neighbour = e.m_source == starting_v ? e.m_target : e.m_source;

            if (vertex_level[neighbour] == level) {
                current_e = e;
                current_v = neighbour;
                break;
            }
        }

        while (current_v != starting_v) {
            component.push_back(current_v);

            int current_e_it = get_edge_it(current_e, current_v);
            auto& edges2 = embedding[current_v];
            int i = (current_e_it - 1 + edges2.size()) % edges2.size();
            for (;i != current_e_it; i = (i - 1 + edges2.size()) % edges2.size()) {
                Edge e = edges2[i];
                int neighbour = e.m_source == current_v ? e.m_target : e.m_source;

                if (vertex_level[neighbour] == level) {
                    current_e = e;
                    current_v = neighbour;
                    break;
                }
            }
        }

//        for (int i = 0, j = component.size() - 1; i < j; i++, j--) {
//            std::swap(component[i], component[j]);
//        }
    }

    ::tree<Problem> build_tree(int v) {
        int level = vertex_level[v];

        std::vector<int> component;
        get_component(component, v);

        if (level > 1) {
            std::reverse(component.begin(), component.end());
        }

        std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
        std::vector<std::vector<int> > vertices_in_face;
        my_visitor<Edge> my_vis(faces, vertices_in_face);
        level_face_traversal(g, embedding, my_vis, level, vertex_level, component);
        ::tree<Problem> t(my_vis.current_face, level);

//        for (Problem p : t.t) {
//            p.my_tree = &t;
//        }

        tree_builder<graph_traits<Graph>::edge_descriptor, Problem, PlanarEmbedding> tree_b(faces, t, g, embedding);
        level_face_traversal(g, embedding, tree_b, level, vertex_level, component);

        for (int i = 1; i < vertices_in_face.size(); i++) {
            auto &face = vertices_in_face[i];

            int v_in_c = check_for_component(face);
            if (v_in_c > -1) {
                t[i].component_tree = build_tree(v_in_c);
                t[i].component_tree.enclosing_tree = &t;
                t[i].component_tree.enclosing_face = i;
            }
        }

        return t;
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

        for (int i = 1; i < t.size(); i++) {
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
        } else if (level > 0) {
            std::vector<int> z_table;
            z_table.push_back(t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[0]].label.first);
            for (int x : t.enclosing_tree->t[t.enclosing_face].children) {
                z_table.push_back(t.enclosing_tree->t[x].label.second);
            }

            int p = t[v].RB;
            for (int i = t[v].LB; i < t[v].RB; i++) {
                if (boost::edge(t[v].label.second, z_table[i], g).second) {
                    p = i;
                    break;
                }
            }

            t[v].create(p, g);

            int j = p - 1;

            while (j >= t[v].LB) {
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.first, g);
                t[v].merge(second.val, t[v].val);
                j--;
            }

            j = p;

            while (j < t[v].RB) {
                Problem second = t.enclosing_tree->t[t.enclosing_tree->t[t.enclosing_face].children[j]].extend(t[v].label.second, g);
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

        tree = build_tree(v);

        root_tree_with_root(tree, 0);

        create_boudaries(tree, tree, 1);

        table(tree, 1);
    }

    int result() {
        return tree[0].val[0];
    }

};

template <typename Problem>
int baker2(const Graph& g) {
    baker_impl < Problem > b(g);
    return b.result();
}

#define TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP

#endif //TECHNIKA_BAKER_BAKER_K_OUTER_PLANAR_HPP
