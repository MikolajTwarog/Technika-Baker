//
// Created by mikolajtwarog on 2021-05-18.
//

#ifndef TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP
#define TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP

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
#include "../utils/cyclic_vector.hpp"

using namespace boost;

typedef subgraph<adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int>
        > > Graph;

typedef std::vector<cyclic_vector< graph_traits<Graph>::edge_descriptor > > PlanarEmbedding;

typedef graph_traits<Graph>::edge_descriptor Edge;
#include "../utils/visitors.hpp"
#include "../utils/level_face_traversal.hpp"
#include "../utils/find_outer_face.hpp"
#include "../utils/name_levels.hpp"

struct tree_decomposition{
    std::vector< std::vector<int> > tree;
    std::vector< std::set<int> > nodes;
    int n;

    std::map<Edge, int> edge_it;

    std::set<int>& operator[] (int x) {
        return nodes[x];
    }

    std::set<int>& operator[] (Edge e) {
        return nodes[edge_it[e]];
    }
};

class create_tree_decomposition {

    std::vector< cyclic_vector<Edge> > embedding;
    Graph graph;
    tree_decomposition tr;
    std::vector<int> vertex_level;
    std::vector< std::vector<Edge> > outer_edges;
    std::vector< cyclic_vector<int> > spanning_forest;
    int k;
    std::vector< std::vector< cyclic_vector<int> > > level_graphs;
    std::vector< std::vector< cyclic_vector<int> > > level_cycles;
    int outer_face;
    std::vector<int> expanded_vertices;

    std::vector<int> visited;
    std::vector<int> grey;
    int cycle_num = 0;
    int cycle_dfs(int v, int p ,int level) {
        visited[v] = 1;

        int cycle = -1;
        for (int nei : level_graphs[level][v]) {
            if (visited[nei] == 1 && nei != p) {
                grey.push_back(nei);
                cycle = cycle_num++;
                level_cycles[level].emplace_back();
            }
            if (visited[nei] == 0) {
                cycle = std::max(cycle, cycle_dfs(nei, v, level));
            }
        }

        visited[v] = 2;

        if (cycle >= 0) {
            level_cycles[level][cycle].push_back(v);
        }

        if (!grey.empty() && grey.back() == v) {
            cycle = false;
            grey.pop_back();
        }

        return cycle;
    }

    bool check_for_edge(int v, int w, std::vector< cyclic_vector<int> >& g) {
        for (int nei : g[v]) {
            if (nei == w) {
                return true;
            }
        }

        return false;
    }

    void get_leaves(int v, int p, int level, std::set<int>& leaves) {
        for (int w : spanning_forest[v]) {
            if (w != p) {
                get_leaves(w, v, level, leaves);
            }
        }

        if (vertex_level[v] == level) {
            leaves.insert(v);
        }
    }

    void outerplanar() {
        std::vector< std::list<int> > g(embedding.size());
        std::queue<int> low_degree;

        for (int i = 0; i < embedding.size(); i++) {
            if (embedding[i].size() <= 2) {
                low_degree.push(i);
            }

            for (Edge& edge : embedding[i]) {
                int v = edge.m_source == i ? edge.m_target : edge.m_source;
                g[i].push_back(v);
            }
        }

        std::queue< std::vector<int> > deleted;

        for (int i = 0; i < g.size() - 1; i++) {
            int v = low_degree.front();
            low_degree.pop();
            deleted.emplace();
            deleted.back().push_back(v);

            if (g[v].size() == 0) {
                continue;
            }

            if (g[v].size() == 1) {
                int nei = g[v].front();
                deleted.back().push_back(nei);

                for (auto it = g[nei].begin(); it != g[nei].end(); ) {
                    if (*it == v) {
                        it = g[nei].erase(it);
                    } else {
                        ++it;
                    }
                }

                if (g[nei].size() == 2) {
                    low_degree.push(nei);
                }

                continue;
            }

            int n_one = g[v].front();
            int n_two = g[v].back();
            deleted.back().push_back(n_one);
            deleted.back().push_back(n_two);
            bool edge_needed = true;

            for (auto it = g[n_one].begin(); it != g[n_one].end(); ) {
                if (*it == n_two) {
                    edge_needed = false;
                }
                if (*it == v) {
                    it = g[n_one].erase(it);
                } else {
                    ++it;
                }
            }

            for (auto it = g[n_two].begin(); it != g[n_two].end(); ) {
                if (*it == v) {
                    it = g[n_two].erase(it);
                } else {
                    ++it;
                }
            }

            if (edge_needed) {
                g[n_one].push_back(n_two);
                g[n_two].push_back(n_one);
            }

            if (g[n_one].size() == 2) {
                low_degree.push(n_one);
            }

            if (g[n_two].size() == 2) {
                low_degree.push(n_two);
            }
        }

        tr.tree.emplace_back();
        tr.nodes.emplace_back();
        int last_v = low_degree.front();
        tr.nodes.back().insert(last_v);
        tr.nodes.back().insert(g[last_v].front());

        while (!deleted.empty()) {
            std::vector<int> front = deleted.front();
            deleted.pop();

            tr.tree.back().push_back(tr.tree.size());
            tr.tree.emplace_back();
            tr.tree.back().push_back(tr.tree.size() - 2);

            tr.nodes.emplace_back();
            for (int v : front) {
                tr.nodes.back().insert(v);
            }
        }
    }

    void expand_vertices() {
        std::map< std::pair<int, int>, std::vector<int> > faces;
        std::vector< std::vector<int> > vertices_in_face;
        face_getter<Edge> my_vis(&faces, vertices_in_face);
        level_face_traversal<Graph>(embedding, my_vis);

        for (int j = 0; j < vertices_in_face.size(); j++) {
            auto& face = vertices_in_face[j];
            Edge current_e;
            Edge next_e(face[0], face[1], nullptr);

            bool res = true;
            for (int i = 1; i < face.size() - 1; i++) {
                current_e = next_e;
                next_e.m_source = face[i];
                next_e.m_target = face[i + 1];

                int dis = (get_edge_it(next_e, face[i], embedding) - get_edge_it(current_e, face[i], embedding)
                           + embedding[face[i]].size()) % embedding[face[i]].size();
                if (dis != 1 && dis != embedding[face[i]].size() - 1) {
                    res = false;
                    break;
                }
            }

            if (res) {
                outer_face = j;
                break;
            }
        }

        std::vector<int> face_level(vertices_in_face.size(), INT16_MAX);
        face_level[outer_face] = 0;

        std::queue<int> q;
        q.push(outer_face);

        while (!q.empty()) {
            int face = q.front();
            q.pop();

            int prev_v = vertices_in_face[face].back();
            for (int v : vertices_in_face[face]) {
                std::pair<int, int> e(std::min(prev_v, v), std::max(prev_v, v));
                std::vector<int>& fs = faces[e];
                int nei = fs[0] == face ? fs[1] : fs[0];

                if (face_level[nei] == INT16_MAX) {
                    q.push(nei);
                } else {
                    face_level[face] = std::min(face_level[face], face_level[nei] + 1);
                }

                prev_v = v;
            }
        }

        int size = embedding.size();
        for (int v = 0; v < size; v++) {
            if (embedding[v].size() <= 3) {
                continue;
            }

            int starting_e_it;
            int lowest_level = INT16_MAX;
            int face_num;

            Edge& prev_e = embedding[v].back();
            std::pair<int, int> prev_edge =std::make_pair(std::min(prev_e.m_source, prev_e.m_target),
                                                          std::max(prev_e.m_source, prev_e.m_target));
            for (int i = 0; i < embedding[v].size(); i++) {
                Edge& e = embedding[v][i];
                int w = e.m_source == v ? e.m_target : e.m_source;
                std::pair<int, int> edge(std::min(v, w), std::max(v, w));
                std::vector<int>& edge_fs = faces[edge];
                std::vector<int>& prev_edge_fs = faces[prev_edge];
                int face;
                for (int f1 : edge_fs) {
                    for (int f2 : prev_edge_fs) {
                        if (f1 == f2) {
                            face = f1;
                            break;
                        }
                    }
                }

                if (face_level[face] < lowest_level ) {
                    lowest_level = face_level[face];
                    face_num = face;
                    starting_e_it = i;
                }

                prev_edge = edge;
            }

            cyclic_vector<Edge> v_emb = embedding[v];
            embedding[v].clear();
            embedding[v].push_back(v_emb[starting_e_it]);
            embedding[v].push_back(v_emb[starting_e_it + 1]);
            embedding[v].emplace_back(v, embedding.size(), nullptr);
            int last = v;

            for (int i = starting_e_it + 2; i < v_emb.size() + starting_e_it - 2; i++) {
                int curr = embedding.size();
                embedding.emplace_back();
                embedding.back().emplace_back(last, curr, nullptr);

                int nei = v_emb[i].m_source == v ? v_emb[i].m_target : v_emb[i].m_source;
                embedding.back().emplace_back(curr, nei, nullptr);
                embedding[nei][get_edge_it(v_emb[i], nei, embedding)] = Edge(curr, nei, nullptr);

                embedding.back().emplace_back(curr, curr + 1, nullptr);
                last = curr;

                expanded_vertices.push_back(v);
            }

            int curr = embedding.size();
            embedding.emplace_back();
            embedding.back().emplace_back(last, curr, nullptr);

            int e_it = starting_e_it - 2;
            int nei = v_emb[e_it].m_source == v ? v_emb[e_it].m_target : v_emb[e_it].m_source;
            embedding.back().emplace_back(curr, nei, nullptr);
            embedding[nei][get_edge_it(v_emb[e_it], nei, embedding)] = Edge(curr, nei, nullptr);

            e_it = starting_e_it - 1;
            nei = v_emb[e_it].m_source == v ? v_emb[e_it].m_target : v_emb[e_it].m_source;
            embedding.back().emplace_back(curr, nei, nullptr);
            embedding[nei][get_edge_it(v_emb[e_it], nei, embedding)] = Edge(curr, nei, nullptr);

            expanded_vertices.push_back(v);
        }

        Graph new_graph;
        for (int v = 0; v < embedding.size(); v++) {
            auto& v_emb = embedding[v];
            for (int i = 0; i < v_emb.size(); i++) {
                Edge edge = v_emb[i];
                int nei = edge.m_source == v ? edge.m_target : edge.m_source;
                if (v < nei) {
                    Edge e = add_edge(edge.m_source, edge.m_target, new_graph).first;
                    embedding[v][i] = e;
                    embedding[nei][get_edge_it(edge, nei, embedding)] = e;
                }
            }
        }

        graph = new_graph;
    }

    void create_spanning_forest() {
        std::vector< std::vector<int> > vertices_on_level(k + 1);

        for (int v = 0; v < vertex_level.size(); v++) {
            vertices_on_level[vertex_level[v]].push_back(v);
        }

        for (int level = 1; level <= k; level++) {
            std::vector< cyclic_vector<int> >& level_graph = level_graphs[level];
            std::vector<Edge>& edges = outer_edges[level];

            for (auto &edge : edges) {
                int v = edge.m_source;
                int w = edge.m_target;

                if (!check_for_edge(v, w, level_graph)) {
                    level_graph[v].push_back(w);
                    level_graph[w].push_back(v);
                }
            }

            for (int v : vertices_on_level[level]) {
                for (auto& edge : embedding[v]) {
                    int w = edge.m_source == v ? edge.m_target : edge.m_source;
                    if (vertex_level[w] == level - 1) {
                        level_graph[v].push_back(w);
                        level_graph[w].push_back(v);
                        edges.emplace_back(v, w, nullptr);
                    }
                }
            }

            for (int& v : visited) {
                v = 0;
            }

            cycle_dfs(vertices_on_level[level][0], -1, level);
        }

        std::vector<bool> vertex_in_tree(embedding.size());


        for (int v : vertices_on_level.back()) {
            for (auto& edge : embedding[v]) {
                int w = edge.m_source == v ? edge.m_target : edge.m_source;
                if (vertex_level[v] == k && vertex_level[w] == k && !check_for_edge(v, w, level_graphs[k])) {
                    spanning_forest[v].push_back(w);
                    vertex_in_tree[v] = true;
                }
            }
        }

        // 0 = gamma, 1 = alpha, 2 = beta
        std::vector<int> category(embedding.size());
        for (int level = k; level > 0; level--) {
            for (int v : vertices_on_level[k]) {
                category[v] = 0;
            }

            std::vector< cyclic_vector<int> >& cycles = level_cycles[level];
            std::vector< cyclic_vector<int> >& level_graph = level_graphs[level];
            std::vector<Edge>& level_edges = outer_edges[level];

            for (auto& cycle : cycles) {
                for (int v : cycle) {
                    if (level_graph[v].size() == 3) {
                        category[v] = 1;
                    } else {
                        if (vertex_in_tree[v]) {
                            category[v] = 2;
                        }
//                        for (auto& edge : embedding[v]) {
//                            int w = edge.m_source == v ? edge.m_target : edge.m_source;
//                            if (vertex_in_tree[w]) {
//                                category[v] = 2;
//                                if (!vertex_in_tree[v]) {
//                                    spanning_forest[v].push_back(w);
//                                    spanning_forest[w].push_back(v);
//                                    vertex_in_tree[v] = true;
//                                }
//                                break;
//                            }
//                        }
                    }
                }
            }

            for (auto& edge : level_edges) {
                int v = edge.m_source;
                int w = edge.m_target;

                if (category[v] == 0 && category[w] == 0 && !check_for_edge(v, w, spanning_forest)) {
                    spanning_forest[v].push_back(w);
                    spanning_forest[w].push_back(v);
                    vertex_in_tree[v] = true;
                    vertex_in_tree[w] = true;
                }
            }

            for (auto& cycle : cycles) {
                std::vector< std::pair<int, int> > edges_to_add;

                int starting_alpha_it = -1;
                for (int i = 0; i < cycle.size(); i++) {
                    int v = cycle[i];

                    if (category[v] == 1) {
                        for (int w : level_graph[v]) {
                            if (w != cycle[i - 1] && w != cycle[i + 1]) {
                                int x = cycle[-1];
                                edges_to_add.emplace_back(v, w);
                                vertex_in_tree[v] = true;
                                vertex_in_tree[w] = true;
                            }
                        }

                        starting_alpha_it = i;
                    }
                }

                if (starting_alpha_it == -1) {
                    for (int i = 0; i < cycle.size(); i++) {
                        int v = cycle[i];

                        if (category[v] == 0) {
                            starting_alpha_it = i;
                        }
                    }
                }

                std::set<int> leaves;

                if (starting_alpha_it == -1) {
                    starting_alpha_it = 0;
                    get_leaves(cycle[0], -1, level, leaves);
                }

                int prev_v = cycle[starting_alpha_it];
                for (int v_it = starting_alpha_it + 1; v_it < starting_alpha_it + cycle.size(); v_it++) {
                    int v = cycle[v_it];

                    if (category[v] == 2 && leaves.find(v) != leaves.end()) {
                        prev_v = v;
                        continue;
                    }

                    if (category[v] == 2) {
                        get_leaves(v, -1, level, leaves);
                    }

                    edges_to_add.emplace_back(v, prev_v);
                    vertex_in_tree[v] = true;
                    vertex_in_tree[prev_v] = true;
                }

                for (auto& e : edges_to_add) {
                    spanning_forest[e.first].push_back(e.second);
                    spanning_forest[e.second].push_back(e.first);
                }
            }
        }
    }

    struct node{
        std::pair<int, Edge>  p;
        std::vector< std::pair<int, Edge> > children;
        std::vector<Edge> face;
    };

    void tr_dfs(int v, int p, std::vector<Edge>& edges, std::vector<node>& tree) {
        Edge missing_e;
        for (auto c : tree[v].children) {
            if (c.first == p) {
                missing_e = c.second;
            } else {
                std::vector<Edge> temp_edges;
                tr_dfs(c.first, v, temp_edges, tree);
                for (auto& e : temp_edges) {
                    edges.push_back(e);
                }
            }
        }

        if (p == -1) {
            return;
        }

        for (auto& e : tree[v].face) {
            edges.push_back(e);
        }

        int end = missing_e.m_source;
        int other_end = missing_e.m_target;

        for (auto& e : edges) {
            tr[e].insert(end);
            if (e.m_source != other_end) {
                tr[e.m_source].insert(end);
            }
            if (e.m_target != other_end) {
                tr[e.m_target].insert(end);
            }
        }
    }

    void build_tree_decomposition() {
        for (int v = 0; v < embedding.size(); v++) {
            tr.nodes.emplace_back();
            tr.tree.emplace_back();
            tr[v].insert(v);
        }

        std::map< std::pair<int, int>, std::vector<int> > faces;
        std::vector< std::vector<int> > vertices_in_face;
        face_getter<Edge> my_vis(&faces, vertices_in_face);
        level_face_traversal<Graph>(embedding, my_vis);

        std::vector<node> face_tree(vertices_in_face.size());

        for (int face = 0; face < vertices_in_face.size(); face++) {
            int prev_v = vertices_in_face[face].back();
            for (int v : vertices_in_face[face]) {
                if (check_for_edge(prev_v, v, spanning_forest)) {
                    face_tree[face].face.push_back(Edge(prev_v, v, nullptr));
                }
                prev_v = v;
            }
        }

        graph_traits<Graph>::edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
            int v = std::min(ei->m_source, ei->m_target);
            int w = std::max(ei->m_source, ei->m_target);
            std::vector<int>& fs = faces[std::make_pair(v, w)];

            if (!check_for_edge(v, w, spanning_forest)) {
                face_tree[fs[0]].children.push_back(std::make_pair(fs[1], *ei));
                face_tree[fs[1]].children.push_back(std::make_pair(fs[0], *ei));
            } else {
                int s = tr.nodes.size();
                tr.nodes.emplace_back();
                tr.tree.emplace_back();
                tr[s].insert(v);
                tr[s].insert(w);
                tr.edge_it[*ei] = s;
                tr.tree[s].push_back(v);
                tr.tree[v].push_back(s);
                tr.tree[s].push_back(w);
                tr.tree[w].push_back(s);
            }
        }

        std::vector<Edge> edges;
        tr_dfs(outer_face, -1, edges, face_tree);
    }

    void contract_vertices(){
        for (auto& node : tr.nodes) {
            std::vector<int> insert_list;
            for(auto it = node.begin(); it != node.end(); ) {
                if(expanded_vertices[*it] != -1) {
                    insert_list.push_back(expanded_vertices[*it]);
                    it = node.erase(it);
                } else {
                    ++it;
                }
            }

            for (int v : insert_list) {
                node.insert(v);
            }
        }
    }

public:
    create_tree_decomposition(Graph& g, std::vector< cyclic_vector<Edge> >& emb): graph(g), embedding(emb),
    expanded_vertices(embedding.size(), -1), vertex_level(embedding.size()) {
        k = name_levels(embedding, vertex_level, outer_edges);
        if (k == 1) {
            outerplanar();
            return;
        }

        expand_vertices();
        vertex_level.resize(embedding.size());
        visited.resize(embedding.size());
        spanning_forest.resize(embedding.size());
        k = name_levels(embedding, vertex_level, outer_edges);
        outer_edges.emplace_back();
        level_graphs = std::vector< std::vector< cyclic_vector<int> > >(k + 1, std::vector< cyclic_vector<int> >(embedding.size()));
        level_cycles.resize(k + 1);


        create_spanning_forest();

        build_tree_decomposition();

        contract_vertices();
    }

    tree_decomposition get_tree_decomposition() {
        return tr;
    }
};

void get_tree_decomposition(Graph& g, PlanarEmbedding& emb, tree_decomposition& tr) {
    create_tree_decomposition ctr(g, emb);
    tr = ctr.get_tree_decomposition();
}


#endif //TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP
