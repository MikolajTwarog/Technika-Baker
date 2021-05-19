//
// Created by mikolajtwarog on 2021-05-18.
//

#ifndef TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP
#define TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP

#include <vector>
#include <list>
#include <queue>
#include <set>
#include "../utils/cyclic_vector.hpp"
#include "../utils/name_levels.hpp"

typedef subgraph<adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
property<edge_index_t, int>
> >
Graph;

typedef graph_traits<Graph>::edge_descriptor Edge;

struct tree_decomposition{
    std::vector< std::vector<int> > tree;
    std::vector< std::vector<int> > nodes;
};

class create_tree_decomposition {

    std::vector< cyclic_vector<int> > embedding;
    Graph graph;
    tree_decomposition tr;
    std::vector<int> vertex_level;
    std::vector< std::vector<Edge> > outer_edges;
    std::vector< std::vector<int> > spanning_forest;
    int k;
    std::vector< std::vector< std::vector<int> > > level_graphs;
    std::vector< std::vector< std::cyclic_vector<int> > > level_cycles;

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
            }
            if (visited[nei] == 0) {
                cycle = std::max(cycle, cycle_dfs(nei, v, level));
            }
        }

        visited[v] = 2;

        if (cycle >= 0) {
            level_cycles[level][cycle].push_back(v);
        }

        if (grey.back() == v) {
            cycle = false;
            grey.pop_back();
        }

        return cycle;
    }

    bool check_for_edge(int v, int w, int level) {
        for (int nei : level_graphs[level][v]) {
            if (nei == w) {
                return true;
            }
        }

        return false;
    }

    void get_leaves(int v, int p, int level, std::set<int>& leaves) {
        for (int w : spanning_forest[v]) {
            if (w != p) {
                get_leaves(w, v, leaves);
            }
        }

        if (vertex_level[v] == level) {
            leaves.insert(v);
        }
    }

    void outerplanar() {
        std::vector< std::list<int> > graph;
        std::queue<int> low_degree;

        for (int i = 0; i < embedding.size(); i++) {
            if (embedding[i].size() <= 2) {
                low_degree.push(i);
            }

            for (int v : embedding[i]) {
                graph[i].push_back(v);
            }
        }

        std::queue< std::vector<int> > deleted;

        for (int i = 0; i < graph.size() - 1; i++) {
            int v = low_degree.front();
            low_degree.pop();
            deleted.emplace();
            deleted.back().push_back(v);

            if (graph[v].size() == 0) {
                continue;
            }

            if (graph[v].size() == 1) {
                int nei = graph[v].front();
                deleted.back().push_back(nei);

                for (auto it = graph[nei].begin(); it != graph[nei].end(); ) {
                    if (*it == v) {
                        it = graph[nei].erase(it);
                    } else {
                        ++it;
                    }
                }

                if (graph[nei].size() == 2) {
                    low_degree.push(nei);
                }

                continue;
            }

            int n_one = graph[v].front();
            int n_two = graph[v].back();
            deleted.back().push_back(n_one);
            deleted.back().push_back(n_two);
            bool edge_needed = true;

            for (auto it = graph[n_one].begin(); it != graph[n_one].end(); ) {
                if (*it == n_two) {
                    edge_needed = false;
                }
                if (*it == v) {
                    it = graph[n_one].erase(it);
                } else {
                    ++it;
                }
            }

            for (auto it = graph[n_two].begin(); it != graph[n_two].end(); ) {
                if (*it == v) {
                    it = graph[n_two].erase(it);
                } else {
                    ++it;
                }
            }

            if (edge_needed) {
                graph[n_one].push_back(n_two);
                graph[n_two].push_back(n_one);
            }

            if (graph[n_one].size() == 2) {
                low_degree.push(n_one);
            }

            if (graph[n_two].size() == 2) {
                low_degree.push(n_two);
            }
        }

        tr.tree.emplace_back();
        tr.nodes.emplace_back();
        int last_v = low_degree.front();
        tr.nodes.back().push_back(last_v);
        tr.nodes.back().push_back(graph[last_v].front());

        while (!deleted.empty()) {
            std::vector<int> front = deleted.front();
            deleted.pop()

            tr.tree.last().push_back(tree.size());
            tr.tree.emplace_back();
            tr.tree.last().push_back(tree.size() - 2);

            tr.nodes.emplace_back();
            for (int v : front) {
                tr.nodes.back().push_back(v);
            }
        }
    }

    void expand_vertices(){}

    void create_spanning_forest() {
        for (int level = 1; level <= k; level++) {
            std::vector<Edge>& edges = outer_edges[level];
            std::vector< std::vector<int> >& level_graph = level_graphs[level];

            for (auto& edge : edges) {
                int v = edge.m_source;
                int w = edge.m_target;

                if (check_for_edge(v, w, level)) {
                    level_graph[v].push_back(w);
                    level_graph[w].push_back(v);
                }
            }

            cycle_dfs(edges[0].m_source, -1, level);
        }

        std::vector< std::vector<int> > vertices_on_level(k + 1);

        for (int v = 0; v < vertex_level.size(); v++) {
            vertices_on_level[vertex_level[v]].push_back(v);
        }

        std::vector<bool> vertex_in_tree(embedding.size());

        for (int v : vertices_on_level.back()) {
            for (auto& edge : embedding[v]) {
                int w = edge.m_sorce == v ? edge.m_target : edge.m_source;
                if (vertex_level[v] == k && vertex_level[w] == k && check_for_edge(v, w, k)) {
                    spanning_forest[v].push_back(w);
                    vertex_in_tree[v] = true;
                }
            }
        }

        // 0 = gamma, 1 = alpha, 2 = beta
        std::vector<int> category(embedding.size());
        for (int level = k; level > 0; level--) {
            std::vector< std::cyclic_vector<int> >& cycles = level_cycles[level];
            std::vector< std::vector<int> >& level_graph = level_graphs[level];
            std::vector<Edge>& level_edges = outer_edges[level];

            for (auto& cycle : cycles) {
                for (int v : cycle) {
                    if (level_graph[v].size() == 3) {
                        category[v] = 1;
                    } else {
                        for (auto& edge : embedding[v]) {
                            int w = edge.m_sorce == v ? edge.m_target : edge.m_source;
                            if (vertex_in_tree[w]) {
                                category[v] = 2;
                                if (!vertex_in_tree[v]) {
                                    spanning_forest[v].push_back(w);
                                    spanning_forest[w].push_back(v);
                                    vertex_in_tree[v] = true;
                                }
                                break;
                            }
                        }
                    }
                }
            }

            for (auto& edge : level_edges) {
                int v = edge.m_source;
                int w = edge.m_target;

                if (category[v] == 0 && category[w] == 0) {
                    spanning_forest[v].push_back(w);
                    spanning_forest[w].push_back(v);
                    vertex_in_tree[v] = true;
                    vertex_in_tree[w] = true;
                }
            }

            for (auto& cycle : cycles) {
                int starting_alpha_it;
                for (int i = 0; i < cycle.size(); i++) {
                    int v = cycle[i];

                    if (category[v] == 1) {
                        for (int w : level_graph[v]) {
                            if (w != cycle[i - 1] && w != cycle[i + 1]) {
                                spanning_forest[v].push_back(w);
                                spanning_forest[w].push_back(v);
                                vertex_in_tree[v] = true;
                                vertex_in_tree[w] = true;
                            }
                        }

                        starting_alpha_it = i;
                    }
                }

                std::set<int> leaves;

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

                    spanning_forest[v].push_back(prev_v);
                    spanning_forest[prev_v].push_back(v);
                    vertex_in_tree[v] = true;
                    vertex_in_tree[prev_v] = true;
                }
            }
        }
    }

    create_tree_decomposition(Graph& g, std::vector< cyclic_vector<int> >& emb): graph(g), embedding(emb),
    visited(embedding.size()) {
        k = name_levels(graph, embedding, vertex_level, outer_edges);
        level_graphs.resize(k, std::vector<int>(embedding.size()));
        level_cycles.resize(k);

        if (k == 1) {
            outerplanar();
            return;
        }

        expand_vertices();
        create_spanning_forest();
    }
};


#endif //TECHNIKA_BAKER_CREATE_TREE_DECOMPOSITION_HPP
