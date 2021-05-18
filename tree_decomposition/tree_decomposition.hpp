//
// Created by mikolajtwarog on 2021-05-18.
//

#ifndef TECHNIKA_BAKER_TREE_DECOMPOSITION_HPP
#define TECHNIKA_BAKER_TREE_DECOMPOSITION_HPP

#include <vector>
#include <list>
#include <queue>
#include "../utils/cyclic_vector.hpp"

class tree_decomposition {

    std::vector< cyclic_vector<int> > embedding;
    std::vector< std::vector > tree;
    std::vector<std::vector > nodes;

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

        tree.emplace_back();
        nodes.emplace_back();
        int last_v = low_degree.front();
        nodes.back().push_back(last_v);
        nodes.back().push_back(graph[last_v].front());

        while (!deleted.empty()) {
            std::vector<int> front = deleted.front();
            deleted.pop()

            tree.last().push_back(tree.size());
            tree.emplace_back();
            tree.last().push_back(tree.size() - 2);

            nodes.emplace_back();
            for (int v : front) {
                nodes.back().push_back(v);
            }
        }

    }
};


#endif //TECHNIKA_BAKER_TREE_DECOMPOSITION_HPP
