//
// Created by mikolajtwarog on 2021-04-29.
//


#ifndef TECHNIKA_BAKER_BAKERS_TECHNIQUE_HPP
#define TECHNIKA_BAKER_BAKERS_TECHNIQUE_HPP

#include <boost/graph/connected_components.hpp>
#include "../k-outerplanar/baker-k-outer-planar.hpp"

int bakers_technique(Graph g, int k) {
    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);

    std::vector<cyclic_vector< graph_traits<Graph>::edge_descriptor > > embedding(num_vertices(g));
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )
            )
        std::cout << "Input graph is planar" << std::endl;
    else {
        std::cout << "Input graph is not planar" << std::endl;
        return -1;
    }

    std::vector<int> vertex_level(num_vertices(g));

    std::vector< std::vector<Edge> > outer_edges;
    std::vector<int> outer_face; // to do
    int max_level = name_levels(embedding, outer_face, vertex_level, outer_edges);

    std::vector<std::vector<int> > levels(max_level + 1);

    for (int i = 0; i < vertex_level.size(); i++) {
        levels[vertex_level[i]].push_back(i);
    }

    int res = 0;

    for (int i = 1; i <= max_level; i += (k + 1)) {
        Graph& sub_g = g.create_subgraph();
        for (int j = i; j < std::min(max_level + 1, i + k + 1); j++) {
            for (int v : levels[j]) {
                add_vertex(v, sub_g);
            }
        }

        Graph sub_g2;

        for (auto ei2 : sub_g.m_local_edge)
            add_edge(ei2.second.m_source, ei2.second.m_target, sub_g2);

        if (num_vertices(sub_g2) == 0) {
            res++;
            continue;
        }

        std::vector<int> v_comp(num_vertices(sub_g2));
        int num = connected_components(sub_g2, &v_comp[0]);

        std::vector<std::vector<int> > components(num);

        for (int v = 0; v < num_vertices(sub_g2); v++) {
            components[v_comp[v]].push_back(v);
        }

        for (int c = 0; c < num; c++) {
            Graph graph = sub_g2.create_subgraph(components[c].begin(), components[c].end());

            if (num_vertices(graph) == 1) {
                res++;
            } else {
                Graph graph2;

                for (auto edge : graph.m_local_edge)
                    add_edge(edge.second.m_source, edge.second.m_target, graph2);

//                res += baker2<independent_set>(graph2);
            }
        }
    }

    return res;
}


#endif //TECHNIKA_BAKER_BAKERS_TECHNIQUE_HPP
