//
// Created by mikolajtwarog on 2021-04-29.
//

#ifndef TECHNIKA_BAKER_RANDOM_GRAPH_HPP
#define TECHNIKA_BAKER_RANDOM_GRAPH_HPP

#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/make_biconnected_planar.hpp>

//typedef adjacency_list
//        <
//                vecS,
//                vecS,
//                undirectedS,
//                property<vertex_index_t, int>,
//property<edge_index_t, int, EdgeProperty>
//>
//Graph;

Graph random_graph(int n, int m) {
    Graph g;
    boost::random::mt19937 gen;
    for (int i = 1; i < n; i++) {
        boost::random::uniform_int_distribution<> dist(0, i-1);
        boost::add_edge(i, dist(gen), g);
    }

    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        put(e_index, *ei, edge_count++);
    }

    std::vector<cyclic_vector< graph_traits<Graph>::edge_descriptor > > embedding(n);
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]);

    make_biconnected_planar(g, &embedding[0]);
    edge_count = 0;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        put(e_index, *ei, edge_count++);
    }
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                 boyer_myrvold_params::embedding =
                                         &embedding[0]);


    make_maximal_planar(g, &embedding[0]);

    std::cout << n << " " << m << std::endl;

//    boost::random::uniform_int_distribution<> dist(0, n-1);
//    while(num_edges(g) > m) {
//        remove_edge(dist(gen), dist(gen), g);
//    }

    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        std::cout << ei->m_source << " " << ei->m_target << "  ";
    }
    std::cout << std::endl;

    return g;
}

#endif //TECHNIKA_BAKER_RANDOM_GRAPH_HPP
