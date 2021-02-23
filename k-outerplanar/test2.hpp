#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "baker-k-outer-planar.hpp"
#include "problems2.hpp"

using namespace boost;

typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

typedef std::vector< graph_traits<Graph>::edge_descriptor > vec_t;

typedef adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int, EdgeProperty>
        >
        Graph;

void get_embedding(Graph &g, std::vector<vec_t> &embedding) {
    // Initialize the interior edge index
    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);


    // Test for planarity - we know it is planar, we just want to
    // compute the planar embedding as a side-effect
//    std::vector<vec_t> embedding(num_vertices(g));
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )
            )
        std::cout << "Input graph is planar" << std::endl;
    else
        std::cout << "Input graph is not planar" << std::endl;
}

void test2(Graph &g, std::string name) {
    std::cout << name << "\n";
    std::vector<vec_t> embedding(num_vertices(g));
    get_embedding(g, embedding);

    std::cout << "Independent set\nresult: " << baker2<independent_set>(g) << "\n";
}

void unit2()
{

    Graph g6(5);
    add_edge(0,1,g6);
    add_edge(1,2,g6);
    add_edge(2,3,g6);
    add_edge(3, 0, g6);

    add_edge(4, 5, g6);
    add_edge(5, 6, g6);
    add_edge(6, 4, g6);

    add_edge(4, 0, g6);
    add_edge(4, 1, g6);
    add_edge(5, 2, g6);
//    add_edge(1, 6, g6);

    test2(g6, "");

//    Graph g7(5);
//    add_edge(0,1,g7);
//    add_edge(1,2,g7);
//    add_edge(2,3,g7);
//    add_edge(3, 0, g7);
//
//    add_edge(4, 0, g7);
//    add_edge(4, 1, g7);
//    add_edge(4, 2, g7);
//    add_edge(4, 3, g7);
//
//    add_edge(0, 5, g7);
//    add_edge(1,5, g7);
//
//    test2(g7, "");
}

