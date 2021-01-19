#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "problems.hpp"
#include "baker-outer-planar.hpp"

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

void test(Graph &g) {
    std::vector<vec_t> embedding(num_vertices(g));
    get_embedding(g, embedding);
    std::cout << "Independent set\nresult: " << baker<independent_set>(g, embedding, 1) << "\n";
    std::cout << "Vertex cover\nresult: " << baker<vertex_cover>(g, embedding, 1) << "\n";
//    std::cout << "Dominating set\nresult: " << baker<dominating_set>(g, embedding, 1) << "\n";
//    std::cout << "Edge dominating set\nresult: " << baker<edge_dominating_set>(g, embedding, 1) << "\n";
}

int main(int argc, char** argv)
{

    Graph g1(9);
    add_edge(0,1,g1);
    add_edge(1,2,g1);
    add_edge(2,3,g1);
    add_edge(3,4,g1);
    add_edge(4,5,g1);
    add_edge(5, 6,g1);
    add_edge(6,7,g1);
    add_edge(7,8,g1);
    add_edge(8,0,g1);

    add_edge(0,2,g1);
    add_edge(2,6,g1);
    add_edge(3,5,g1);
    add_edge(6,8,g1);

    test(g1);

    Graph g2(5);
    add_edge(0,1,g2);
    add_edge(1,2,g2);
    add_edge(2,3,g2);
    add_edge(3,4,g2);
    add_edge(0,2,g2);
    add_edge(0,3,g2);
    add_edge(0,4,g2);

    test(g2);


    return 0;
}

