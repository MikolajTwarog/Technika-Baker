#include <iostream>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/biconnected_components.hpp>

#include <vector>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "../k-outerplanar/baker-k-outer-planar.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kouter
#include <boost/test/unit_test.hpp>
#include <boost/graph/subgraph.hpp>

using namespace boost;

typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

typedef std::vector< graph_traits<Graph>::edge_descriptor > vec_t;

//typedef adjacency_list
//        <
//                vecS,
//                vecS,
//                undirectedS,
//                property<vertex_index_t, int>,
//                property<edge_index_t, int, EdgeProperty>
//        >
//        Graph;

struct file_reader{
    std::ifstream file;

    file_reader(std::string fn):
        file("/home/mikolajtwarog/Desktop/licencjat/Technika-Baker/Boost_tests/test_graphs/" + fn){}

    bool next_graph(Graph& g) {
        if(!file.is_open()) {
            return false;
        }

        int n, m;

        if (!(file >> n >> m)) {
            return false;
        }

        int a, b;
        while (m--) {
            file >> a >> b;
            add_edge(a, b, g);
        }

        return true;
    }
};

int independent_set_(Graph& g) {
    const int n = num_vertices(g);
    int count = 1 << n;

    int mx = 0;

    for (int i = 0; i < count; i++) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        std::bitset<32> u(i);
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (u[e.m_source] && u[e.m_target]) {
                res = false;
                break;
            }
        }

        if (res) {
            int ones = 0;
            for (int j = 0; j < n; j++) {
                ones += u[j];
            }

            mx = std::max(mx, ones);
        }
    }

    return mx;
}

BOOST_AUTO_TEST_SUITE(kouter)
    BOOST_AUTO_TEST_CASE(smiple)
    {
        std::vector<independent_set> as;
        as.emplace_back(1, nullptr);

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

        BOOST_CHECK_EQUAL(baker2<independent_set>(g6), 3);
    }

    BOOST_AUTO_TEST_CASE(two_bicomps) {
        Graph g;
        add_edge(0, 1, g);
        add_edge(1, 2, g);
        add_edge(2, 0, g);
        add_edge(2, 3, g);
        add_edge(3, 4, g);
        add_edge(4, 2, g);

        add_edge(5, 6, g);
        add_edge(6, 7, g);
        add_edge(7, 8, g);
        add_edge(8, 5, g);

        add_edge(0, 5, g);
        add_edge(1, 6, g);
        add_edge(2, 7, g);

        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 4);
    }

    BOOST_AUTO_TEST_CASE(bridge) {
        Graph g;
        add_edge(0, 1, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 1);
    }

    BOOST_AUTO_TEST_CASE(point_componet) {
        Graph g;
        add_edge(0, 1, g);
        add_edge(0, 2, g);
        add_edge(0, 3, g);
        add_edge(1, 2, g);
        add_edge(1, 3, g);
        add_edge(2, 3, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 1);
    }

//    5 6
//    0 3  0 4  1 3  1 4  2 3  2 4
    BOOST_AUTO_TEST_CASE(five) {
        Graph g;
        add_edge(0, 3, g);
        add_edge(0, 4, g);
        add_edge(1, 3, g);
        add_edge(1, 4, g);
        add_edge(2, 3, g);
        add_edge(2, 4, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 3);
    }

//    6 10
//    0 1  0 2  0 4  1 3  1 5  2 4  2 5  3 4  3 5  4 5
    BOOST_AUTO_TEST_CASE(six) {
        Graph g;
        add_edge(0, 1, g);
        add_edge(0, 2, g);
        add_edge(0, 4, g);
        add_edge(1, 3, g);
        add_edge(1, 5, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        add_edge(4, 5, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 2);
    }

//    6 8
//    0 4  0 5  1 4  1 5  2 3  2 5  3 4  3 5
    BOOST_AUTO_TEST_CASE(six2) {
        Graph g;
        add_edge(0, 4, g);
        add_edge(0, 5, g);
        add_edge(1, 4, g);
        add_edge(1, 5, g);
        add_edge(2, 3, g);
        add_edge(2, 5, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 3);
    }

//    7 10
//    0 5  0 6  1 5  1 6  2 5  2 6  3 5  3 6  4 5  4 6
    BOOST_AUTO_TEST_CASE(seven) {
        Graph g;
        add_edge(0, 5, g);
        add_edge(0, 6, g);
        add_edge(1, 5, g);
        add_edge(1, 6, g);
        add_edge(2, 5, g);
        add_edge(2, 6, g);
        add_edge(3, 5, g);
        add_edge(3, 6, g);
        add_edge(4, 5, g);
        add_edge(4, 6, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 5);
    }

//    7 13
//    0 4  0 5  0 6  1 4  1 5  1 6  2 3  2 5  2 6  3 5  3 6  4 5  4 6
    BOOST_AUTO_TEST_CASE(seven2) {
        Graph g;
        add_edge(0, 4, g);
        add_edge(0, 5, g);
        add_edge(0, 6, g);
        add_edge(1, 4, g);
        add_edge(1, 5, g);
        add_edge(1, 6, g);
        add_edge(2, 3, g);
        add_edge(2, 5, g);
        add_edge(2, 6, g);
        add_edge(3, 5, g);
        add_edge(3, 6, g);
        add_edge(4, 5, g);
        add_edge(4, 6, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 5);
    }

//   7 11
//0 5  1 6  2 4  2 5  2 6  3 4  3 5  3 6  4 5  4 6  5 6
    BOOST_AUTO_TEST_CASE(seven3) {
        Graph g;
        add_edge(0, 5, g);
//        add_edge(0, 6, g);
//        add_edge(1, 4, g);
//        add_edge(1, 5, g);
        add_edge(1, 6, g);
//        add_edge(2, 3, g);
        add_edge(2, 4, g);
        add_edge(2, 5, g);
        add_edge(2, 6, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        add_edge(3, 6, g);
        add_edge(4, 5, g);
        add_edge(4, 6, g);
        add_edge(5, 6, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 4);
    }

//  8 11
//0 6  0 7  1 6  1 7  2 6  2 7  3 4  3 5  4 5  4 6  5 7
    BOOST_AUTO_TEST_CASE(eight) {
        Graph g;
        add_edge(0, 6, g);
        add_edge(0, 7, g);
        add_edge(1, 6, g);
        add_edge(1, 7, g);
//        add_edge(2, 5, g);
        add_edge(2, 6, g);
        add_edge(2, 7, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
//        add_edge(3, 6, g);
//        add_edge(3, 7, g);
        add_edge(4, 5, g);
        add_edge(4, 6, g);
//        add_edge(4, 7, g);
        add_edge(5, 7, g);
//        add_edge(6, 7, g);
        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 5);
    }

    BOOST_AUTO_TEST_CASE(four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            std::cout << 2*i + 1 << std::endl;
            i++;
            int result = baker2<independent_set>(g);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(five_vertices) {
        file_reader f("5vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            std::cout << 2*i + 1 << std::endl;
            i++;
            int result = baker2<independent_set>(g);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(six_vertices) {
        file_reader f("6vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            std::cout << 2*i + 1 << std::endl;
            i++;
            int result = baker2<independent_set>(g);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(seven_vertices) {
        file_reader f("7vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            std::cout << 2*i + 1 << std::endl;
            i++;
            int result = baker2<independent_set>(g);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(eight_vertices) {
        file_reader f("8vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            std::cout << 2*i + 1 << std::endl;
            i++;
            int result = baker2<independent_set>(g);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(art) {
        Graph g;
        add_edge(0, 5, g);
        add_edge(0, 6, g);
        add_edge(1, 2, g);
        add_edge(1, 5, g);
        add_edge(1, 6, g);
        add_edge(2, 5, g);
        add_edge(2, 6, g);
        add_edge(3, 4, g);
        add_edge(3, 5, g);
        add_edge(3, 6, g);
        add_edge(4, 6, g);
        std::vector<int> sub{0, 1, 5};
        Graph& g2 = g.create_subgraph(sub.begin(), sub.end());
        add_edge(5, 0, g2);
    }

BOOST_AUTO_TEST_SUITE_END()

