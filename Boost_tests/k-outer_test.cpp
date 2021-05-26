#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unistd.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/biconnected_components.hpp>

#include <vector>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "../Baker's Technique/bakers_technique.hpp"
#include "../tree_decomposition/bodlaender.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kouter
#include <boost/test/unit_test.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/dynamic_bitset.hpp>

#include "../utils/random_graph.hpp"

using namespace boost;

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

        std::cout << n << " " << m << std::endl;

        int a, b;
        while (m--) {
            file >> a >> b;
            add_edge(a, b, g);
            std::cout << a << " " << b << "  ";
        }

        std::cout << std::endl;

        return true;
    }
};

void make_graph(Graph& g, int m, std::string edges) {
    std::stringstream stream(edges);
    int a, b;
    while (m--) {
        stream >> a >> b;
        add_edge(a, b, g);
    }
}

void get_embedding(Graph& g, PlanarEmbedding& embedding, std::vector<int>& outer_face) {
    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);


    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )) {}
//        std::cout << "Input graph is planar" << std::endl;
    else {
        std::cout << "Input graph is not planar" << std::endl;
        return;
    }

    std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    face_getter<Edge> my_vis(&faces, vertices_in_face);
    level_face_traversal<Graph>(embedding, my_vis);

    for (int v : vertices_in_face[0]) {
        outer_face.push_back(v);
    }
}

int independent_set_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mx = 0;

    int last = n;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (s[e.m_source] && s[e.m_target]) {
                res = false;
                break;
            }
        }

        if (res) {
            int ones = s.count();
            mx = std::max(mx, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }

//        int f = s.find_first();
//        if (f < last) {
//            last = f;
//            std::cout << f << std::endl;
//        }
    }

    return mx;
}

int vertex_cover_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mn = INT16_MAX;

    int last = n;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        bool res = true;
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (!s[e.m_source] && !s[e.m_target]) {
                res = false;
                break;
            }
        }

        if (res) {
            int ones = s.count();
            mn = std::min(mn, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }

//        int f = s.find_first();
//        if (f < last) {
//            last = f;
//            std::cout << f << std::endl;
//        }
    }

    return mn;
}

int dominating_set_(Graph& g) {
    const int n = num_vertices(g);
    boost::dynamic_bitset<> s(n, 0);

    int mn = INT16_MAX;

    while (s.count() != n) {
        graph_traits<Graph>::edge_iterator ei, ei_end;
        boost::dynamic_bitset<> dom(n, 0);
        for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (s[e.m_source] || s[e.m_target]) {
                dom[e.m_source] = true;
                dom[e.m_target] = true;
            }
        }

        if (dom.count() == n) {
            int ones = s.count();
            mn = std::min(mn, ones);
        }

        for(int i = s.size() - 1; i >= 0; --i) {
            if ((s[i] ^= 0x1) == 0x1) {
                break;
            }
        }
    }

    return mn;
}

BOOST_AUTO_TEST_SUITE(kouter)
    BOOST_AUTO_TEST_CASE(przyklad) {
        Graph g;
        make_graph(g, 22, "0 1  0 2  0 4  0 7  0 8  1 2  1 8  2 3  2 8  3 4  3 5  4 5  5 6  5 7  6 7  6 12  7 11  8 9  9 10  11 12  11 13  12 13");
//        BOOST_CHECK_EQUAL(baker2<independent_set>(g), 5);
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<independent_set>(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

//    BOOST_AUTO_TEST_CASE(random) {
//        Graph g = random_graph(30, 120);
//        int result = baker2<independent_set>(g);
//        int expected = independent_set_(g);
//        BOOST_CHECK_EQUAL(result, expected);
//    }
//
//    BOOST_AUTO_TEST_CASE(technique) {
//        Graph g = random_graph(500, 700);
//        std::cout << num_edges(g) << std::endl;
//        int result = bakers_technique(g, 5);
//        int expected = baker2<independent_set>(g);
//        BOOST_CHECK_EQUAL(result, expected);
//    }

    BOOST_AUTO_TEST_CASE(vc_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_five_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_six_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_seven_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_eight_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(vc_nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<vertex_cover>(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

//    BOOST_AUTO_TEST_CASE(vc_random) {
//        Graph g = random_graph(30, 120);
//        int result = baker2<vertex_cover>(g);
//        int expected = vertex_cover_(g);
//        BOOST_CHECK_EQUAL(result, expected);
//    }

//    BOOST_AUTO_TEST_CASE(vc_technique) {
//        Graph g = random_graph(500, 700);
//        std::cout << num_edges(g) << std::endl;
//        int result = bakers_technique(g, 5);
//        int expected = baker2<independent_set>(g);
//        BOOST_CHECK_EQUAL(result, expected);
//    }

    BOOST_AUTO_TEST_CASE(ds_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

//    BOOST_AUTO_TEST_CASE(ds_four) {
//        Graph g;
//        make_graph(g, 7, "0 4  1 2  1 3  1 4  2 3  2 4  3 4");
//        BOOST_CHECK_EQUAL(baker2<dominating_set>(g), 1);
//    }

    BOOST_AUTO_TEST_CASE(ds_five_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_six_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_seven_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_eight_vertices) {
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
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(ds_nine_vertices) {
        file_reader f("9vertices");
//        std::string s = get_current_dir_name();

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = baker2<dominating_set>(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_simple)
    {
        Graph g;
        make_graph(g, 12, "0 1  0 2  0 3  0 4  1 2  1 3  1 5  2 4  2 5  3 4  3 5  4 5");

        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);

        bodlaender_vertex_cover(g, embedding, outer_face);
    }

    BOOST_AUTO_TEST_CASE(tr_four_vertices) {
        file_reader f("4vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_five_vertices) {
        file_reader f("5vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_six_vertices) {
        file_reader f("6vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_seven_vertices) {
        file_reader f("7vertices");

        int i = 0;
        bool res = true;
        while (true) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_vertex_cover(g, embedding, outer_face);
            int expected = vertex_cover_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }


BOOST_AUTO_TEST_SUITE_END()

