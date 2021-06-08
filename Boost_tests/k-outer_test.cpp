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
#include <chrono>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include "../Baker's Technique/bakers_technique.hpp"
#include "../tree_decomposition/bodlaender_impl.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE kouter
#include <boost/test/unit_test.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/dynamic_bitset.hpp>

#include "../utils/random_graph.hpp"

using namespace boost;

struct file_reader{
    std::ifstream file;
    std::string edges;

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

        edges = std::to_string(n) + " " + std::to_string(m) + "\n";

        int a, b;
        while (m--) {
            file >> a >> b;
            add_edge(a, b, g);
            edges += std::to_string(a) + " " + std::to_string(b) + "  ";
        }

        edges += "\n";

        return true;
    }

    bool next_graph2(Graph& g) {
        edges.clear();
        std::string s;
        std::set <std::pair<int, int> > es;
        while (getline(file, s)) {
            if (s.empty()) {
                break;
            }
            else {
                std::istringstream tmp(s);
                int a, b;
                char nic;
                tmp >> a >> nic >> b;
                a--; b--;
                std::pair<int, int> e(std::min(a,b), std::max(a,b));
                if (a != b && es.find(e) == es.end()) {
                    add_edge(a, b, g);
                    es.insert(e);
                    edges += std::to_string(a) + " " + std::to_string(b) + "  ";
                }
            }
        }
        edges = std::to_string(num_vertices(g)) + " " + std::to_string(num_edges(g)) + "\n" + edges;
        return true;
    }

    std::string get_last_edges(){
        return edges;
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
    BOOST_AUTO_TEST_CASE(technique) {
        Graph g = random_graph(100, 200);
//        Graph g;
//        make_graph(g, 14, "0 7  1 6  2 5  2 6  2 8  3 4  3 7  3 8  4 7  4 8  5 6  5 7  5 8  6 8");
        std::cout << num_edges(g) << std::endl;
        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);
        int result = bakers_technique(g, embedding, outer_face, 3, Baker, vc);
//        int expected = baker2<independent_set>(g, embedding, outer_face);
        BOOST_CHECK_EQUAL(result, 1);
    }

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

    BOOST_AUTO_TEST_CASE(ds_four) {
        Graph g;
        make_graph(g, 14, "0 7  1 6  2 5  2 6  2 8  3 4  3 7  3 8  4 7  4 8  5 6  5 7  5 8  6 8");
        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);
        int result = baker2<dominating_set>(g, embedding, outer_face);
        BOOST_CHECK_EQUAL(result, dominating_set_(g));
    }

//    0 5  0 6  1 5  1 6  2 5  2 6  3 5  3 6  4 5  4 6  5 6
    BOOST_AUTO_TEST_CASE(ds_four_2) {
        Graph g;
        make_graph(g, 13, "0 8  1 5  2 6  2 7  3 7  3 8  4 6  4 7  4 8  5 7  5 8  6 8  7 8");
        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);
        int result = baker2<dominating_set>(g, embedding, outer_face);
        BOOST_CHECK_EQUAL(result, dominating_set_(g));
    }

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
            if (result != expected) {
                break;
            }
        }
    }

    BOOST_AUTO_TEST_CASE(tr_simple)
    {
        Graph g;
        make_graph(g, 11, "0 1  1 2  2 3  3 4  4 5  5 6  6 0  1 3  1 4  1 5  1 6");

        PlanarEmbedding embedding(num_vertices(g));
        std::vector<int> outer_face;
        get_embedding(g, embedding, outer_face);

        bodlaender_vertex_cover(g, embedding, outer_face);
    }

    BOOST_AUTO_TEST_CASE(tr_vc_seven_vertices) {
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

    BOOST_AUTO_TEST_CASE(tr_vc_eight_vertices) {
        file_reader f("8vertices");

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

    BOOST_AUTO_TEST_CASE(tr_vc_nine_vertices) {
        file_reader f("9vertices");

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

    BOOST_AUTO_TEST_CASE(tr_is_four_vertices) {
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
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_seven_vertices) {
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
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_eight_vertices) {
        file_reader f("8vertices");

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
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_is_nine_vertices) {
        file_reader f("9vertices");

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
            int result = bodlaender_independent_set(g, embedding, outer_face);
            int expected = independent_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_lcc_seven_vertices) {
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
            int result = bodlaender_dominating_set_lcc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_lcc_eight_vertices) {
        file_reader f("8vertices");

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
            int result = bodlaender_dominating_set_lcc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_ecc_seven_vertices) {
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
            int result = bodlaender_dominating_set_ecc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_ecc_eight_vertices) {
        file_reader f("8vertices");

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
            int result = bodlaender_dominating_set_ecc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(tr_ds_ecc_nine_vertices) {
        file_reader f("9vertices");

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
            int result = bodlaender_dominating_set_ecc(g, embedding, outer_face);
            int expected = dominating_set_(g);
            BOOST_CHECK_EQUAL(result, expected);
        }
    }

    BOOST_AUTO_TEST_CASE(cos) {
        file_reader f("performance_test_graphs/7-outer");

        int i = 10;
        bool res = true;
        while (i--) {
            Graph g;
            res = f.next_graph(g);
            if (!res) {
                break;
            }
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            int result = bodlaender_independent_set(g, embedding, outer_face);
        }
    }

    BOOST_AUTO_TEST_CASE(graph_sorter) {
        file_reader f("performance_test_graphs/graphs");

        int z = 10000;
        bool res = true;
        std::vector<std::ofstream> files(15);
        for (int i = 1; i < 15; i++) {
            files[i].open("/home/mikolajtwarog/Desktop/licencjat/Technika-Baker/Boost_tests/test_graphs/performance_test_graphs/" + std::to_string(i) + "-outer",
                          std::ofstream::out | std::ofstream::app);
            if (!files[i].is_open()) {
                return;
            }
        }
        std::vector<int> graphs(10);
        while (z--) {
            Graph g;
            res = f.next_graph2(g);
            if (num_vertices(g) == 0) {
                continue;
            }
//            std::cout << f.get_last_edges() << std::endl;
            PlanarEmbedding embedding(num_vertices(g));
            std::vector<int> outer_face;
            get_embedding(g, embedding, outer_face);
            std::vector<int> vertex_level(num_vertices(g));
            std::vector< std::vector<Edge> > aaaa;
            int level = name_levels(embedding, outer_face, vertex_level, aaaa);
            files[level] << f.get_last_edges() << "\n";
        }
    }

    BOOST_AUTO_TEST_CASE(graph_sorter2) {
        file_reader f("performance_test_graphs/graphs");

        int z = 2000;
        bool res = true;
        std::ofstream file;
        file.open("/home/mikolajtwarog/Desktop/licencjat/Technika-Baker/Boost_tests/test_graphs/performance_test_graphs/small_graphs",
                      std::ofstream::out | std::ofstream::app);
        std::vector< std::vector<std::string> > graphs(4000);
        while (z--) {
            Graph g;
            res = f.next_graph2(g);
            if (num_vertices(g) == 0) {
                continue;
            }
            graphs[num_vertices(g)].push_back(f.get_last_edges());
        }
        for (auto& n : graphs) {
            for (auto& g : n) {
                file << g << "\n";
            }
        }
    }

    BOOST_AUTO_TEST_CASE(k_outer_baker_performance) {
        bool res = true;
        std::vector< std::pair<double, int> > results(15);
        for (int i = 1; i < 15; i++) {
            file_reader f("performance_test_graphs/" + std::to_string(i) + "-outer");
            int z=10;
            while (z--) {
                Graph g;
                res = f.next_graph(g);
                if (!res) {
                    break;
                }
                PlanarEmbedding embedding(num_vertices(g));
                std::vector<int> outer_face;
                get_embedding(g, embedding, outer_face);
                auto start = std::chrono::steady_clock::now();
                baker2<independent_set>(g, embedding, outer_face);
                auto stop = std::chrono::steady_clock::now();
                results[i].first += std::chrono::duration<double, std::milli>(stop - start).count();
                results[i].second++;
            }
        }
        std::ofstream file("/home/mikolajtwarog/Desktop/licencjat/Technika-Baker/Boost_tests/results/k-outer_baker_performance");
        for (int i = 1; i < 15; i++) {
            file << i << " " << results[i].first / results[i].second << "\n";
        }
    }


BOOST_AUTO_TEST_SUITE_END()

