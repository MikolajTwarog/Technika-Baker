#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/ref.hpp>
#include <vector>

#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>


using namespace boost;

namespace boost {
    enum edge_faces_t { edge_faces };

    BOOST_INSTALL_PROPERTY(edge, faces);
}

template <typename Edge>
struct my_visitor : public planar_face_traversal_visitor
{
    std::map<Edge, std::vector<int> > &map;
    int current_face = 0;
    std::vector<std::vector<Edge> > &faces;

    my_visitor(std::map<Edge, std::vector<int> >  &m, std::vector<std::vector<Edge> > &f)
    : map(m), faces(f) { }

    void end_face() {
        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
        std::cout << e << " ";
        map[e].push_back(current_face);
        faces[current_face].push_back(e);
    }
};

struct node
{
    int parent;
    std::vector<int> children;
    std::vector<int> val;
    node(): parent(-1){
        val.push_back(0);
        val.push_back(1);
        val.push_back(1);
        val.push_back(-1);
    }
};

void merge(std::vector<int> &one, std::vector<int> &two) {
    std::vector<int> copy(one);

    std::cout << "one\n";
    for(auto j : one) {
        std::cout << j << " ";
    }
    std::cout << "\ntwo\n";
    for(auto j : two) {
        std::cout << j << " ";
    }
    std::cout << "\n";

    for(int i=0; i<4; i++){
        for(int j=0; j<2; j++) {
            if(one[i] == -1 || two[((i&1)*2) + j] == -1)
                continue;
            one[i] = std::max(one[i], copy[(i&2) + j] + two[((j<<1) + (i&1))] - j);
//            std::cout << copy[i] << " " << two[((i&1)*2) + j] << " " << (i&1) << " " << (copy[i] + two[((i&1)*2) + j] - (i&1)) << "\n";
        }
//        std::cout << one[i] << "\n";
    }
}

void calculate(std::vector<node> &tree, int root){
    if(tree[root].children.empty()) {
        return;
    }

    tree[root].val = tree[tree[root].children[0]].val;

    for(int i=1; i < tree[root].children.size(); i++) {
        calculate(tree, tree[root].children[i]);
        merge(tree[root].val, tree[tree[root].children[i]].val);
        std::cout << "result\n";
        for(auto j : tree[root].val) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
}

int main(int argc, char** argv)
{

    typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

    typedef adjacency_list
            <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int, EdgeProperty>
            >
            graph;

    graph g(9);
    add_edge(0,1,g);
    add_edge(1,2,g);
    add_edge(2,3,g);
    add_edge(3,4,g);
    add_edge(4,5,g);
    add_edge(5, 6,g);
    add_edge(6,7,g);
    add_edge(7,8,g);
    add_edge(8,0,g);

    add_edge(0,2,g);
    add_edge(2,6,g);
    add_edge(3,5,g);
    add_edge(6,8,g);


    // Initialize the interior edge index
    property_map<graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<graph>::edges_size_type edge_count = 0;
    graph_traits<graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);


    // Test for planarity - we know it is planar, we just want to
    // compute the planar embedding as a side-effect
    typedef std::vector< graph_traits<graph>::edge_descriptor > vec_t;
    std::vector<vec_t> embedding(num_vertices(g));
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )
            )
        std::cout << "Input graph is planar" << std::endl;
    else
        std::cout << "Input graph is not planar" << std::endl;


    std::map<graph_traits<graph>::edge_descriptor, std::vector<int> > m;
    std::vector<std::vector<graph_traits<graph>::edge_descriptor> > faces(10);

    my_visitor<graph_traits<graph>::edge_descriptor> my_vis(m, faces);
    planar_face_traversal(g, &embedding[0], my_vis);

    for(auto elem : m)
    {
        std::cout << elem.first << " " << elem.second[0] << " " << elem.second[1] << "\n";
    }

    std::vector<node> tree(my_vis.current_face);

    for(int face_num = 1; face_num < faces.size(); face_num++){

        int it = 0;
        int parent_it = -1;
        for(auto edge : faces[face_num]) {
            int neighbor;
            if (m[edge][0] == face_num)
                neighbor = m[edge][1];
            else
                neighbor = m[edge][0];

            if(neighbor == 0) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[face_num].children.push_back(last);
                tree[last].parent = face_num;
            } else if(tree[face_num].parent != neighbor && tree[neighbor].parent == -1) {
                tree[face_num].children.push_back(neighbor);
                tree[neighbor].parent = face_num;
            } else if((tree[face_num].parent == neighbor || tree[neighbor].parent > -1) && parent_it == -1) {
                parent_it = it;
            }
            it++;
        }

        if(parent_it != -1)
            std::rotate(tree[face_num].children.begin(), tree[face_num].children.begin()+parent_it, tree[face_num].children.end());

    }

    int k=0;
    for(auto i : tree) {
        std::cout << k++ << " p " << i.parent << " c ";
        for(auto j : i.children) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }

    int y = 1;
    calculate(tree, y);
    std::cout << tree[y].val[0] << " " << tree[y].val[1] << " " << tree[y].val[2] << " " << tree[y].val[3] << "\n";

    return 0;
}

//    graph g(9);
//    add_edge(0,1,g);
//    add_edge(1,2,g);
//    add_edge(2,3,g);
//    add_edge(3,4,g);
//    add_edge(4,5,g);
//    add_edge(5, 6,g);
//    add_edge(6,7,g);
//    add_edge(7,8,g);
//    add_edge(8,0,g);
//
//    add_edge(0,2,g);
//    add_edge(2,6,g);
//    add_edge(3,5,g);
//    add_edge(6,8,g);