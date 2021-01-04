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

typedef property<edge_faces_t, std::vector<int> > EdgeProperty;

typedef adjacency_list
        <
                vecS,
                vecS,
                undirectedS,
                property<vertex_index_t, int>,
                property<edge_index_t, int, EdgeProperty>
        >
        Graph;

template <typename Edge>
struct my_visitor : public planar_face_traversal_visitor
{
    property_map<Graph, edge_faces_t>::type &faces;
    int current_face = 0;

    my_visitor(property_map<Graph, edge_faces_t>::type &f)
    : faces(f) { }

    void end_face() {
        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
        std::cout << e << " ";
        faces[e].push_back(current_face);
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

template <typename Edge>
struct tree_builder : public planar_face_traversal_visitor
{
    const property_map<Graph, edge_faces_t>::type &faces;
    std::vector<node> &tree;
    int current_face = 0;
    int parent_it = -1;
    int edge_it = 0;

    tree_builder(property_map<Graph, edge_faces_t>::type &f, std::vector<node> &t)
            : faces(f), tree(t){ }

    void end_face() {
        if(parent_it != -1)
            std::rotate(tree[current_face].children.begin(),
                        tree[current_face].children.begin()+parent_it,
                        tree[current_face].children.end());
        current_face++;
        parent_it = -1;
        edge_it = 0;
    }

    void next_edge(Edge e)
    {
        if(current_face > 0) {
            int neighbor = faces[e][0];
            if (neighbor == current_face)
                neighbor = faces[e][1];

            if (neighbor == 0) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
            } else if (tree[current_face].parent != neighbor && tree[neighbor].parent == -1) {
                tree[current_face].children.push_back(neighbor);
                tree[neighbor].parent = current_face;
            } else if ((tree[current_face].parent == neighbor || tree[neighbor].parent > -1) && parent_it == -1) {
                parent_it = edge_it;
            }

            edge_it++;
        }
    }
};


void merge(std::vector<int> &one, std::vector<int> &two) {
    std::vector<int> copy(one);

//    std::cout << "one\n";
//    for(auto j : one) {
//        std::cout << j << " ";
//    }
//    std::cout << "\ntwo\n";
//    for(auto j : two) {
//        std::cout << j << " ";
//    }
//    std::cout << "\n";

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
//        std::cout << "result\n";
        for(auto j : tree[root].val) {
//            std::cout << j << " ";
        }
//        std::cout << "\n";
    }
}

int baker(Graph &g, std::vector<std::vector< graph_traits<Graph>::edge_descriptor > > &embedding, int root) {
    property_map<Graph, edge_faces_t>::type faces = get(edge_faces_t(), g);
    my_visitor<graph_traits<Graph>::edge_descriptor> my_vis(faces);
    planar_face_traversal(g, &embedding[0], my_vis);
    std::vector<node> tree(my_vis.current_face);

    tree_builder<graph_traits<Graph>::edge_descriptor> tree_b(faces, tree);
    planar_face_traversal(g, &embedding[0], tree_b);

    calculate(tree, root);

//    int k=0;
//    for(auto i : tree) {
//        std::cout << k++ << " p " << i.parent << " c ";
//        for(auto j : i.children) {
//            std::cout << j << " ";
//        }
//        std::cout << "\n";
//    }
//
//    int y = 1;
//    calculate(tree, y);
//    std::cout << tree[y].val[0] << " " << tree[y].val[1] << " " << tree[y].val[2] << " " << tree[y].val[3] << "\n";

    return tree[root].val[0];
}

int main(int argc, char** argv)
{

    Graph g(9);
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
    property_map<Graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<Graph>::edges_size_type edge_count = 0;
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for(boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);


    // Test for planarity - we know it is planar, we just want to
    // compute the planar embedding as a side-effect
    typedef std::vector< graph_traits<Graph>::edge_descriptor > vec_t;
    std::vector<vec_t> embedding(num_vertices(g));
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::embedding =
                                             &embedding[0]
    )
            )
        std::cout << "Input graph is planar" << std::endl;
    else
        std::cout << "Input graph is not planar" << std::endl;

    std::cout << baker(g, embedding, 1) << "\n";

    return 0;
}