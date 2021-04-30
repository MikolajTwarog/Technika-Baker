//
// Created by mikolajtwarog on 2021-04-30.
//

#ifndef TECHNIKA_BAKER_VISITORS_HPP
#define TECHNIKA_BAKER_VISITORS_HPP

template <typename Edge>
struct face_getter : public boost::planar_face_traversal_visitor
{
    std::map<Edge, std::vector<int> > &faces;
    std::vector<std::vector<int> >& vertices_in_face;
    int current_face = 0;

    face_getter(std::map<Edge, std::vector<int> >& f, std::vector<std::vector<int> >& o)
            : faces(f), vertices_in_face(o){ }

    void begin_face() {
        vertices_in_face.emplace_back();
    }

    void end_face() {
//        std::cout << "\n";
        current_face++;
    }

    void next_edge(Edge e)
    {
//        std::cout << e << " " << current_face << "\n";
        faces[e].push_back(current_face);
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        vertices_in_face[current_face].push_back(v);
    }
};


template <typename Edge, typename Problem, typename PlanarEmbedding>
struct tree_builder : public boost::planar_face_traversal_visitor
{
    std::map<Edge, std::vector<int> > &faces;
    ::tree<Problem> &tree;
    int current_face = 0;
    int last_vertex;
    Graph graph;
    PlanarEmbedding embedding;

    tree_builder(std::map<Edge, std::vector<int> >& f, ::tree<Problem>& t, Graph& g, PlanarEmbedding& emb)
            : faces(f), tree(t), graph(g), embedding(emb){ }

    void end_face() {
        current_face++;
    }

    template <typename Vertex>
    void next_vertex(Vertex v) {
        tree[current_face].face.push_back(graph.local_to_global(v));
        last_vertex = graph.local_to_global(v);
    }

    void next_edge(Edge e)
    {
        if(current_face != tree.outer_face) {
            int global_source = e.m_source;
            int global_target = e.m_target;

            int neighbor = faces[e][0] == current_face ? faces[e][1] : faces[e][0];

//            std::cout << "\n" << e << "\n";
//            for (int i : faces[e]) {
//                std::cout << i << "\n";
//            }

            if (neighbor == tree.outer_face) {
                tree.emplace_back();
                int last = tree.size() - 1;
                tree[current_face].children.push_back(last);
                tree[last].parent = current_face;
                tree[last].label.first = last_vertex;
                tree[last].label.second = last_vertex == global_source ? global_target : global_source;
//                tree[last].my_tree = &tree;

//                if (e != *embedding[target(e, graph)].begin())
//                    std::swap(tree[last].label.first, tree[last].label.second);
            } else {
                tree[current_face].children.push_back(neighbor);
            }
        }
    }
};

#endif //TECHNIKA_BAKER_VISITORS_HPP
