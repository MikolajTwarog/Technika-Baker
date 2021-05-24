//
// Created by mikolajtwarog on 2021-04-30.
//

#ifndef TECHNIKA_BAKER_FIND_OUTER_FACE_HPP
#define TECHNIKA_BAKER_FIND_OUTER_FACE_HPP

void find_outer_face(PlanarEmbedding& embedding, std::vector<int>& outer_face) {
    std::map<graph_traits<Graph>::edge_descriptor, std::vector<int> > faces;
    std::vector<std::vector<int> > vertices_in_face;
    face_getter<Edge> my_vis(&faces, vertices_in_face);
    level_face_traversal<Graph>(embedding, my_vis);

    for (const auto& face : vertices_in_face) {
        Edge current_e;
        Edge next_e(face[0], face[1], nullptr);

        bool res = true;
        for (int i = 1; i < face.size() - 1; i++) {
            current_e = next_e;
            next_e.m_source = face[i];
            next_e.m_target = face[i + 1];

            int dis = (get_edge_it(next_e, face[i], embedding) - get_edge_it(current_e, face[i], embedding)
                       + embedding[face[i]].size()) % embedding[face[i]].size();
            if (dis != 1 && dis != embedding[face[i]].size() - 1) {
                res = false;
                break;
            }

        }

        if (res) {
            for (int v : face) {
                outer_face.push_back(v);
            }
            return;
        }
    }
}

#endif //TECHNIKA_BAKER_FIND_OUTER_FACE_HPP
