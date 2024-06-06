#include "BGAL/Model/ManifoldModel.h"
#include "BGAL/Model/Model_Iterator.h"

namespace BGAL {

_ManifoldModel::_MMEdge::_MMEdge() : _Segment3(), _id_left_vertex(-1), _id_right_vertex(-1), _id_opposite_vertex(-1), _id_left_edge(-1), _id_right_edge(-1), _id_reverse_edge(-1), _id_face(-1) {}

_ManifoldModel::_MMEdge::_MMEdge(const _Point3 &in_s, const _Point3 &in_t)
    : _Segment3(in_s, in_t), _id_left_vertex(-1), _id_right_vertex(-1), _id_opposite_vertex(-1), _id_left_edge(-1), _id_right_edge(-1), _id_reverse_edge(-1), _id_face(-1) {}

_ManifoldModel::_ManifoldModel() {}

_ManifoldModel::_ManifoldModel(const std::string &in_file_name) : _Model(in_file_name) { preprocess_model_(); }

_ManifoldModel::_ManifoldModel(const std::vector<_Point3> &in_vertices, const std::vector<_Model::_MFace> &in_faces) : _Model() {
  _vertices = in_vertices;
  _faces = in_faces;
  compute_normal_boundingbox_();
  preprocess_model_();
}

_ManifoldModel::_ManifoldModel(const _ManifoldModel &in_mmodel) {
  _vertices = in_mmodel._vertices;
  _faces = in_mmodel._faces;
  compute_normal_boundingbox_();
  preprocess_model_();
}

void _ManifoldModel::preprocess_model_() {
  creat_edges_from_vertices_faces_();
  arrange_neighs_of_vertex_face_();
}

_Edge_Iterator _ManifoldModel::edge_begin() const { return _Edge_Iterator(this, 0); }

_FE_Iterator _ManifoldModel::fe_begin(const int &fid) const {
  if (fid < 0 || fid >= number_faces_()) {
    throw std::runtime_error("Beyond the index!");
  }
  return _FE_Iterator(this, fid, 0);
}

_FF_Iterator _ManifoldModel::ff_begin(const int &fid) const {
  if (fid < 0 || fid >= number_faces_()) {
    throw std::runtime_error("Beyond the index!");
  }
  return _FF_Iterator(this, fid, 0);
}

_VV_Iterator _ManifoldModel::vv_begin(const int &vid) const {
  if (vid < 0 || vid >= number_vertices_()) {
    throw std::runtime_error("Beyond the index!");
  }
  return _VV_Iterator(this, vid, 0);
}

_VE_Iterator _ManifoldModel::ve_begin(const int &vid) const {
  if (vid < 0 || vid >= number_vertices_()) {
    throw std::runtime_error("Beyond the index!");
  }
  return _VE_Iterator(this, vid, 0);
}

_VF_Iterator _ManifoldModel::vf_begin(const int &vid) const {
  if (vid < 0 || vid >= number_vertices_()) {
    throw std::runtime_error("Beyond the index!");
  }
  return _VF_Iterator(this, vid, 0);
}

void _ManifoldModel::creat_edges_from_vertices_faces_() {
  _edges.reserve(800);
  std::map<std::pair<int, int>, int> vertids_to_edgeid;
  size_t id_face = -1;
  for (auto &face : _faces) {
    id_face++;
    std::vector<int> eid_in_f(3, -1);
    for (int j = 0; j < 3; ++j) {
      int post = (j + 1) % 3;
      int pre = (j + 2) % 3;

      int left_vertex = face[pre];
      int right_vertex = face[j];

      std::pair<int, int> vert_ids({left_vertex, right_vertex});

      auto it = vertids_to_edgeid.find(vert_ids);
      if (it != vertids_to_edgeid.end()) {
        int loc_in_edges = it->second;
        if (_edges[loc_in_edges]._id_opposite_vertex != -1) {
          // std::cout << "--Repeated edges!" << std::endl;
          // throw std::runtime_error("Repeated edges!");
          continue;
        }
        eid_in_f[j] = loc_in_edges;
        _edges[loc_in_edges]._id_opposite_vertex = face[post];
        _edges[loc_in_edges]._id_face = id_face;
      } else {
        _MMEdge e1(_vertices[left_vertex], _vertices[right_vertex]);
        e1._id_left_vertex = left_vertex;
        e1._id_right_vertex = right_vertex;
        e1._id_face = id_face;
        e1._id_opposite_vertex = face[post];
        e1._id_reverse_edge = _edges.size() + 1;
        _edges.push_back(e1);

        eid_in_f[j] = _edges.size() - 1;

        vertids_to_edgeid[vert_ids] = eid_in_f[j];

        _MMEdge e2(_vertices[right_vertex], _vertices[left_vertex]);
        e2._id_left_vertex = right_vertex;
        e2._id_right_vertex = left_vertex;
        e2._id_face = -1;
        e2._id_opposite_vertex = -1;
        e2._id_reverse_edge = _edges.size() - 1;
        _edges.push_back(e2);
        vertids_to_edgeid[{right_vertex, left_vertex}] = _edges.size() - 1;
      }
    }
    for (int j = 0; j < 3; ++j) {
      int index = eid_in_f[j];
      if (index < 0 || index >= _edges.size()) {
        continue;
      }
      _edges[eid_in_f[j]]._id_left_edge = _edges[eid_in_f[(j + 2) % 3]]._id_reverse_edge;
      _edges[eid_in_f[j]]._id_right_edge = _edges[eid_in_f[(j + 1) % 3]]._id_reverse_edge;
    }
  }
}

void _ManifoldModel::arrange_neighs_of_vertex_face_() {
  _neight_edge_of_vertices.clear();
  _neight_edge_of_vertices.resize(_vertices.size(), -1);
  _neigh_edge_of_faces.clear();
  _neigh_edge_of_faces.resize(_faces.size(), -1);
  _degree_of_vertices.clear();
  _degree_of_vertices.resize(_vertices.size(), 0);
  for (int i = 0; i < _edges.size(); ++i) {
    if (_edges[i]._id_face == -1)
      continue;
    if (_neight_edge_of_vertices[_edges[i]._id_left_vertex] == -1 || _edges[_edges[i]._id_reverse_edge]._id_face == -1) {
      _neight_edge_of_vertices[_edges[i]._id_left_vertex] = i;
    }
    ++_degree_of_vertices[_edges[i]._id_left_vertex];
    if (_faces[_edges[i]._id_face][0] == _edges[i]._id_opposite_vertex) {
      _neigh_edge_of_faces[_edges[i]._id_face] = i;
    }
  }
  _isolated_vertices.resize(0);
  for (int i = 0; i < _neight_edge_of_vertices.size(); ++i) {
    if (_neight_edge_of_vertices[i] == -1) {
      _isolated_vertices.push_back(i);
    } else {
      int degree = 1;
      int start_edge = _neight_edge_of_vertices[i];
      int cur_edge = _edges[start_edge]._id_left_edge;
      while (_edges[cur_edge]._id_face != -1 && cur_edge != start_edge) {
        degree++;
        cur_edge = _edges[cur_edge]._id_left_edge;
      }
      if (degree != _degree_of_vertices[i]) {
        throw std::runtime_error("complex vertex: " + std::to_string(i));
      }
    }
  }
}
} // namespace BGAL