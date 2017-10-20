#include "Model.h"
#include "util.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void parse_levels(vector<int>& levels, const vector<int>& fields_ind, const vector<string>& fields)
{
  levels.resize(fields_ind.size());
  for (unsigned int i = 0; i<fields_ind.size(); i++)
    levels[i] = atoi(fields[fields_ind[i]].c_str());
}

void init_field_indices(vector<int>& field_ind, const vector<string> &fields, const vector<string>& titles)
{
  field_ind.resize(fields.size(), -1);
  for (unsigned int i = 0; i<titles.size(); i++)
    for (unsigned int j = 0; j<fields.size(); j++)
      if (titles[i] == fields[j])
	field_ind[j] = i;
  for (unsigned int i = 0; i<field_ind.size(); i++)
    massert(field_ind[i] != -1, "field %s not found", fields[i].c_str());
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
//////////////////////////////////////////////////////////////////////////////////////////////////

void Model::read_fends_table(const string& fn,
			     const string& fend_field,
			     const string& contig_field,
			     const string& coord_field,
			     const string& anchor_field,
			     int& n_anchors)
{
  cerr << "reading fends table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  vector<int> m_features_ind(m_features.size, -1);

  // parse header
  split_line(in, fields, delim);
  init_field_indices(m_features_ind, m_features.fields, fields);

  int fend_ind = get_field_index(fend_field, fields);
  int contig_ind = get_field_index(contig_field, fields);
  int coord_ind = get_field_index(coord_field, fields);
  int anchor_ind = get_field_index(anchor_field, fields);

  n_anchors = 0;
  int n_fends = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig_id = fields[contig_ind];
    int anchor = atoi(fields[anchor_ind].c_str());

    // skip 0: no anchor
    if (anchor == 0)
      continue;

    n_fends++;
    Fend fend;
    fend.index = atoi(fields[fend_ind].c_str());
    fend.contig_index = contig2int(contig_id);
    fend.coord = atoi(fields[coord_ind].c_str());

    // classify fend according to its model levels
    parse_levels(fend.levels, m_features_ind, fields);

    // add fend
    m_fends_contig[fend.contig_index].push_back(fend);
    m_fends_anchor[anchor-1].push_back(fend);
    m_fend_contig_map[fend.index] = contig_id;

    // add contig once
    if (m_contigs_flag.find(contig_id) == m_contigs_flag.end()) {
      m_contigs_flag[contig_id] = true;
      m_contigs.push_back(Contig(contig_id, anchor-1));
    }

    n_anchors = max(n_anchors, anchor);
  }

  in.close();

  cerr << "number of fends: " << n_fends << endl;
}

void Model::read_feature_table(string fn, vector<double>& result, int mat_dim)
{
  cout << "reading feature table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';
  const unsigned int width = 3;
  bool first = true;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      return;
    if (first)
      {
	first = false;
	continue;
      }

    assert(fields.size() == width);

    int i = atoi(fields[0].c_str());
    int j = atoi(fields[1].c_str());
    double prob = atof(fields[2].c_str());

    massert(i <= mat_dim, "field1, value %d out of matrix range: %d", i, mat_dim);
    massert(j <= mat_dim, "field2, value %d out of matrix range: %d", j, mat_dim);

    int v1 = MATRIX_INDEX(i, j, mat_dim);
    int v2 = MATRIX_INDEX(j, i, mat_dim);

    result[v1] = prob;
    result[v2] = prob;
  }
  in.close();
}

void Model::read_feature_tables(const vector<string>& feature_fns,
				const vector<int>& sizes)
{
  m_model_matrices.resize(feature_fns.size());
  for (unsigned int i=0; i < feature_fns.size(); i++)
    {
      vector<double>& matrix = m_model_matrices[i];
      matrix.resize(sizes[i]*sizes[i]);
      read_feature_table(feature_fns[i], matrix, sizes[i]);
    }
}

double Model::compute_simple(vector<Fend>& fends_x, vector<Fend>& fends_y)
{
  double total = 0;
  vector<Fend>::const_iterator jt_x, jt_y;
  vector<int>& msizes = m_features.sizes;

  for(jt_x = fends_x.begin(); jt_x != fends_x.end(); jt_x++)
    for(jt_y = fends_y.begin(); jt_y != fends_y.end(); jt_y++) {
      const Fend& fend_x = *jt_x;
      const Fend& fend_y = *jt_y;

      // skip if self-interaction or if within contig to itself
      if (fend_x.index == fend_y.index)
	continue;

      double value = m_prior;
      for(unsigned int mi = 0; mi < m_model_matrices.size(); mi++)
	value *= m_model_matrices[mi][MATRIX_INDEX(fend_x.levels[mi], fend_y.levels[mi], msizes[mi])];
      total += value;
    }
  return total;
}

void Model::intersect_fends(vector<Fend>& fends_x, vector<Fend>& fends_y, vector<Fend>& fends_common)
{
  sort(fends_x.begin(), fends_x.end(), comp_fends());
  sort(fends_y.begin(), fends_y.end(), comp_fends());

  // X && Y
  set_intersection(fends_x.begin(), fends_x.end(),
		   fends_y.begin(), fends_y.end(), back_inserter(fends_common), comp_fends());
}


double Model::compute(vector<Fend>& fends_x, vector<Fend>& fends_y)
{
  vector<Fend> fends_common;
  intersect_fends(fends_x, fends_y, fends_common);

  double total_common = compute_simple(fends_common, fends_common);
  double total = compute_simple(fends_x, fends_y);
  total -= total_common;
  return total;
}

double Model::compute(string contig, int anchor)
{
  massert(m_contig_index_map.find(contig) != m_contig_index_map.end(), "index of contig %s not found", contig.c_str());
  int contig_index = m_contig_index_map[contig];
  massert(m_fends_contig.find(contig_index) != m_fends_contig.end(), "contig not found in m_fends_contig");

  massert(m_fends_anchor.find(anchor) != m_fends_anchor.end(), "anchor not found in m_fends_anchor");

  vector<Fend>& fends_x = m_fends_contig[contig_index];
  vector<Fend>& fends_y = m_fends_anchor[anchor];

  return compute(fends_x, fends_y);
}

double Model::compute(string contig1, string contig2)
{
  massert(m_contig_index_map.find(contig1) != m_contig_index_map.end(), "index of contig %s not found", contig1.c_str());
  massert(m_contig_index_map.find(contig2) != m_contig_index_map.end(), "index of contig %s not found", contig2.c_str());
  int contig_index1 = m_contig_index_map[contig1];
  int contig_index2 = m_contig_index_map[contig2];
  massert(m_fends_contig.find(contig_index1) != m_fends_contig.end(), "contig not found in m_fends_contig");
  massert(m_fends_contig.find(contig_index2) != m_fends_contig.end(), "contig not found in m_fends_contig");

  vector<Fend>& fends_x = m_fends_contig[contig_index1];
  vector<Fend>& fends_y = m_fends_contig[contig_index2];
  return compute(fends_x, fends_y);
}

Model::Model(ModelFeatures features, double prior)
  : m_features(features), m_prior(prior)
{
  read_feature_tables(m_features.fns, m_features.sizes);
}
