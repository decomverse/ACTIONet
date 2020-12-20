#include <ACTIONet.h>

namespace ACTIONet {
double *l1, *l2, *w;
int *match1, *match2, *v1, *v2;
int *s, *tt, *deg, *offset, *list;

/**
 * n the number of nodes
 * m the number of nodes
 * nedges the number of edges
 * vv1 is the source for each of the nedges
 * vv2 is the target for each of the nedges
 * weight is the weight of each of the nedges
 * out1 is a vector of length at most min(n,m),
 * out2 is a vector of length at most min(n,m),
 * noutedges is the number of out edges
 */
double MWM_driver(int n, int m, int nedges, double *vv1, double *vv2,
                  double *weight, double *out1, double *out2, int *noutedges) {
  double ret, al;
  int i, j, k, p, q, r, t1, t2;

  for (i = 0; i < nedges; i++) {
    v1[i] = (int)(vv1[i] + .5);
    v2[i] = (int)(vv2[i] + .5);
  }
  for (i = 0; i < n; i++) {
    offset[i] = 0;
    deg[i] = 1;
  }
  for (i = 0; i < nedges; i++)
    deg[v1[i]]++;
  for (i = 1; i < n; i++)
    offset[i] = offset[i - 1] + deg[i - 1];
  for (i = 0; i < n; i++)
    deg[i] = 0;
  for (i = 0; i < nedges; i++) {
    list[offset[v1[i]] + deg[v1[i]]] = v2[i];
    w[offset[v1[i]] + deg[v1[i]]] = weight[i];
    deg[(int)v1[i]]++;
  }
  for (i = 0; i < n; i++) {
    list[offset[i] + deg[i]] = m + i;
    w[offset[i] + deg[i]] = 0;
    deg[i]++;
  }
  for (i = 0; i < n; i++) {
    l1[i] = 0;
    for (j = 0; j < deg[i]; j++) {
      if (w[offset[i] + j] > l1[i])
        l1[i] = w[offset[i] + j];
    }
  }
  for (i = 0; i < n; i++) {
    match1[i] = -1;
  }
  for (i = 0; i < n + m; i++) {
    l2[i] = 0;
    match2[i] = -1;
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n + m; j++)
      tt[j] = -1;
    s[p = q = 0] = i;
    for (; p <= q; p++) {
      if (match1[i] >= 0)
        break;
      k = s[p];
      for (r = 0; r < deg[k]; r++) {
        if (match1[i] >= 0)
          break;
        j = list[offset[k] + r];
        if (w[offset[k] + r] < l1[k] + l2[j] - 1e-8)
          continue;
        if (tt[j] < 0) {
          s[++q] = match2[j];
          tt[j] = k;
          if (match2[j] < 0) {
            for (; j >= 0;) {
              k = match2[j] = tt[j];
              p = match1[k];
              match1[k] = j;
              j = p;
            }
          }
        }
      }
    }
    if (match1[i] < 0) {
      al = 1e20;
      for (j = 0; j < p; j++) {
        t1 = s[j];
        for (k = 0; k < deg[t1]; k++) {
          t2 = list[offset[t1] + k];
          if (tt[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al) {
            al = l1[t1] + l2[t2] - w[offset[t1] + k];
          }
        }
      }
      for (j = 0; j < p; j++)
        l1[s[j]] -= al;
      for (j = 0; j < n + m; j++)
        if (tt[j] >= 0)
          l2[j] += al;
      i--;
      continue;
    }
  }
  ret = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < deg[i]; j++) {
      if (list[offset[i] + j] == match1[i]) {
        ret += w[offset[i] + j];
      }
    }
  }
  *noutedges = 0;
  for (i = 0; i < n; i++) {
    if (match1[i] < m)
      (*noutedges)++;
  }
  *noutedges = 0;
  for (i = 0; i < n; i++) {
    if (match1[i] < m) {
      out1[*noutedges] = i;
      out2[*noutedges] = match1[i];
      (*noutedges)++;
    }
  }

  return ret;
}

mat MWM_hungarian(mat &G) {
  int n = G.n_rows;
  int m = G.n_cols;
  int N = m + n;
  l1 = new double[N];
  l2 = new double[N];
  v1 = new int[m * n];
  v2 = new int[m * n];
  s = new int[N];
  tt = new int[N];
  match1 = new int[N];
  match2 = new int[N];
  offset = new int[N];
  deg = new int[N];
  list = new int[m * n + N];
  w = new double[m * n + N];

  mat G_matched = zeros(size(G));

  uvec idx = find(G);
  if (idx.n_elem == 0)
    return G_matched;

  int nedges = idx.n_elem;

  double *vv1 = new double[nedges];
  double *vv2 = new double[nedges];
  double *weight = new double[nedges];

  umat subs = ind2sub(size(G), idx);
  for (int i = 0; i < nedges; i++) {
    weight[i] = G(idx(i));
    vv1[i] = subs(0, i);
    vv2[i] = subs(1, i);
  }

  int match_no = std::min(m, n);
  double *ii = new double[match_no];
  double *jj = new double[match_no];

  int matched_edge_no;

  MWM_driver(n, m, nedges, vv1, vv2, weight, ii, jj, &matched_edge_no);

  for (int k = 0; k < matched_edge_no; k++) {
    G_matched(ii[k], jj[k]) = G(ii[k], jj[k]);
  }

  delete[] l1;
  delete[] l2;
  delete[] v1;
  delete[] v2;
  delete[] s;
  delete[] tt;
  delete[] match1;
  delete[] match2;
  delete[] offset;
  delete[] deg;
  delete[] list;
  delete[] w;

  return G_matched;
}

// From: Low Rank Spectral Network Alignment
// (https://dl.acm.org/citation.cfm?doid=3178876.3186128)
umat MWM_rank1(vec u, vec v, double u_threshold, double v_threshold) {
  int pair_no = std::min(u.n_elem, v.n_elem);

  // Rank-1 matching for each paired-archetypes
  vec u_sorted = sort(u, "descend");
  uvec u_sorted_idx = sort_index(u, "descend");

  vec v_sorted = sort(v, "descend");
  uvec v_sorted_idx = sort_index(v, "descend");

  /*
          // To randomly break ties
          u_sorted = u_sorted + stddev(u)*randn(size(u_sorted))/30;
          u_sorted.transform( [](double val) { return (val < 0? 0:val); } );

          v_sorted = v_sorted + stddev(v)*randn(size(v_sorted))/30;
          v_sorted.transform( [](double val) { return (val < 0? 0:val); } );
  */

  // To reduce false positves by aligning only high-quality pairs
  // vec uv_prod = sqrt(u_sorted(span(0, pair_no-1)) % v_sorted(span(0,
  // pair_no-1)));

  int top_rank = min(sum(u > u_threshold),
                     sum(v > v_threshold)); // sum(uv_prod > threshold);
  umat subs(2, top_rank);

  if (top_rank == 0)
    return subs;

  // Rank-based matching
  uvec rows = u_sorted_idx(span(0, top_rank - 1));
  uvec cols = v_sorted_idx(span(0, top_rank - 1));

  subs.row(0) = (rows).t();
  subs.row(1) = (cols).t();

  return subs;
}

} // namespace ACTIONet
