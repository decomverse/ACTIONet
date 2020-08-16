#include <ACTIONet.h>

#include <GraphHelper.h>
#include <RBConfigurationVertexPartition.h>
#include <CPMVertexPartition.h>
#include <Optimiser.h>

namespace ACTIONet {
	vec signed_cluster(sp_mat A, double resolution_parameter = 1.0, uvec initial_clusters = uvec(), int seed = 0) {
		int nV = A.n_rows;
		int nE = A.n_nonzero;
		
		vec clusters = zeros(nV);
						
		igraph_t g;
		igraph_vector_t v;
		igraph_vector_init(&v, 2*nE);
		vector<double> edge_weights(nE);


		int idx = 0;

		sp_mat::iterator it     = A.begin();
		sp_mat::iterator it_end = A.end();

		for(; it != it_end; ++it) {
			edge_weights[idx] = (*it);
			VECTOR(v)[2*idx] = it.row();
			VECTOR(v)[2*idx+1] = it.col();
			
			idx++;
		}
		igraph_create(&g, &v, nV, 0);

		Graph *G = new Graph(&g, edge_weights);
		//printf("idx = %d, V = %d, E = %d, total_weight = %f, weighted = %d, directed = %d\n", idx, G->vcount(), G->ecount(), G->total_weight(), G->is_weighted(), G->is_directed());
		
		CPMVertexPartition* partition;
		if(initial_clusters.n_elem == nV) {
			
			vector<size_t> membership(nV);
			for(int i = 0; i < nV; i++) {
				membership[i] = initial_clusters(i);
			}
			partition = new CPMVertexPartition(G, membership, resolution_parameter);		
		}
		else {
		 partition = new CPMVertexPartition(G, resolution_parameter);
		}

		Optimiser *opt = new Optimiser(seed);
		opt->optimise_partition(partition);
	
		
		for(int i = 0; i < nV; i++) {
			clusters(i) = partition->membership(i)+1;
		}
		
	    partition->destructor_delete_graph = true;
		delete(G);
		igraph_vector_destroy(&v);
		igraph_destroy(&g);
		delete(opt);
		
		return(clusters);
	}
	
	vec unsigned_cluster(sp_mat A, double resolution_parameter = 1.0, uvec initial_clusters = uvec(), int seed = 0) {
		int nV = A.n_rows;
		int nE = A.n_nonzero;
		
		vec clusters = zeros(nV);
						
		igraph_t g;
		igraph_vector_t v;
		igraph_vector_init(&v, 2*nE);
		vector<double> edge_weights(nE);


		int idx = 0;

		sp_mat::iterator it     = A.begin();
		sp_mat::iterator it_end = A.end();

		for(; it != it_end; ++it) {
			edge_weights[idx] = (*it);
			VECTOR(v)[2*idx] = it.row();
			VECTOR(v)[2*idx+1] = it.col();
			
			idx++;
		}
		igraph_create(&g, &v, nV, 0);

		Graph *G = new Graph(&g, edge_weights);
		//printf("idx = %d, V = %d, E = %d, total_weight = %f, weighted = %d, directed = %d\n", idx, G->vcount(), G->ecount(), G->total_weight(), G->is_weighted(), G->is_directed());
		

		RBConfigurationVertexPartition* partition;
		if(initial_clusters.n_elem == nV) {			
			vector<size_t> membership(nV);
			for(int i = 0; i < nV; i++) {
				membership[i] = initial_clusters(i);
			}
			partition = new RBConfigurationVertexPartition(G, membership, resolution_parameter);		
		}
		else {
		 partition = new RBConfigurationVertexPartition(G, resolution_parameter);
		}



		Optimiser *opt = new Optimiser(seed);
		opt->optimise_partition(partition);
	
		for(int i = 0; i < nV; i++) {
			clusters(i) = partition->membership(i)+1;
		}
		
	    partition->destructor_delete_graph = true;
		delete(G);
		igraph_vector_destroy(&v);
		igraph_destroy(&g);
		delete(opt);
		
		return(clusters);
	}
	

	mat unsigned_cluster_batch(sp_mat A, vec resolutions, uvec initial_clusters = uvec(), int seed = 0) {
		int nV = A.n_rows;
		int nE = A.n_nonzero;
								
		igraph_t g;
		igraph_vector_t v;
		igraph_vector_init(&v, 2*nE);
		vector<double> edge_weights(nE);


		int idx = 0;

		sp_mat::iterator it     = A.begin();
		sp_mat::iterator it_end = A.end();

		for(; it != it_end; ++it) {
			edge_weights[idx] = (*it);
			VECTOR(v)[2*idx] = it.row();
			VECTOR(v)[2*idx+1] = it.col();
			
			idx++;
		}
		igraph_create(&g, &v, nV, 0);

		Graph *G = new Graph(&g, edge_weights);
		//printf("idx = %d, V = %d, E = %d, total_weight = %f, weighted = %d, directed = %d\n", idx, G->vcount(), G->ecount(), G->total_weight(), G->is_weighted(), G->is_directed());
		

		RBConfigurationVertexPartition* partition;
		if(initial_clusters.n_elem == nV) {			
			vector<size_t> membership(nV);
			for(int i = 0; i < nV; i++) {
				membership[i] = initial_clusters(i);
			}
			partition = new RBConfigurationVertexPartition(G, membership);		
		}
		else {
		 partition = new RBConfigurationVertexPartition(G);
		}


		Optimiser *opt = new Optimiser(seed);
		mat clusters;
		for(int j = 0; j < resolutions.n_elem; j++) {
			partition->resolution_parameter = resolutions(j);
			opt->optimise_partition(partition);


			vec cur_clusters(nV);
			for(int i = 0; i < nV; i++) {
				cur_clusters(i) = partition->membership(i);
			}						
			
			
			int kk = (int)arma::max(cur_clusters);
			mat M = zeros(nV, kk+1);			
			for(int i = 0; i < M.n_cols; i++) {
				uvec idx = find(cur_clusters == i);				

				vec v = M.col(i);
				v(idx).ones();
				M.col(i) = v;
			}
			if(j == 0) {
				clusters = M;
			} else {
				clusters = join_horiz(clusters, M);
			}			
		}
		

		
	    partition->destructor_delete_graph = true;
		delete(G);
		igraph_vector_destroy(&v);
		igraph_destroy(&g);
		delete(opt);
		
		return(clusters);
	}


}
