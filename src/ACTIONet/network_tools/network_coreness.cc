#include <ACTIONet.h>
#include <set>

namespace ACTIONet {
    uvec compute_core_number(sp_mat &G) {        		
		unsigned int i, j = 0;
		unsigned int no_of_nodes = G.n_rows;   


		// Construct node neighborhood sets
		vector<vector <unsigned int>> N(no_of_nodes);				 
		sp_mat::const_iterator it     = G.begin();
		sp_mat::const_iterator it_end = G.end();
		
		uvec cores(no_of_nodes); cores.zeros();
		for(; it != it_end; ++it) {
			N[it.row()].push_back((unsigned int)it.col());
			cores[it.row()]++;
		}	 
		unsigned int maxdeg = max(cores);

		/* degree histogram */
		uvec bin(maxdeg+1); bin.zeros();	
		for (i = 0; i < no_of_nodes; i++) {
			bin[ (unsigned int)cores[i] ] ++;
		}

		/* start pointers */
		j = 0;
		for (i = 0; i <= maxdeg; i++) {
			unsigned int k = bin[i];
			bin[i] = j;
			j += k;
		}

		/* sort in vert (and corrupt bin) */
		uvec pos(no_of_nodes); pos.zeros();
		uvec vert(no_of_nodes); vert.zeros();
		for (i = 0; i < no_of_nodes; i++) {
			pos[i] = bin[(unsigned int)cores[i]];
			vert[pos[i]] = i;
			bin[(unsigned int)cores[i]] ++;
		}

		/* correct bin */
		for (i = maxdeg; i > 0; i--) {
			bin[i] = bin[i - 1];
		}
		bin[0] = 0;
  

		/* this is the main algorithm */
		for (i = 0; i < no_of_nodes; i++) {
			unsigned int v = vert[i];
			
			for (j = 0; j < N[v].size(); j++) {
				unsigned int u = (N[v])[j];
				
				if (cores[u] > cores[v]) {
					unsigned int du = (unsigned int) cores[u];
					unsigned int pu = pos[u];
					unsigned int pw = bin[du];
					unsigned int w = vert[pw];
					if (u != w) {
						pos[u] = pw;
						pos[w] = pu;
						vert[pu] = w;
						vert[pw] = u;
					}
					bin[du] ++;
					cores[u] --;
				}
			}
		}

		return(cores);
	}
	
    uvec compute_induced_core_number(sp_mat &G, uvec mask) {        		
		unsigned int i, j = 0;
		unsigned int no_of_nodes = G.n_rows;   


		// Construct node neighborhood sets
		vector<vector <unsigned int>> N(no_of_nodes);				 
		sp_mat::const_iterator it     = G.begin();
		sp_mat::const_iterator it_end = G.end();
		
		uvec cores(no_of_nodes); cores.zeros();
		for(; it != it_end; ++it) {
			if(mask[it.row()] && mask[it.col()]) {
				N[it.row()].push_back((unsigned int)it.col());
				cores[it.row()]++;
			}
		}	 
		unsigned int maxdeg = max(cores);

		/* degree histogram */
		uvec bin(maxdeg+1); bin.zeros();	
		for (i = 0; i < no_of_nodes; i++) {
			bin[ (unsigned int)cores[i] ] ++;
		}

		/* start pointers */
		j = 0;
		for (i = 0; i <= maxdeg; i++) {
			unsigned int k = bin[i];
			bin[i] = j;
			j += k;
		}

		/* sort in vert (and corrupt bin) */
		uvec pos(no_of_nodes); pos.zeros();
		uvec vert(no_of_nodes); vert.zeros();
		for (i = 0; i < no_of_nodes; i++) {
			pos[i] = bin[(unsigned int)cores[i]];
			vert[pos[i]] = i;
			bin[(unsigned int)cores[i]] ++;
		}

		/* correct bin */
		for (i = maxdeg; i > 0; i--) {
			bin[i] = bin[i - 1];
		}
		bin[0] = 0;
  

		/* this is the main algorithm */
		for (i = 0; i < no_of_nodes; i++) {
			unsigned int v = vert[i];
			
			for (j = 0; j < N[v].size(); j++) {
				unsigned int u = (N[v])[j];
				
				if (cores[u] > cores[v]) {
					unsigned int du = (unsigned int) cores[u];
					unsigned int pu = pos[u];
					unsigned int pw = bin[du];
					unsigned int w = vert[pw];
					if (u != w) {
						pos[u] = pw;
						pos[w] = pu;
						vert[pu] = w;
						vert[pw] = u;
					}
					bin[du] ++;
					cores[u] --;
				}
			}
		}

		return(cores);
	}
	
	vec compute_archetype_core_centrality(sp_mat &G, uvec sample_assignments) {
		
		vec conn = zeros(G.n_cols);
		for(int i = 0; i <= max(sample_assignments); i++) {
			uvec mask = conv_to<uvec>::from(sample_assignments == i);
			
			if(sum(mask) == 0) {
				continue;
			}
			uvec induced_coreness = compute_induced_core_number(G, mask);
			
			uvec idx = find(mask > 0);
			vec v = conv_to<vec>::from(induced_coreness(idx));
			double sigma = stddev(v);
			vec z = v/ (sigma == 0? 1:sigma);
			
			conn(idx) = z;
		}
		
		return(conn);
	}
}






