#include "bmrf_functions.h"

#include "bmrf_constants.h"
#include "Timer.h"

/**
 * Take ppiArray and produce a subset based on the list of gids
 */
void getppisubnet(Uint2dArray &ppi, Uint1dArray &gid, Uint2dArray &sppi) {

	// Build a list of matching IDs (produces two list of indecies)
	Uint1dArray index1(MAX_SET_SIZE);
	Uint1dArray index2(MAX_SET_SIZE);
	findi(ppi, 0, gid, index1);
	findi(ppi, 1, gid, index2);
//	printf("index1: ");
//	index1.print();
//	printf("index2: ");
//	index2.print();

// Keep intersection of 1 and 2
	Uint1dArray index3(index1.x + index2.x);
	intersectSets(index1, index2, index3);
//	sortSet(index3);
//	printf("index3: ");
//	index3.print(size_index3);

// Extract those indicies and their connections from PPI data
	for (int k = 0; k < index3.x; ++k) {
		sppi(k, 0) = ppi(index3[k], 0);
		sppi(k, 1) = ppi(index3[k], 1);
	}
	sppi.shrink(index3.x);

}

/**
 * ttest, some fixes, and icdf, normalize (will reduce a NxM matrix to 1xM)
 */
void genescore(Float2dArray &x, Float2dArray &y, Double1dArray &zscore) {

	// Equal variances
	unsigned int size = x.x;
	float nx = x.y;
	float ny = y.y;
	float dfe = x.y + y.y - 2;
	Float1dArray mean_x(size);
	Float1dArray mean_y(size);
	mean(x, mean_x);
	mean(y, mean_y);
	Float1dArray difference(size);
	subtract(mean_x, mean_y, difference);
	Float1dArray s2x(size);
	Float1dArray s2y(size);
	var(x, mean_x, s2x);
	var(y, mean_y, s2y);

	Float1dArray test(size);
	double temp = 0.0f;
	for (int i = 0; i < size; ++i) {
		temp = sqrt((double(nx - 1) * s2x[i] + double(ny - 1) * s2y[i]) / dfe)
				* sqrt(1.0f / double(nx) + 1.0f / double(ny));
		test[i] = -1.0f * fabs(difference[i] / temp);
	}

	// Compute the correct p-value for the test
	Double1dArray p(size);
	tcdf(test, dfe, p);

	// Iterate through p and fix upper and lower values and invert
	float machine_epsilon = 1e-16;
	double ptemp;
	for (int i = 0; i < size; ++i) {
		ptemp = 1.0 - (2.0 * p[i]);
		if (ptemp >= 1.0f)
			ptemp = 1.0 - machine_epsilon;
		if (ptemp < machine_epsilon)
			ptemp = machine_epsilon;
		p[i] = ptemp;
	}

	// inverse CDF
	icdf_normal(p, zscore, 0.0f, 1.0f);

	// Scale
	double score_mean = mean(zscore);
	double score_std = stddev(zscore, score_mean);
	for (int i = 0; i < size; ++i) {
		zscore[i] = (zscore[i] - score_mean) / score_std;
	}

}

/**
 * Get scores from genes based on unique sppi, gid should be unique
 */
void extractScores(Uint2dArray &sppi, Uint1dArray &gid, Double1dArray &zscore,
		Uint1dArray &ugeneid, Double1dArray &subscore) {

	getUniqueGeneIds(sppi, ugeneid);

	// Find indices of the unique genes in gid (which is in a parallel structure
	// to zscore) to then extract zscores.
	unsigned int size_scores = 0;
	for (int uid = 0; uid < ugeneid.x; ++uid) {
		for (int id = 0; id < gid.x; ++id) {
			if (ugeneid[uid] == gid[id]) {
				subscore[size_scores] = zscore[id];
				size_scores += 1;
			}
		}
	}

	subscore.shrink(size_scores);
}

/**
 *
 */
void getUniqueGeneIds(Uint2dArray &sppi, Uint1dArray &ugeneid) {
	// Create unique set of gene IDs (this is passed back)
	Uint1dArray cgeneid(sppi.length());
	convert2dTo1d(sppi, cgeneid);
	sortSet(cgeneid);
	uniqueSet(cgeneid, ugeneid);
}

/**
 * find connections
 */
void g_conn(Uint2dArray &ppi, Uint1dArray &sgeneid, unsigned int index,
		Uint1dArray &igconn) {

	// Build a list of matching IDs (produces two list of indecies)
	/*
	Uint1dArray gid(1);
	gid[0] = sgeneid[index];

	Uint1dArray value1(MAX_SET_SIZE);
	Uint1dArray value2(MAX_SET_SIZE);
	findv(ppi, 0, gid, value1);
	findv(ppi, 1, gid, value2);

	// Keep intersection of 1 and 2
	Uint1dArray ids(value1.x + value2.x);
	unionSets(value1, value2, ids); // calls unique
	*/

	// Fast version
	Uint1dArray gid(1);
	gid[0] = sgeneid[index];
	Uint1dArray ids(MAX_SET_SIZE);
	finddoublev(ppi, gid, ids); // ids is non unique, nor sorted
	sortSet(ids);
	uniqueSetInplaceMustBeSorted(ids);

	// Extract the indicies from sgeneid
	intersectSetsReturnIndexAlreadyUnique(sgeneid, ids, igconn);
}

/**
 * x = pinv( (coef1*a1) + (coef2*a2) ). a1, a2 and x must be of same size.
 */
void pinvArrayAddWithCoefMul(const double &coef1, Double2dArray &a1,
		const double &coef2, Double2dArray &a2, Double2dArray &x) {

	// Matrix add with coef
	Double2dArray temp(x.x, x.y);
	for (int i = 0; i < a1.x; ++i) {
		for (int j = 0; j < a1.y; ++j) {
			temp(i, j) = coef1 * a1(i, j) + coef2 * a2(i, j);
		}
	}

	//pinv
	pinv(temp, x);
}

/**
 *
 */
void multiplyVectorMatrix(Double1dArray &a, Double2dArray &b,
		Double1dArray &c) {
	memset(c.data, 0, c.size());
	for (int col = 0; col < b.y; ++col) {
		for (int inner = 0; inner < a.x; ++inner) {
			c[col] += a[inner] * b(inner, col);
		}
	}
}

/**
 *
 */
void multiplyMatrixVector(Double2dArray &a, Double1dArray &b,
		Double1dArray &c) {
	memset(c.data, 0, c.size());
	for (int row = 0; row < a.x; ++row) {
		for (int inner = 0; inner < a.y; ++inner) {
			c[row] += a(row, inner) * b[inner];
		}
	}
}

/**
 *
 */
void scoreUpdate(Double2dArray &x2d, float gamma, Double1dArray &score,
		Double1dArray& x) {
	// Create new score
	Double1dArray temp(score.x);
	for (int i = 0; i < score.x; ++i) {
		temp[i] = (1.0f + (score[i] * gamma)) / float(x2d.x);
	}

	// matrix [ng,ng] x vector [ng,1] = vector [ng,1]
	multiplyMatrixVector(x2d, temp, x);
}

/**
 * (a-b)'*(a-b)
 */
float squareTransposeVectorSubtract(Double1dArray &f, Double1dArray &score) {

	// Subtract
	Double1dArray temp(f.x);
	subtract(f, score, temp);

	// transpose and vector multiply and subtract
	double ans = 0.0f;
	for (int row = 0; row < f.x; ++row) {
		ans += (temp[row] * temp[row]);
	}
	return ans;
}

/**
 * a'*b*a
 */
float squareTransposeMatrixVectorMultiply(Double1dArray &a, Double2dArray &b) {

	// b*a
	Double1dArray temp(a.x);
	multiplyVectorMatrix(a, b, temp);

	// transpose and vector multiply
	double ans = 0.0f;
	for (int row = 0; row < temp.x; ++row) {
		ans += (temp[row] * a[row]);
	}
	return ans;

}

/**
 * Compute the subnetwork score
 * Inputs:
 * gid: gene id vector
 * cand_network_id: gene ids in the interested network
 * zscore: z-score of gene ids in geneid
 * ppi: protein-portein interaction network
 * @return network potentials, negative of network score
 */
double mrfnetscore(Uint1dArray &geneid, Uint1dArray &cand_network_id,
		Double1dArray &zscore, Uint2dArray &ppi) {
	float gamma = 1.0f;
	float lambda = 1.0f;

	// Extract subnetwork
	Uint2dArray sppi(MAX_SET_SIZE, ppi.y);
	getppisubnet(ppi, cand_network_id, sppi);
//	printf("sppi = "); sppi.print();

	// Get scores from genes based on unique sppi
	Uint1dArray uniquegeneid(geneid.x);
	Double1dArray scoretemp(zscore.x);
	extractScores(sppi, geneid, zscore, uniquegeneid, scoretemp);
//	printf("scoretemp = "); scoretemp.print(scoretemp.x);

	// _shrug_
	unsigned int ng = uniquegeneid.x;
	unsigned int ne = std::max(sppi.x, sppi.y);
	Double1dArray a(ng);
	Double2dArray di(ng, ng);
	memset(di.data, 0, di.length() * sizeof(double));
	Uint1dArray igconn(sppi.x);
	for (int i = 0; i < ng; ++i) {
		unsigned int size_a = 0;
		g_conn(sppi, uniquegeneid, i, igconn);
		for (int conn = 0; conn < igconn.x; ++conn) {
			di(i, igconn[conn]) = -1.0f;
			di(igconn[conn], i) = -1.0f;
		}
		di(i, i) = igconn.x;
	}

	Double2dArray di2(ng, ng);
	memcpy(di2.data, di.data, di.size());
	for (int i = 0; i < ng; ++i) {
		for (int j = 0; j < ng; ++j) {
			di2(i, j) = di2(i, j) / sqrt(di(i, i));
		}
	}
	for (int i = 0; i < ng; ++i) {
		for (int j = 0; j < ng; ++j) {
			di2(j, i) = di2(j, i) / sqrt(di(i, i));
		}
	}
//	printf("di2 = "); di2.print();

	// Create eye matrix
	Double2dArray x_eye(ng, ng);
	memset(x_eye.data, 0, x_eye.length() * sizeof(double));
	for (int i = 0; i < ng; ++i) {
		x_eye(i, i) = 1.0f;
	}

	// Pseudo inverse
	Double2dArray x2d(ng, ng);
	pinvArrayAddWithCoefMul(2.0f * lambda / double(ne), di2, gamma / double(ng),
			x_eye, x2d);

	// Not sure what this is
	Double1dArray f(scoretemp.x);
	scoreUpdate(x2d, gamma, scoretemp, f);
//	printf("f = "); f.print(f.x);

	double part1 = -1.0f * mean(f);
	double part2 = lambda / double(ne)
			* squareTransposeMatrixVectorMultiply(f, di2);
	double part3 = gamma / double(2 * ng)
			* squareTransposeVectorSubtract(f, scoretemp);

//	printf("part1=%f + part2=%f + part3=%f \n", part1, part2, part3);
	return part1 + part2 + part3;
}

/**
 *
 */
void populateClasses(Float2dArray &labels, Uint1dArray &index_class1,
		Uint1dArray &index_class2, Float1dArray &classes) {
	unsigned int size_index1 = 0;
	unsigned int size_index2 = 0;
	for (int i = 0; i < labels.x; ++i) {
		if (labels(i, 0) == classes[0]) {
			index_class1[size_index1] = i;
			size_index1 += 1;
		}
		if (labels(i, 0) == classes[1]) {
			index_class2[size_index2] = i;
			size_index2 += 1;
		}
	}
	index_class1.shrink(size_index1);
	index_class2.shrink(size_index2);
}

/**
 * Random resampling
 */
/**
 * get portion of data based on classes (no random resampling)
 * data [ 367 x 40 ], class_labels [ 20 ] or [ 40 ], sampled_data [ 367 x <class_label_size> ];
 * data.x and sampled_data.x must match
 * class_labels.x <= data.y and class_labels.x == sampled_data.y
 */
void randsample(Float2dArray &data, Uint1dArray &class_labels,
		Float2dArray &randsampled_data) {
	Uint1dArray random_index_array(randsampled_data.y);
	unsigned int random_index = 0;
	for (int i = 0; i < randsampled_data.y; ++i) {
		random_index = class_labels[rrand(class_labels.x)];
		random_index_array[i] = random_index;

		// copy
		for (int j = 0; j < randsampled_data.x; ++j)
			randsampled_data(j, i) = data(j, random_index);
	}
//	random_index_array.print();
}

/**
 * get portion of data based on classes (no random resampling)
 * data [ 367 x 40 ], class_labels [ 20 ] or [ 40 ], sampled_data [ 367 x 20 ];
 * data.x and sampled_data.x must match
 * class_labels.x <= data.y and class_labels.x == sampled_data.y
 */
void classsample(Float2dArray &data, Uint1dArray &class_labels,
		Float2dArray &sampled_data) {
	for (int i = 0; i < sampled_data.y; ++i) {
		// copy
		for (int j = 0; j < sampled_data.x; ++j)
			sampled_data(j, i) = data(j, class_labels[i]);
	}
}

/**
 * Generate a candidate network
 */
void netcand(Uint2dArray &ppi, Uint1dArray &geneid,
		Uint1dArray &prev_network_ids, Uint1dArray &prev_network_distances,
		const unsigned int &upper_distance,
		Uint1dArray &candidate_network_id_Array,
		Uint1dArray &candidate_network_distance_Array) {

	Uint1dArray temp_net_id_array(MAX_SET_SIZE);
	temp_net_id_array.shrink(0);
	Uint1dArray temp_net_distance_array(MAX_SET_SIZE);
	temp_net_distance_array.shrink(0);

	// Generate the candidate nodes in the network with one step further
	bool sign = brand();
	if ((sign) && (prev_network_ids.x < 20)) {
		Uint1dArray value1(MAX_SET_SIZE);
		Uint1dArray value2(MAX_SET_SIZE);
		Uint1dArray ids(MAX_SET_SIZE);
		Uint1dArray cand_id_temp(MAX_SET_SIZE);
		Uint1dArray seed_id(1);
		Uint1dArray cand_id_temp2(MAX_SET_SIZE);
		Uint1dArray cand_distance_temp(MAX_SET_SIZE);
		Uint1dArray cand_id_expand(MAX_SET_SIZE);
		Uint1dArray cand_distance_expand(MAX_SET_SIZE);

		// Add node
		for (int i = 0; i < prev_network_ids.x; ++i) {

			// Build a list of matching IDs (produces two list of indecies)
			seed_id[0] = prev_network_ids[i];
			finddoublev(ppi, seed_id, ids);
			sortSet(ids);
			uniqueSetInplaceMustBeSorted(ids);

			// Extract the indicies from sgeneid
			intersectSetsAlreadyUnique(ids, geneid, cand_id_temp);

			// Since we don't consider topology change at this time, all previously
			// included nodes are removed
			diffSets(cand_id_temp, prev_network_ids, cand_id_temp2);

			for (int k = 0; k < cand_id_temp2.x; ++k)
				cand_distance_temp[k] = prev_network_distances[i] + 1;

			unionSets(temp_net_id_array, cand_id_temp2, cand_id_expand);

			long idx = -1;
			for (int j = 0; j < cand_id_expand.x; ++j) {
				unsigned int dist_temp1 = INT_MAX;
				unsigned int dist_temp2 = INT_MAX;
				idx = findi(temp_net_id_array, cand_id_expand[j]);
				if (idx >= 0)
					dist_temp1 = temp_net_distance_array[idx];

				idx = findi(cand_id_temp2, cand_id_expand[j]);
				if (idx >= 0)
					dist_temp2 = cand_distance_temp[idx];

				cand_distance_expand[j] = std::min(dist_temp1, dist_temp2);
			}

			unsigned int count_max = 0;
			for (int i = 0; i < cand_id_expand.x; ++i) {
				temp_net_id_array[count_max] = cand_id_expand[i];
				temp_net_distance_array[count_max] = cand_distance_expand[i];
				count_max++;
			}
			temp_net_id_array.shrink(count_max);
			temp_net_distance_array.shrink(count_max);
		}

		// Prune all the nodes whose distance is larger than upper_distance
		unsigned int count = 0;
		for (int i = 0; i < temp_net_distance_array.x; ++i) {
			if (temp_net_distance_array[i] <= upper_distance) {
				candidate_network_id_Array[count] = temp_net_id_array[i];
				candidate_network_distance_Array[count] =
						temp_net_distance_array[i];
				count++;
			}
		}
		candidate_network_id_Array.shrink(count);
		candidate_network_distance_Array.shrink(count);

	} else {

		// Remove Node
		// get the outer bound of genes into candidate gene list
		unsigned int count_max = 0;
		int max_distance = max(prev_network_distances);
		if (max_distance > 0) {
			for (int i = 0; i < prev_network_distances.x; ++i) {
				if (prev_network_distances[i] == max_distance) {
					candidate_network_id_Array[count_max] = prev_network_ids[i];
					candidate_network_distance_Array[count_max] =
							prev_network_distances[i];
					count_max++;
				}
			}
		}
		candidate_network_id_Array.shrink(count_max);
		candidate_network_distance_Array.shrink(count_max);
	}

//	printf("candidate_network_id_Array(%i) = ", candidate_network_id_Array.length()); candidate_network_id_Array.print();
//	printf("candidate_network_distance_Array(%i) = ", candidate_network_distance_Array.length()); candidate_network_distance_Array.print();

}

/**
 *
 */
void mrfsearchnet(Uint2dArray &ppi, Uint1dArray &geneids, Double1dArray &zscore,
		const unsigned int &seedgene, const unsigned int &distance,
		float temperature, Uint1dArray &sub_network, double &sub_netscore) {

	// Begin network with seed gene
	Uint1dArray prev_network_id_Array(sub_network.x);
	prev_network_id_Array[0] = seedgene;
	prev_network_id_Array.shrink(1);
	Uint1dArray prev_network_distance_Array(sub_network.x);
	prev_network_distance_Array[0] = 0;
	prev_network_distance_Array.shrink(1);

	Uint1dArray ind_gene(geneids.x);
	intersectSetsReturnIndexAlreadyUnique(geneids, prev_network_id_Array, ind_gene);

	// inital network score, the z-score of seed node
	Float1dArray netscore(MAX_ITERATIONS);
	netscore.shrink(0);
	netscore.append(-1.0f * zscore[ind_gene[0]]);

	Uint1dArray candidate_network_id_Array(sub_network.x);
	Uint1dArray candidate_network_distance_Array(sub_network.x);
	Uint1dArray sub_network_id_Array(sub_network.x);
	Uint1dArray sub_network_distance_Array(sub_network.x);
	candidate_network_id_Array.shrink(0);
	candidate_network_distance_Array.shrink(0);
	sub_network_id_Array.shrink(0);
	sub_network_distance_Array.shrink(0);
	copyArray(prev_network_id_Array, sub_network_id_Array);
	copyArray(prev_network_distance_Array, sub_network_distance_Array);

	// Iterate
	unsigned int iter = 1;
	unsigned int last_netchange = 0;
	while (iter < MAX_ITERATIONS) {
//		prettyPrintArray("sub_network_id_Array", sub_network_id_Array);
//		fflush(stdout);

		// Generate the candidate nodes in the network with one step further
		netcand(ppi, geneids, prev_network_id_Array,
				prev_network_distance_Array, distance,
				candidate_network_id_Array, candidate_network_distance_Array);

//		prettyPrintArray("candidate_network_id_Array", candidate_network_id_Array);

		// Check for an empty canidate set
		if (isempty(candidate_network_id_Array) == true) {
			iter++;
			temperature = 0.9f * temperature;
			continue;
		}

		// Get indicies matching seed genes, again if empty goto next iteration
		Uint1dArray candidate_indices(candidate_network_id_Array.x);
		intersectSetsReturnIndexAlreadyUnique(candidate_network_id_Array, geneids,
				candidate_indices);

		// assume prior probability of selecting genes are equal
		// randomly select one for candidate genes
		unsigned int idx_gene = candidate_indices[rrand(candidate_indices.x)];

		// given the candidate gene, calculate potential energy of new sub-network with and without new gene
		Uint1dArray new_id(1);
		new_id[0] = candidate_network_id_Array[idx_gene];
		Uint1dArray cand_network_id(sub_network.x);
		cand_network_id.shrink(0);
		intersectSets(prev_network_id_Array, new_id,
				cand_network_id);

		bool add_node = false;
		if (isempty(cand_network_id)) {
			// Add
			add_node = true;
			copyArray(prev_network_id_Array, cand_network_id);
			cand_network_id.append(new_id[0]);
		} else {
			// Remove
			diffSets(prev_network_id_Array, new_id, cand_network_id);
		}
		double netscore_cand = netscore[0];
		if (cand_network_id.x > 1)
			netscore_cand = mrfnetscore(geneids, cand_network_id, zscore, ppi);

		// simulated annealing to ADD/DELETE a new node
		// if new_energy < old_energy then ADD new node
		// else ADD new node with probability of exp(-(Difference of energy)/T);
		double netscore_diff = netscore_cand - double(netscore[netscore.x - 1]);
		if (netscore_diff < 0.0f) {
			// add/delete new node into current sub-network
			if (add_node) {
				// tack on end
				if (!sub_network_id_Array.append(new_id[0]))
					printf("failed to append\n");
				if (!sub_network_distance_Array.append(
						candidate_network_distance_Array[idx_gene]))
					printf("failed to append\n");
			} else {
				// remove node
				unsigned int resized_size = 0;
				for (int b = 0; b < sub_network_id_Array.x; ++b) {
					if (sub_network_id_Array[b] != new_id[0]) {
						sub_network_id_Array[resized_size] =
								sub_network_id_Array[b];
						sub_network_distance_Array[resized_size] =
								sub_network_distance_Array[b];
						resized_size++;
					}
				}
				sub_network_id_Array.shrink(resized_size);
				sub_network_distance_Array.shrink(resized_size);
			}
			netscore.append(float(netscore_cand));
			last_netchange = iter;
		} else {
			// ADD/DELETE new node with probability of exp(-(Difference of energy)/T);
			double sprob = exp(-1.0f * netscore_diff / temperature);
			bool fix = brand(sprob);
//			printf("netscore_diff %e, sprob %e, temp %e, fix %i\n",
//					netscore_diff, sprob, temperature, fix);
			if (fix == true) {
				if (add_node) {
					// tack on end
					if (!sub_network_id_Array.append(new_id[0]))
						printf("failed to append\n");
					if (!sub_network_distance_Array.append(
							candidate_network_distance_Array[idx_gene]))
						printf("failed to append\n");
				} else {
					// Remove node
					unsigned int resized_size = 0;
					for (int b = 0; b < sub_network_id_Array.x; ++b) {
						if (sub_network_id_Array[b] != new_id[0]) {
							sub_network_id_Array[resized_size] =
									sub_network_id_Array[b];
							sub_network_distance_Array[resized_size] =
									sub_network_distance_Array[b];
							resized_size++;
						}
					}
					sub_network_id_Array.shrink(resized_size);
					sub_network_distance_Array.shrink(resized_size);
				}
				netscore.append(float(netscore_cand));
				last_netchange = iter;
			}
		}

		// Get ready for next iteration
		iter++;
		temperature = 0.9 * temperature;
//		printf("iteration %i: network @ %f\n", iter, netscore[netscore.x-1]);
//	    prettyPrintArray("sub_network_id_Array", sub_network_id_Array);

		// Speed optimization
//		if (last_netchange + 400 < iter) {
//			iter = MAX_ITERATIONS;
//		}

		// Copy sub network for next iteration
		// This might look like a memory leak but we allocated them to the same size
		copyArray(sub_network_id_Array, prev_network_id_Array);
		copyArray(sub_network_distance_Array, prev_network_distance_Array);
	}

	// Return results
	sub_netscore = netscore[netscore.x - 1];
	copyArray(prev_network_id_Array, sub_network);
//	printf("subnetwork @ %f = ", sub_netscore); sub_network.print(sub_network.x);
//	printf("terminate mrfsearchnet\n"); fflush(stdout);


}

void printProgress(unsigned int count) {

	if (count % 10 || count == 0)
		printf(".");
	else
		printf(" ");

	fflush(stdout);
}

void printProgress2(unsigned int count, unsigned int seed) {

	if (seed == 0)
		printProgress(count);

}
