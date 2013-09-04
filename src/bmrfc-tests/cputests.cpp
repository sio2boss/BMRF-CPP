#include "cputests.h"

/**
 *
 */
void test_matvec() {

	printf("\n===test_matvec====\n\n");

	// Matrix Multiply
	unsigned int n = 3;
	unsigned int m = 4;

	printf("matrixMultipy:\n");
	Float2dArray m1(n, m);
	Float2dArray m2(m, n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			m1(i, j) = frand();
		}
	}
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			m2(i, j) = frand();
		}
	}
	Float2dArray r(n, n);
	printf("  m1 = ");
	m1.print();
	printf("  m2 = ");
	m2.print();

	multiplyMatrix(m1, m2, r);
	printf("  r = ");
	r.print();
	printf("\n");

	// Allocate the components
	n = 5;
	m = 5;

	printf("pinvArrayAddWithCoefMul:\n");
	Double2dArray a(n, m);
	Double2dArray b(n, m);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			a(i, j) = frand();
			b(i, j) = frand();
		}
	}
	Double2dArray result(n, m);
	double coef1 = frand();
	double coef2 = frand();

	// Debug
	printf("  a1 = ");
	a.print();
	printf("  a2 = ");
	b.print();
	printf("  coef1 = %f\n", coef1);
	printf("  coef2 = %f\n", coef2);

	// do it
	pinvArrayAddWithCoefMul(coef1, a, coef2, b, result);
	printf("  x = ");
	result.print();

	// Another case
	Double2dArray s1(2, 2);
	Double2dArray s2(2, 2);
	s2(0, 0) = 1.0f;
	s2(0, 1) = 0.0f;
	s2(1, 0) = 0.0f;
	s2(1, 1) = 1.0f;
	s1(0, 0) = 1.0f;
	s1(0, 1) = -1.0f;
	s1(1, 0) = -1.0f;
	s1(1, 1) = 1.0f;
	Double2dArray x(2, 2);

	printf("  s1 = ");
	s1.print();
	printf("  s2 = ");
	s2.print();
	printf("  coef1 = 1\n");
	printf("  coef2 = 0.5\n");

	pinvArrayAddWithCoefMul(1, s1, 0.5, s2, x);
	printf("  x = ");
	x.print();

}

/**
 * Test unique
 */
void test_unique() {
	printf("\n===test_unique====\n\n");
	Uint1dArray a(15);
	Uint1dArray c(a.max);

// 1
	a[0] = 3;
	a.shrink(1);
	copyArray(a, c);
	uniqueSetInplaceMustBeSorted(c);
	printf("just one element:\n");
	printf("  A was: ");
	a.print();
	printf("  A is: ");
	c.print();

	// 2
	a[0] = 3;
	a[1] = 3;
	a.shrink(2);
	copyArray(a, c);
	uniqueSetInplaceMustBeSorted(c);
	printf("just 2 element, duplicates:\n");
	printf("  A was: ");
	a.print();
	printf("  A is: ");
	c.print();

	// 2
	a[0] = 3;
	a[1] = 4;
	a.shrink(2);
	copyArray(a, c);
	uniqueSetInplaceMustBeSorted(c);
	printf("just 2 element, no duplicates:\n");
	printf("  A was: ");
	a.print();
	printf("  A is: ");
	c.print();

	// 2
	a[0] = 3;
	a[1] = 4;
	a[2] = 5;
	a.shrink(3);
	copyArray(a, c);
	uniqueSetInplaceMustBeSorted(c);
	printf("just 3 element, no duplicates:\n");
	printf("  A was: ");
	a.print();
	printf("  A is: ");
	c.print();

	// many
	for (int i = 0; i < a.max; ++i) {
		a[i] = rand() % a.max;
	}
	a.shrink(a.max);
	sortSet(a);
	copyArray(a, c);
	uniqueSetInplaceMustBeSorted(c);
	printf("man elements, duplicates:\n");
	printf("  A was: ");
	a.print(a.x);
	printf("  A is: ");
	c.print(c.x);

}

/**
 * Test union, must handle no duplicates and duplicates
 */
void test_union() {
	printf("\n===test_union====\n\n");
	Uint1dArray a(5);
	Uint1dArray b(3);
// resulting set has a maximum of length(a)+length(b)
	Uint1dArray c(5 + 3);
	unsigned int size_c = 0;

// duplicates
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 9;
	}
	unionSets(a, b, c);
	printf("duplicates:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A union B is: ");
	c.print();

// no overlap
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 30;
	}
	unionSets(a, b, c);
	printf("no duplicates:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A union B is: ");
	c.print();
}

/**
 * Test intersect, three cases: full, partial and no overlap
 */
void test_intersect() {
	printf("\n===test_intersect====\n\n");
	Uint1dArray a(5);
	Uint1dArray b(3);
	Uint1dArray c(3);
	unsigned int size_c = 0;

// full overlap
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 10;
	}
	intersectSets(a, b, c);
	printf("full overlap:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A intersect B is: ");
	c.print(c.x);
	intersectSetsReturnIndex(a, b, c);
	printf("  index: ");
	c.print(c.x);
	intersectSetsReturnIndexAlreadyUnique(a, b, c);
	printf("  indexAlreadyUnique: ");
	c.print(c.x);

// partial
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 12;
	}
	intersectSets(a, b, c);
	printf("partial overlap:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A intersect B is: ");
	c.print(c.x);
	intersectSetsReturnIndex(a, b, c);
	printf("  index: ");
	c.print(c.x);
	intersectSetsReturnIndexAlreadyUnique(a, b, c);
	printf("  indexAlreadyUnique: ");
	c.print(c.x);


// no overlap
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 30;
	}
	intersectSets(a, b, c);
	printf("no overlap:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A intersect B is: ");
	c.print(c.x);
	intersectSetsReturnIndex(a, b, c);
	printf("  index: ");
	c.print(c.x);
	intersectSetsReturnIndexAlreadyUnique(a, b, c);
	printf("  indexAlreadyUnique: ");
	c.print(c.x);

}

/**
 * Test union, must handle no duplicates and duplicates
 */
void test_setdiff() {
	printf("\n===test_setdiff====\n\n");
	Uint1dArray a(5);
	Uint1dArray b(3);
	Uint1dArray c(5);
	unsigned int size_c = 0;

// duplicates
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 9;
	}
	diffSets(a, b, c);
	printf("duplicates:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A setdiff B is: ");
	c.print(c.x);

// no overlap
	for (int i = 0; i < a.x; ++i) {
		a[i] = i + 8;
	}
	for (int i = 0; i < b.x; ++i) {
		b[i] = i + 30;
	}
	diffSets(a, b, c);
	printf("no duplicates:\n");
	printf("  A is: ");
	a.print();
	printf("  B is: ");
	b.print();
	printf("  A setdiff B is: ");
	c.print(c.x);
}

/**
 *
 */
void test_getppisubnet() {

	printf("\n===test_getppisubnet====\n\n");

// Load ppi
	Uint2dArray ppi(8, 2);
	for (int i = 0; i < ppi.x; ++i) {
		ppi(i, 0) = i;
		ppi(i, 1) = i + 1;
	}
	printf("ppi network:\n");
	ppi.print();

// Load gene ids
	Uint1dArray gid(4);
	gid[0] = 6;
	gid[1] = 7;
	gid[2] = 2;
	gid[3] = 3;
	printf("gid list: ");
	gid.print();
	printf("\n");

// ppi subnet
	Uint2dArray sppi(8, 2);
	getppisubnet(ppi, gid, sppi);
	printf("extracted subnetwork:\n");
	sppi.print();

}

/**
 * SVD and pinv
 */
void test_matops() {

	Double2dArray a(5, 5);

	for (int i = 0; i < a.x; ++i) {
		for (int j = 0; j < a.y; ++j) {
			a(i, j) = (float) round(rand() / 10e6);
		}
	}

	// PINV
	printf("a=");
	a.print();
	Double2dArray a_pinv(5, 5);
	pinv(a, a_pinv);
	printf("a_pinv=");
	a_pinv.print();

}

/**
 *
 */
void test_basicstats() {

	printf("\n===test_basicstats====\n\n");

	printf("One-dim:");
	Float1dArray oned(5);
	oned[0] = 4.0f;
	oned[1] = -2.0f;
	oned[2] = 1.0f;
	oned[3] = 9.0f;
	oned[4] = 5.0f;
	printf("  using array: ");
	oned.print();

// mean
	float oned_mean = mean(oned);
	printf("  mean = %f, should be 3.4\n", oned_mean);

// var
	float oned_var = var(oned, oned_mean);
	printf("  var = %f, should be 17.3000\n", oned_var);

// std
	float oned_std = stddev(oned, oned_mean);
	printf("  std = %f, should be 4.1593\n", oned_std);
	printf("\n");
	printf("\n");

	printf("Two-dim:");
	Float2dArray twod(2, 3);
	twod(0, 0) = 4.0f;
	twod(0, 1) = -2.0f;
	twod(0, 2) = 1.0f;
	twod(1, 0) = 9.0f;
	twod(1, 1) = 5.0f;
	twod(1, 2) = 7.0f;
	printf("  using array: ");
	twod.print();
	printf("\n");

// mean
	Float1dArray twod_mean(2);
	mean(twod, twod_mean);
	printf("  mean =    ");
	twod_mean.print();
	printf("  should be [ 1.0000 7.0000 ]\n\n");

// var
	Float1dArray temp(2);
	var(twod, twod_mean, temp);
	printf("  var =     ");
	temp.print();
	printf("  should be [ 9.0000 4.0000 ]\n\n");

// std
	stddev(twod, twod_mean, temp);
	printf("  std =     ");
	temp.print();
	printf("  should be [ 3.0000  2.0000 ]\n\n");

}

/**
 *
 */
void test_cdfs() {

	printf("\n===test_cdfs====\n\n");
	Float1dArray in1(1);
	Double1dArray out1(1);
	in1[0] = -0.4205;
	float dfe = 9;

	// test tcdf
	printf("tcdf:\n");
	tcdf(in1, dfe, out1);
	printf("  input:     ");
	in1.print();
	printf("  dfe:       %f\n", dfe);
	printf("  output:    ");
	out1.print();
	printf("  should be  [ 0.3420 ]\n\n");

	// test icdf
	Double1dArray in2(1);
	Double1dArray out2(1);
	printf("icdf:\n");
	in2[0] = 0.3420;
	icdf_normal(in2, out2, 0, 1);
	printf("  input:     ");
	in2.print();
	printf("  output:    ");
	out2.print();
	printf("  should be  [ -0.4071 ]\n\n");

	//make cdf
	Float1dArray cdf_vals(100);
	Uint1dArray range(50);
	for (int i = 0; i < range.x; ++i)
		range[i] = rrand(100);
	make_cdf(range, cdf_vals, 100);
	printf("cdf_vals = ");
	cdf_vals.print(cdf_vals.x);

}

/**
 *
 */
void test_genescore(Float2dArray &normgenes, Float2dArray &geneLabelArray,
		Double1dArray &zscore_genes) {

	printf("\n===test_genescore====\n\n");

	// test vs matlab
	Float2dArray class1(2, 10);
	Float2dArray class2(2, 10);
	Double1dArray zscore(2);

	for (int i = 0; i < class1.x; ++i) {
		for (int j = 0; j < class1.y; ++j) {
			class1(i, j) = frand();
			class2(i, j) = frand();
		}
	}

	printf("class1 = ");
	class1.print();
	printf("class2 = ");
	class2.print();
	genescore(class1, class2, zscore);
	printf("c_zscore = ");
	zscore.print();
	printf("should_be = [ 0.707107 -0.707107 ]\n");

	// Have to fix genes
	for (unsigned int i = 0; i < normgenes.x; ++i) {
		for (unsigned int j = 0; j < normgenes.y; ++j) {
			normgenes(i, j) = log2(normgenes(i, j) + 4);
		}
	}

	Float1dArray geneMean(normgenes.x);
	Float1dArray geneStddev(normgenes.x);
	mean(normgenes, geneMean);
	stddev(normgenes, geneMean, geneStddev);
	for (int i = 0; i < normgenes.x; ++i) {
		for (int j = 0; j < normgenes.y; ++j) {
			normgenes(i, j) = (normgenes(i, j) - geneMean[i]) / geneStddev[i];
		}
	}

//	MatlabFile mf(WRITE);
//	mf.open("zscore-cpu-after-fix.mat");
//	mf.writeArray("normgenes", &normgenes);

// Separate classes
	Float1dArray classes(2);
	classes[0] = 1.0f;
	classes[1] = 2.0f;
	Uint1dArray index_class1(geneLabelArray.x);
	Uint1dArray index_class2(geneLabelArray.x);
	populateClasses(geneLabelArray, index_class1, index_class2, classes);

	Float2dArray normgenes_class1(normgenes.x, index_class1.length());
	classsample(normgenes, index_class1, normgenes_class1);

	Float2dArray normgenes_class2(normgenes.x, index_class2.length());
	classsample(normgenes, index_class2, normgenes_class2);

//	mf.writeArray("normgenes_class1", &normgenes_class1);
//	mf.writeArray("normgenes_class2", &normgenes_class2);

// Zscore
	genescore(normgenes_class1, normgenes_class2, zscore_genes);

	Float2dArray zscore_cpu(zscore_genes.x, 1);
	for (int i = 0; i < zscore_cpu.x; ++i) {
		float temp = float(zscore_genes[i]);
		zscore_cpu(i, 0) = temp;
	}
//	mf.writeArray("zscore_cpu_after_fix", &zscore_cpu);
//	mf.close();

}

/**
 *
 */
void test_g_conn() {

	printf("\n===test_g_conn====\n\n");

	Uint2dArray ppi(8, 2);
	for (int i = 0; i < ppi.x; ++i) {
		ppi(i, 0) = i;
		ppi(i, 1) = i + 1;
	}
	printf("ppi = ");
	ppi.print();

// Load gene ids
	Uint1dArray geneid(2);
	for (int i = 0; i < geneid.x; ++i) {
		geneid[i] = rrand(ppi.x);
	}
	printf("geneid = ");
	geneid.print();

// Run
	Uint1dArray igconn(2);
	printf("index = 0\n");
	g_conn(ppi, geneid, 0, igconn);
	printf("  igconn = ");
	igconn.print();
	printf("\n");

	printf("index = 1\n");
	g_conn(ppi, geneid, 1, igconn);
	printf("  igconn = ");
	igconn.print();
	printf("\n");

}

/**
 *
 */
void test_rands() {
	printf("\n===test_rands====\n\n");

	// Run brand 1000 times...should be 50% yes and 50% no.
	unsigned int true_times = 0;
	unsigned int false_times = 0;
	for (int i = 0; i < 1000; ++i) {
		bool val = brand();
		if (val)
			true_times++;
		else
			false_times++;
	}
	printf("brand 1000 iterations: true %i false %i\n", true_times,
			false_times);

	// Run brand 1000 times...should be 20% yes and 80% no.
	true_times = 0;
	false_times = 0;
	for (int i = 0; i < 1000; ++i) {
		bool val = brand(0.2);
		if (val)
			true_times++;
		else
			false_times++;
	}
	printf("brand 1000 iterations: true %i false %i\n", true_times,
			false_times);

}

/**
 *
 */
void test_scoreNetwork(Uint2dArray &ppi, Uint1dArray &geneid,
		Double1dArray & zscore) {

	printf("\n===test_scoreNetwork====\n\n");

	// Create candidate network
	Uint1dArray cand_network_id(ppi.x);
	cand_network_id.shrink(0);
	cand_network_id.append(2099);
	cand_network_id.append(627);

	printf("cand_network_id = ");
	cand_network_id.print(cand_network_id.x);

	// Evaluate
	double netscore = mrfnetscore(geneid, cand_network_id, zscore, ppi);
	printf("Netscore is: %f, should be -0.3238\n", netscore);

	// Second network
	cand_network_id.shrink(0);
	cand_network_id.append(2099);
	cand_network_id.append(6908);
	cand_network_id.append(1387);
	cand_network_id.append(6197);
	printf("cand_network_id = ");
	cand_network_id.print(cand_network_id.x);

	netscore = mrfnetscore(geneid, cand_network_id, zscore, ppi);
	printf("Netscore is: %f, should be -0.3745\n", netscore);

	// Third network
	cand_network_id.shrink(0);
	cand_network_id.append(2099);
	cand_network_id.append(8841);
	cand_network_id.append(10014);
	printf("cand_network_id = ");
	cand_network_id.print(cand_network_id.x);

	netscore = mrfnetscore(geneid, cand_network_id, zscore, ppi);
	printf("Netscore is: %f, should be -3.9732\n", netscore);

	// Fourth network
	cand_network_id.shrink(0);
	cand_network_id.append(190);
	cand_network_id.append(595);
	cand_network_id.append(672);
	cand_network_id.append(1387);
	cand_network_id.append(1956);
	cand_network_id.append(2099);
	cand_network_id.append(2100);
	cand_network_id.append(3065);
	cand_network_id.append(3066);
	cand_network_id.append(3265);
	cand_network_id.append(3320);
	cand_network_id.append(5295);
	cand_network_id.append(5469);
	cand_network_id.append(5594);
	cand_network_id.append(5595);
	cand_network_id.append(5604);
	cand_network_id.append(5894);
	cand_network_id.append(6714);
	cand_network_id.append(6908);
	cand_network_id.append(7157);
	cand_network_id.append(8204);
	cand_network_id.append(8841);
	cand_network_id.append(9112);
	cand_network_id.append(9612);
	cand_network_id.append(9759);
	cand_network_id.append(10014);
	cand_network_id.append(10891);
	cand_network_id.append(23013);
	cand_network_id.append(27043);
	cand_network_id.append(83933);
	printf("cand_network_id = ");
	cand_network_id.print(cand_network_id.x);

	netscore = mrfnetscore(geneid, cand_network_id, zscore, ppi);
	printf("Netscore is: %f, should be -2.6717\n", netscore);

}

void test_randsample(Float2dArray &genes, Float2dArray &labels) {

	printf("\n===test_randsample====\n\n");

	// get classes
	Float1dArray classes(2);
	classes[0] = 1.0f;
	classes[1] = 2.0f;
	Uint1dArray index_class1(labels.x);
	Uint1dArray index_class2(labels.x);
	populateClasses(labels, index_class1, index_class2, classes);

	printf("labels are: %i x %i\n", labels.x, labels.y);

	// Rand sample
	printf(
			"class 1: valid options, followed by randomized selection with resample...\n");
	index_class1.print();
	Float2dArray rand_class1(genes.x, index_class1.length());
	randsample(genes, index_class1, rand_class1);

	printf(
			"class 2: valid options, followed by randomized selection with resample...\n");
	index_class2.print();
	Float2dArray rand_class2(genes.x, index_class2.length());
	randsample(genes, index_class2, rand_class2);

}

void test_netcand(Uint2dArray &ppi, Uint1dArray &geneids,
		const unsigned int &seedgene) {

	printf("\n===test_netcand====\n\n");

	Uint1dArray prev_network_id_Array(ppi.x);
	prev_network_id_Array[0] = seedgene;
	prev_network_id_Array.shrink(1);
	Uint1dArray prev_network_distance_Array(ppi.x);
	prev_network_distance_Array[0] = 0;
	prev_network_distance_Array.shrink(1);

	Uint1dArray candidate_network_id_Array(ppi.x);
	Uint1dArray candidate_network_distance_Array(candidate_network_id_Array.x);
	candidate_network_id_Array.shrink(0);
	candidate_network_distance_Array.shrink(0);

	// Call this in a loop because netcand has rand in it
	for (int i = 0; i < 4; ++i) {
		// Generate the candidate nodes in the network with one step further
		netcand(ppi, geneids, prev_network_id_Array,
				prev_network_distance_Array, 1, candidate_network_id_Array,
				candidate_network_distance_Array);

		candidate_network_id_Array.print();
		candidate_network_distance_Array.print();
	}
}

void test_mrfsearchnet(Uint2dArray &ppi, Uint1dArray &geneids,
		const unsigned int &seedgene, Double1dArray &zscore_genes) {

	printf("\n===test_mrfsearchnet====\n\n");

	for (int i = 0; i < 20; ++i) {
		Uint1dArray sub_network(ppi.x);
		double sub_netscore;

		mrfsearchnet(ppi, geneids, zscore_genes, seedgene, 2, 1.0f, sub_network,
				sub_netscore);

		printf("subnetwork score: %f\n", sub_netscore);
		printf("sub_network = ");
		sub_network.print(sub_network.x);
		fflush(stdout);
	}

}
