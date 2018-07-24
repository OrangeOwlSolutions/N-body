#include <thrust\host_vector.h>
#include <thrust\unique.h>
#include <thrust\merge.h>
//#include <thrust/execution_policy.h>

/********/
/* MAIN */
/********/
int main() {

	thrust::host_vector<int> h_a(6);
	thrust::host_vector<int> h_b(7);
	thrust::host_vector<int> h_c(13);
	h_a[0] = 1; h_a[1] = 3; h_a[2] = 5; h_a[3] = 7; h_a[4] = 9; h_a[5] = 11;
	h_b[0] = 1; h_b[1] = 1; h_b[2] = 2; h_b[3] = 3; h_b[4] = 5; h_b[5] = 8; h_b[6] = 13;

	thrust::merge(h_a.begin(), h_a.end(), h_b.begin(), h_b.end(), h_c.begin());

	//h_b.insert(h_b.end(), h_a.begin(), h_a.end());
	
	printf("After merging\n");
	for (int k = 0; k < 13; k++) printf("h_c[%d] = %d\n", k, h_c[k]);

	auto new_end = thrust::unique(h_c.begin(), h_c.end());
	
	int new_elements = new_end - h_c.begin();
	
	printf("\nAfter unique\n");
	for (int k = 0; k < new_elements; k++) printf("h_c[%d] = %d\n", k, h_c[k]);

	return 0;
}
