#include "TMalign_wrapper.h"

using namespace std;
using namespace TMalignC;

int main(int argc, char *argv[]){
	TMalign_wrapper tmalign;
	tmalign.run_main(argc, argv);

	return 0;
}
