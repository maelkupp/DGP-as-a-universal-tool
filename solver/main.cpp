#include "dgp_bnb.h"


int main(int argc, char* argv[]){
    if(argc != 2){
        std::cerr << "Usage " << argv[0] << "<dat_file.dat>";
        return 1;
    }

    auto error = solve_minerr_dgp1(argv[1]);
    std::cerr << "Error " << error << "\n" << std::endl;
    std::cout << error;
}