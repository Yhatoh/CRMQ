// c++ includes
#include <iostream>
#include <string>

// sdsl includes
#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_mcl.hpp>
#include "sdsl/util.hpp"

// local includes
#include "scst.hpp"
//#include "cbp.hpp"

int main() {
        
  //std::string bp_example = "((()(()()()))(()()())(()()())(()(()()())))";
  //std::string bp_example = "((((())())((())())())((())())())";
  std::string bp_example = "((((())())((())())(())())((())())(())())";

  sdsl::bit_vector bv = sdsl::bit_vector(bp_example.length(), 0);
  sdsl::select_support_mcl<0> select_bv;
  std::cout << bv.size() << std::endl;
  for(uint64_t i = 0; i < bv.size(); i++) {
    bv[i] = (bp_example[i] == '(' ? 1 : 0);
  }
  sdsl::util::init_support(select_bv, &bv);

  scst<4> cbv(bv, 2);
  std::cout << "Amount of Maximal Identical Subtrees: " << cbv.max_trees() << std::endl;
  std::cout << "Length compressed bp: " << cbv.nodes() * 2 << std::endl;
  std::cout << "Amount of nodes: " << cbv.nodes() << std::endl;

  std::cout << "LCA QUERIES" << std::endl;
  std::cout << " CASE 1: " << cbv.lca(select_bv(1), select_bv(3)) << std::endl;
  std::cout << " CASE 2: " <<  cbv.lca(select_bv(5), select_bv(9)) << std::endl;
  std::cout << " CASE 3: " <<  cbv.lca(select_bv(5), select_bv(11)) << std::endl;
  std::cout << " CASE 4: " <<  cbv.lca(select_bv(11), select_bv(13)) << std::endl;

  std::cout << "INORDER QUERIES" << std::endl;
  std::cout << " u = 5: " << cbv.inorder(5, 0) << std::endl;
  std::cout << " u = 6: " << cbv.inorder(6, 0) << std::endl;
  std::cout << " u = 8: " << cbv.inorder(8, 0) << std::endl;
  std::cout << " u = 9: " << cbv.inorder(9, 0) << std::endl;
  std::cout << " u = 13: " << cbv.inorder(13, 1) << std::endl;
  std::cout << " u = 14: " << cbv.inorder(14, 2) << std::endl;
  std::cout << " u = 16: " << cbv.inorder(16, 3) << std::endl;
  std::cout << " u = 17: " << cbv.inorder(17, 0) << std::endl;
  std::cout << " u = 19: " << cbv.inorder(19, 0) << std::endl;
  std::cout << " u = 20: " << cbv.inorder(20, 0) << std::endl;
  std::cout << " u = 24: " << cbv.inorder(24, 1) << std::endl;
  std::cout << " u = 25: " << cbv.inorder(25, 2) << std::endl;
  std::cout << " u = 27: " << cbv.inorder(26, 3) << std::endl;
  std::cout << " u = 28: " << cbv.inorder(27, 4) << std::endl;
  std::cout << " u = 30: " << cbv.inorder(30, 5) << std::endl;

  std::cout << "SELECT 0 QUERIES" << std::endl;
  std::cout << " select(B, 1) expected: " << select_bv(1) << " obtained: " << cbv.select0(1) << std::endl;
  std::cout << " select(B, 2) expected: " << select_bv(2) << " obtained: " << cbv.select0(2) << std::endl;
  std::cout << " select(B, 3) expected: " << select_bv(3) << " obtained: " << cbv.select0(3) << std::endl;
  std::cout << " select(B, 4) expected: " << select_bv(4) << " obtained: " << cbv.select0(4) << std::endl;
  std::cout << " select(B, 5) expected: " << select_bv(5) << " obtained: " << cbv.select0(5) << std::endl;
  std::cout << " select(B, 6) expected: " << select_bv(6) << " obtained: " << cbv.select0(6) << std::endl;
  std::cout << " select(B, 7) expected: " << select_bv(7) << " obtained: " << cbv.select0(7) << std::endl;
  std::cout << " select(B, 8) expected: " << select_bv(8) << " obtained: " << cbv.select0(8) << std::endl;
  std::cout << " select(B, 9) expected: " << select_bv(9) << " obtained: " << cbv.select0(9) << std::endl;
  std::cout << " select(B, 10) expected: " << select_bv(10) << " obtained: " << cbv.select0(10) << std::endl;
  std::cout << " select(B, 11) expected: " << select_bv(11) << " obtained: " << cbv.select0(11) << std::endl;
  std::cout << " select(B, 12) expected: " << select_bv(12) << " obtained: " << cbv.select0(12) << std::endl;
  std::cout << " select(B, 13) expected: " << select_bv(13) << " obtained: " << cbv.select0(13) << std::endl;
  std::cout << " select(B, 14) expected: " << select_bv(14) << " obtained: " << cbv.select0(14) << std::endl;
  std::cout << " select(B, 15) expected: " << select_bv(15) << " obtained: " << cbv.select0(15) << std::endl;
  return 0;
}
