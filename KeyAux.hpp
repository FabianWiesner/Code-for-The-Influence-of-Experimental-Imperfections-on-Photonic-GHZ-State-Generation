/**
 * @file KeyAux.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Some helper functions for Key.hpp.
 * @version 0.1
 * @date 2024-06-06
 * 
 * @copyright None, this file only contains trivial common knowledge code. 
 * 
 */

 /*
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
 */
 
#ifndef KEYAUX_HPP
#define KEYAUX_HPP
#include <vector>
#include <math.h> 
#include <string>
#include <algorithm>
#include <sstream>
#include <complex>

/**
 * @brief Pre-computed faculties up to 12! as this is the maximal number of photons in the circuit.
 * 
 */
const int FACUT12[13]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

/**
 * @brief Short handle for pi.
 * 
 */
namespace numbers{
 const double pi = 3.14159265358979323846;
}

/**
 * @brief Recrusive definition of the faculty. Up to 12, this is just a look-up.
 * 
 * @param i Input to compute the faculty of.
 * @return int The fuculty of i, i.e. i!.
 */
int facut(int i){
	if (i<= 12){return FACUT12[i];}
	else std::cout << i << " facut" << std::endl;
	return i*facut(i-1);
}

/**
 * @brief Returns the binomial coefficient n over k. Shouldn't be used for large numbers and is only used for n<=12.
 * 
 * @tparam R Real-type for the output. 
 * @tparam I Integer number type used for input.
 * @param n Upper number of the binomial coefficient, 0<=n<=12.
 * @param k Lower number of the binomial coefficient, 0<=k<=n.
 * @return R n over k.
 */
template<class R, class I>
R binomialCoeff(I n, I k)
{	
    return (R) (facut((int) n)/(facut((int) k)*facut((int) (n-k))));
}
#endif