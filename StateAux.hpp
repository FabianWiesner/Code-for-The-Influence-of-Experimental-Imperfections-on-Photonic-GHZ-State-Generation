/**
 * @file StateAux.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Some helper functions for State.hpp
 * @version 0.1
 * @date 2024-06-06
 * 
 * @copyright Copyright (c) 2024, provided under CC BY-NC 4.0. license
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

#ifndef STATEAUX_HPP
#define STATEAUX_HPP
#include <vector>
#include <math.h> 
#include <string>
#include <algorithm>
#include <sstream>
#include <complex>

/**
 * @brief Performs the usual complex conjugate.
 * 
 * @tparam T Number type
 * @param t Number to perform the complex conj on.
 * @return T Complex conjugate
 */
template<class T>
T conj(T t){
    return (T) std::conj(t);
}

/**
 * @brief Template specialization. Complex conjugate is identity for floats.
 * 
 * @param t Number to perform the complex conj on.
 * @return float Complex conjugate
 */
template<>
float conj<float>(float f){return f;};

/**
 * @brief Template specialization. Complex conjugate is identity for doubles.
 * 
 * @param t Number to perform the complex conj on.
 * @return double Complex conjugate
 */
template<>
double conj<double>(double f){return f;};

/**
 * @brief Template specialization. Complex conjugate is identity for ints.
 * 
 * @param t Number to perform the complex conj on.
 * @return int Complex conjugate
 */
template<>
int conj<int>(int f){return f;};

/**
 * @brief Computes the inner product of a OWF, and WF according to get_ovlp (cf. State).
 * 
 * @tparam Val Value-type, cf. State.
 * @tparam Real Real-type, cf. State
 * @param b Orthogonal wave function (OWF).
 * @param wf Wave function WF
 * @param get_ovlp Function that provides the overlap given two wave functions.
 * @param waves Wave functions that are used to express b. b provides the coefficients.
 * @return Val Overlap of b and wf. 
 */
template<class Val, class Real>
Val ovlpH(const std::vector<Val>& b, const std::vector<Real>& wf, Val (*get_ovlp)(const std::vector<Real>&, const std::vector<Real>&), const std::vector<std::vector<Real>>& waves){
    Val res = 0;
    for (int i = 0; i<b.size(); i++){
        res += conj(b[i])*get_ovlp(waves[i], wf);
    }
    return res;
}

/**
 * @brief Computes the inner product of two proto-orthogonal wave functions according to get_ovlp. Used in addBasisElem() (cf. State). 
 * 
 * @tparam Val Value-type, cf. State.
 * @tparam Real Real-type, cf. State.
 * @param b proto-orthogonal wave function
 * @param c proto-orthogonal wave function
 * @param Function that provides the overlap given two wave functions.
 * @param waves Wave functions that are used to express b. b provides the coefficients.
 * @return Val Overlap of b and c.
 */
template<class Val, class Real>
Val ovlp(const std::vector<Val>&  b, const std::vector<Val>&  c, Val (*get_ovlp)(const std::vector<Real>&, const std::vector<Real>&), const std::vector<std::vector<Real>>& waves){
    Val res = 0;
    for (int i = 0; i<c.size(); i++){
        res += c[i]*ovlpH<Val, Real>(b, waves[i], get_ovlp, waves);
    }
    return res;
}


#endif