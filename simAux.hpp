/**
 * @file simAux.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Some helper functions for the simulation
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

#ifndef SIMAUX_HPP
#define SIMAUX_HPP
#include "State.hpp"
#include "Key.hpp"
#include <iomanip>

template<class T>
using CNum = std::complex<T>;

/**
 * @brief Generate the rotation unitaries used for the wave-plates including the errors
 * 
 * @tparam V Value-type, cf. State
 * @tparam R Real-type, cf. State
 * @param angleErrs Rotation errors
 * @return std::array<std::vector<V>, 15> Unitaries for the rotations in line-form
 */
template<class V, class R>
std::array<std::vector<V>, 15> genRotationsBasic(const std::vector<R>& angleErrs){
    std::array<std::vector<V>, 15> Ops;
    std::vector<V> U = {0, 0, 0, 0};
    V emix, epix;
    for (int i= 0;i<15;i++){
		emix = std::cos((45.0+angleErrs[i])/180.0*numbers::pi);
        epix = std::sin((45.0+angleErrs[i])/180.0*numbers::pi);
        U[0] = emix;
        U[1] = epix;
        U[2] = epix;
        U[3] = -emix;
        Ops[i] = U;}
    return Ops;
}

/**
 * @brief Writes the results of the simulation to a file
 * 
 * @tparam R Real number type that should be used, e.g. float.
 * @param LossPositions Positions of loss in the circuit
 * @param doublePrep Spatial&Polarization modes with two-photon preparation
 * @param angleErrs Rotation error for wave plated
 * @param ovl pairwise overlap of wave functions
 * @param path Path-prefix where to save
 * @param rank Rank of the process, used for saving
 * @param res Result of the simulation
 */
template<class R>
inline void write(const std::vector<int>& LossPositions, const std::vector<int>& doublePrep, const std::vector<R>& angleErrs, const R& ovl, const std::string& path, int rank, const std::vector<R>& res){
	std::ostringstream sstream;
	std::ofstream myfile;
	sstream << ovl << " ";
	for (R i: angleErrs) sstream << i << "|";
	sstream <<" ";
	for (int i: doublePrep) sstream << i << "|";
	sstream <<" ";
	for (int i: LossPositions) sstream << i << "|";
	sstream <<" ";
	for (R i: res) sstream << std::setprecision(12) << i << " ";
	sstream <<"\n";
	std::string r = std::to_string(rank);
    myfile.open(path+r+".txt", std::ios_base::app); 
	myfile << sstream.str(); 
	myfile.close();
}

/**
 * @brief Given to float-vectors that parametrize wave functions, it returs the overlap, which is trivially saved in the first vector.
 * 
 * @param V First vector, at least len 2 with identifier on 0 and overlap on 1
 * @param W Second vector, at least len 1 with identifier on 0
 * @return float Overlap of the wave functions
 */
float trivOvlF(const std::vector<float>& V, const std::vector<float>& W){
    if (V[0]!= W[0]) 
        return V[1]; 
    return 1.0;
}

/**
 * @brief Applies loss on modes in S if the current position pos is in lossPos
 * 
 * @tparam K Key type that should be used. 
 * @tparam V Amplitude type that should be used, e.g. float or std::complex<float>
 * @tparam R Real number type that should be used, e.g.
 * @param S State to apply loss on
 * @param pos current position in circuit
 * @param modes modes affected by loss
 * @param lossPos potitions in circuit where loss happens
 */
template<class K, class V, class R>
void detloss(State<K, V, R>& S, int pos, const std::vector<int>& modes, const std::vector<int>& lossPos){
    if (std::find(lossPos.cbegin(), lossPos.cend(), pos)!= lossPos.cend())
        S.loss(modes);
}

#endif