/**
 * @file simFid.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Defines the function used for the simulation of GHZ state generation with imperfections.
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

#ifndef SIMFID_HPP
#define SIMFID_HPP
#include "State.hpp"
#include "Key.hpp"
#include "simAux.hpp"
#include <iomanip>




/**
 * @brief Cleans input such that only the keys with the same DMode in one of the parts of the GHZ state remain.
 * 
 * @tparam K Key-type, cf. State
 * @tparam V Value-type, cf. State
 * @tparam R Real-type, cf. State
 * @param S State to compute to overlap for
 * @param a phase for GHZ state - either 1 or -1
 */
template<class K, class V, class R>
inline void cleanOvlGHZ(State<K, V, R>& S, int a){
    V f = 1.0/std::sqrt(2), f2 = ((V) a)/std::sqrt(2);
    S.sameDModeDel({{0, 6, 10}, {1, 7, 11}}, {f, f2}); //Projection where all 
}

/**
 * @brief Computes the fidelity for pre-computed data and writes it to a file.
 * 
 * @param PreData Vector of Arrays (one for every input combination of DModes) of the remaining states after the measurement already projected onto the spatial and polarization modes of the GHZ state
 * @param Compl Same data structure as PreData. These contains the complement of the states after the measurement, i.e. those parts orthogonal to the GHZ state.
 * @param ovl Pairwise overlap
 * @param angErrs Angle errors in the setup
 * @param doublePrep Events/positions of double-preparation
 * @param lossPos Events/positions of loss
 * @param path The path to the folder where the data should be saves.
 * @param rank The rank of the current process. Used for saving the data to the file.
 */
inline void fid(const std::vector<std::array<State<Key<int>, float, float>, 8>>& PreData, const std::vector<std::array<State<Key<int>, float, float>, 8>>& Compl, const float& ovl, const std::vector<float>& angErrs, const std::vector<int>& doublePrep, const std::vector<int>& lossPos, const std::string& path, int rank){
    State<Key<int>, float, float> S, S2, S3;
    S.set(&trivOvlF);
    S.set((float) std::pow(10, -8));
    std::array<State<Key<int>, float, float>, 8> StV, StV2;
    for (int i=0; i<8; i++) {StV[i] = S; StV2[i] = S;};
    for (int i=0; i<6; i++){
            S.addPhoton({(float) i, ovl}, 2*i, 1);
    }
    int i=0;
    boost::container::flat_map<int, int> M({{0, 1}, {1, -1}, {2, -1}, {3, 1}, {4, -1}, {5, 1}, {6, 1}, {7, -1}});
    for (typename State<Key<int>, float, float>::iterator it = S.begin(); it!=S.end();it++){
        for (int j = 0; j<8; j++){
            S2 = PreData[i][j];
            cleanOvlGHZ(S2, M[j]);
            StV[j].add(S2, it->second);
            StV2[j].add(PreData[i][j], it->second);
            StV2[j].add(Compl[i][j], it->second);}
        i++;
    }
    std::vector<float> res;
    for (int i=0; i<8; i++) {
        res.push_back(std::pow(StV2[i].norm(), 2));
        res.push_back(std::pow(StV[i].norm(), 2)/std::pow(StV2[i].norm(), 2));}
    write(lossPos, doublePrep, angErrs, ovl, path, rank, res);
}

/**
 * @brief The photonic circuit we considered to create a GHZ state.
 * 
 * @tparam K Key-type, cf. State
 * @tparam V Value-type, cf. State
 * @tparam R Real-type, cf. State
 * @param S State to perform the circuit on.
 * @param lossPos Positions where loss happens
 * @param apl Rotations as unitaries repr. as single line unitaries 
 */
template<class K, class V, class R>
inline void circuitFid(State<K, V, R>& S, const std::vector<int>& lossPos, const std::array<std::vector<V>, 15>& apl){
    State<K, V, R> S2;
    for (int i=0; i<6; i++) detloss(S, i, {2*i}, lossPos);
    for (int i=0; i<6; i++) detloss(S, i+6, {2*i}, lossPos);
    for (int i=0; i<6; i++)
        S.apply(apl[i], {2*i, 2*i+1});
    for (int i=0; i<6; i++) detloss(S, i+12, {2*i, 2*i+1}, lossPos);
    for (int i=0; i<3; i++)
        S.swap(4*i+1, 4*i+3);
    for (int i=0; i<6; i++) detloss(S, i+18, {2*i, 2*i+1}, lossPos);
    for (int i=0; i<6;i++)
        S.apply(apl[i+6], {2*i, 2*i+1});
    detloss(S, 24, {2, 3}, lossPos);
    detloss(S, 25, {4, 5}, lossPos);
    S.swap(3, 5);
    detloss(S, 26, {4, 5}, lossPos);
    detloss(S, 27, {8, 9}, lossPos);
    S.swap(5, 9);
    detloss(S, 28, {2, 3}, lossPos);
    detloss(S, 29, {4, 5}, lossPos);
    detloss(S, 30, {8, 9}, lossPos);
    S.apply(apl[12], {2, 3});
    S.apply(apl[13], {4, 5});
    S.apply(apl[14], {8, 9});
    detloss(S, 32, {2, 3}, lossPos);
    detloss(S, 33, {4, 5}, lossPos);
    detloss(S, 35, {8, 9}, lossPos);
}

/**
 * @brief Maps a the perfectly distinguishable configuration to a partially distinguishability conf.
 * 
 * @param K Key that encodes the desired distinguishability configuration
 * @param SVec States that should be used for the result, completly overlapping with GHZ in spatial and polarization and acceptable measurement results.
 * @param comp States that are used for normalization, i.e. part of the measurement result but orthogonal to GHZ, as map is not norm preserving
 * @param S Part of the state that is orthogonal to all acceptable measurement results
 */
void collapseRenorm(const Key<int>& K, std::array<State<Key<int>, float, float>, 8>& SVec, std::array<State<Key<int>, float, float>, 8>& comp, State<Key<int>, float, float>& S){
    float n=0.0;
    for (int i = 0; i<8; i++){
        SVec[i].collapse(K);
        comp[i].collapse(K);
        n += std::pow(SVec[i].norm(),2);
        n += std::pow(comp[i].norm(),2);
    }
    S.collapse(K);
    n += std::pow(S.norm(),2);
    if (n!=0.0){
    n = 1/std::sqrt(n);
    for (int i = 0; i<8; i++){
        SVec[i].mul(n);
        comp[i].mul(n);
    }}
}

/**
 * @brief Computes the fidelity for a given parameters.
 * 
 * @param ovls Overlaps, for all of them the fidelity is computed
 * @param doublePrep Spatial modes with two-photon preparation
 * @param lossPos Positions where loss happens
 * @param angErrs Rotation-errors for wave-plates
 * @param apl Rotations as unitaries repr. as single line unitaries (already including the rotation errors)
 * @param path Pathsuffix where to save the outcome
 * @param rank Rank of the process (used for saving the outcome)
 */
void fidsim(const std::vector<float>& ovls, const std::vector<int>& doublePrep, const std::vector<int>& lossPos, const std::vector<float>& angErrs, const std::array<std::vector<float>, 15>& apl, const std::string& path, int rank){
    State<Key<int>, float, float> SFullDist, SKeyIter, STemp, SComplement;

    SFullDist.set(12);
    SFullDist.set(&trivOvlF);
    for (int i=0; i<6; i++){
        if (std::find(doublePrep.begin(), doublePrep.end(), i)!=doublePrep.end())
            SFullDist.addPhoton({(float) i, 0.0}, 2*i, 2);
        else
            SFullDist.addPhoton({(float) i, 0.0}, 2*i, 1);
    }
    circuitFid(SFullDist, lossPos, apl);

    std::array<std::pair<State<Key<int>, float, float>, State<Key<int>, float, float>>, 8> SV;
    std::vector<std::array<State<Key<int>, float, float>, 8>> SarS, SarComp;
    std::vector<boost::container::flat_map<int, int>> 
    g = {{{0, 1}, {1, 0}, {6, 1}, {7, 0}, {10, 1}, {11, 0}}, {{0, 0}, {1, 1}, {6, 0}, {7, 1}, {10, 0}, {11, 1}}};
    std::vector<std::vector<int>> occModes={{0, 6, 10}, {1, 7, 11}};
    std::array<State<Key<int>, float, float>, 8> SVec, compVec, SVecTemp, compVecTemp;
    std::vector<boost::container::flat_map<int,int>> FirstMeasTargets;
    for (int p1: {0, 1})
        for (int p2: {0, 1})
            for (int p4: {0, 1}){
                STemp = SFullDist;
                SComplement = STemp.overlapWithFilter({{2, 1-p1}, {3, p1}, {4, 1-p2}, {5, p2}, {8, 1-p4}, {9, p4}}, occModes, g);
                FirstMeasTargets.push_back({{2, 1-p1}, {3, p1}, {4, 1-p2}, {5, p2}, {8, 1-p4}, {9, p4}});
                SVec[p4+2*p2+4*p1] = std::move(STemp);
                compVec[p4+2*p2+4*p1] = std::move(SComplement);
                }
    SFullDist.overlapCompl(FirstMeasTargets, occModes);
    SKeyIter.set(12);
    SKeyIter.set(&trivOvlF);
    for (int i=0; i<6; i++){
        if (std::find(doublePrep.begin(), doublePrep.end(), i)!=doublePrep.end())
            SKeyIter.addPhoton({(float) i, 0.7}, 2*i, 2);
        else
            SKeyIter.addPhoton({(float) i, 0.7}, 2*i, 1);
    }
    for (typename State<Key<int>, float, float>::iterator it = SKeyIter.begin(); it!=SKeyIter.end();it++){
        STemp = SFullDist;
        SVecTemp = SVec;
        compVecTemp = compVec;
        collapseRenorm(it->first, SVecTemp, compVecTemp, temp);
        SarS.push_back(std::move(SVecTemp));
        SarComp.push_back(std::move(compVecTemp));
    }
    for (float o : ovls)
        fid(SarS, SarComp, o, angErrs, doublePrep, lossPos, path, rank);
}

/**
 * @brief This function iterates over most likely 10214 combinations of loss and two-photon creation and saves the fidelity and the probailities for all of them.
 * 
 * @param ovls Overlaps, for all of them the fidelity is computed
 * @param angErrs Rotation-errors for wave-plates
 * @param path Pathsuffix where to save the outcome
 * @param global_lower Lower end of the intervall, between 0 and 10214
 * @param global_upper Upper end of the intervall, between 0 and 10214
 * @param rank_off Offset for the rank, which is used for saving the result
 * @param rank Rank of the process (used for saving the outcome)
 * @param size Number of processes
 * @param shuffle_path Path to a file where all 10214 combinations are shuffeled
 */
void schedulerGHZshuffled(const std::vector<float>& ovls, std::vector<float>& angErrs, std::string path, int global_lower, int global_upper, int rank_off, int rank, int size, std::string shuffle_path){
    int n=global_upper-global_lower;
    int count = 0;
    std::vector<int> todo;
    std::ifstream infile(shuffle_path);
    int i = 0, p;
    while (infile >> p){
        if (i%size==rank && p>=global_lower && p<global_upper){
            todo.push_back(p);
        }
        i++;
    }
    std::vector<int> p2 = {}, pl={};
    std::array<std::vector<float>, 15> apl = genRotationsBasic<float, float>(angErrs);
    if (std::find(todo.begin(), todo.end(), count)!= todo.end()) {fidsim(ovls, p2, pl, angErrs, apl, path, rank+rank_off);
    ;std::cout << count << std::endl;} //no error
	if (count>global_upper) return;
    for (int i=0;i<6;i++){ //single error comb.
    	p2 = {i};
    	for (int j=0;j<37;j++){
    		pl = {j};
    		count++;
    		if (std::find(todo.begin(), todo.end(), count)!= todo.end()) {fidsim(ovls, p2, pl, angErrs, apl, path, rank+rank_off);
    ;std::cout << count << std::endl;}
			if (count>global_upper) return;
    	}
    }

    for (int i0=0;i0<5;i0++){ //two error combs.
    	for (int i1=i0+1;i1<6;i1++){
    		p2 = {i0, i1};
    		for (int j0=0;j0<36;j0++){
    			for (int j1=j0+1;j1<37;j1++){
    				pl = {j0, j1};
    				count++;
    				if (std::find(todo.begin(), todo.end(), count)!= todo.end()) {fidsim(ovls, p2, pl, angErrs, apl, path, rank+rank_off);
    ;std::cout << count << std::endl;}
					if (count>global_upper) return;
    			}
    		}
    	}
    }

    for (int i0=0;i0<4;i0++){
    	for (int i1=i0+1;i1<5;i1++){
    		for (int i2=i1+1;i2<6;i2++){
    			p2 = {i0, i1, i2};
    			for (int j0=0;j0<35;j0++){
    				for (int j1=j0+1;j1<36;j1++){
    					for (int j2=j1+1;j2<37;j2++){
    						pl = {j0, j1, j2};
		    				count++;
		    				if (std::find(todo.begin(), todo.end(), count)!= todo.end()) {fidsim(ovls, p2, pl, angErrs, apl, path, rank+rank_off);
    ;std::cout << count << std::endl;}
							if (count>global_upper) return;
    					}
    				}
    			}    		
    		}
    	}
    }

}

#endif