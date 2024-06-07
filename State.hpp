/**
 * @file State.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Defines the central data structure State.
 * @version 0.1
 * @date 2024-06-06
 * 
 * @copyright Copyright (c) 2024, provided under CC BY-NC 4.0. 
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
 
#ifndef STATE_HPP
#define STATE_HPP

#include <vector>
#include <math.h> 
#include <string>
#include <algorithm>
#include <sstream>
#include <complex>
#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <initializer_list>
#include <array>
#include <utility>
#include <functional>
#include <stdio.h>
#include <boost/container/flat_map.hpp>
#include <boost/algorithm/string.hpp>
#include "StateAux.hpp" 

/**
 * @brief Definition of the data structure used to represent states. Inherits from boost::container::flat_map<Key, Val>.
 * 
 * @tparam Key Key type used in the State 
 * @tparam Val Amplitude type used in the State, e.g. float or std::complex<float> 
 * @tparam Real Real number type that should be used, e.g. float. 
 */
template<class Key, class Val, class Real>
class State : public boost::container::flat_map<Key, Val>{
    /**
     * @brief Short handle for the parent type that provides the actual data structure.
     * 
     */
    using Par = boost::container::flat_map<Key, Val>;

    /**
     * @brief Wave functions are parametrized as vectors of Reals. The overlap function needs to interprete these as intended by the user.
     * 
     */
    using WF = std::vector<Real>;

    /**
     * @brief The orthogonal wave functions are represented are vectors of Vals, as they are amplitudes for the non-orthogonal wave functions after the Gram-Schmidt procedure which is implemented in addBasisElem(). 
     * 
     */
    using OWF = std::vector<Val>;

    /**
     * @brief Handle for the user integer number type from Key.
     * 
     */
    using Int = typename Key::basetype;

    /**
     * @brief The collections of non-orthogonal wave functions.
     * 
     */
    std::vector<WF> waves = {};
    
    /**
     * @brief The collections of orthogonal wave functions,  build up using addBasisElem() and addPhoton().
     * 
     */
    std::vector<OWF> basis = {};        

    /**
     * @brief The function that returns the overlap given two WF. Can be set by the user using set()
     * 
     */
    Val (*get_ovlp)(const WF&, const WF&);      

    /**
     * @brief The tolerance for the amplitudes. Amplitudes with lower absolute value are deleted in clean().
     * 
     */
    Real tol = (Real) std::pow(10, -9);

    /**
     * @brief Loss is modelled as a map to a new Spatial&Polarization mode. This mode should be always higher than modes used for computation.
     * 
     */
    Int lossMode = 0;

    public:

        /**
         * @brief Construct a new (empty) State object.
         * 
         */
        State(){}

        /**
         * @brief Construct a new State object from a single key with amplitude 1.0.
         * 
         * @param k Key for the State (k:1.0). 
         */
        State(Key k){Par::insert(std::make_pair(k, (Val) 1.0));}

        /**
         * @brief Sets the overlap function to f.
         * 
         * @param f function that gets two vector<Val> and returns a Val, that is the overlap.
         */
        inline void set(Val (*f)(const WF&, const WF&)){get_ovlp = f;}

        /**
         * @brief Sets the tolerance to t.
         * 
         * @param t Real the tolerance is set to.
         */
        inline void set(Real t){tol = t;}

        /**
         * @brief Sets the lossMode to n.
         * 
         * @param n Int the lossMode is set to.
         */
        inline void set(Int n){lossMode = n;}

        /**
         * @brief Sets the actual State data to p.
         * 
         * @param p The data is set to. Moved input.
         */
        inline void set(Par&& p) noexcept {Par::operator=(p);}

        /**
         * @brief Sets the actual State data to p.
         * 
         * @param p The data is set to.
         */
        inline void set(const Par& p){Par::operator=(p);}

        /**
         * @brief Inserts/assigns a new Key Value pair.
         * 
         * @param k Key for new entry. If k is already in the State the amplitude is overwritten with v.
         * @param v new amplitude.
         */
        inline void set(const Key& k, const Val& v){Par::insert_or_assign(k, v);}

        /**
         * @brief Inserts/assigns a new Key Value pair.
         * 
         * @param k Key for new entry. If k is already in the State the amplitude is overwritten with v. Key is moved.
         * @param v new amplitude.
         */
        inline void set(Key&& k, const Val& v){Par::insert_or_assign(k, v);}

        /**
         * @brief Inserts/assigns a new Key Value pair. The key is implicitly given by (a,b,c).
         * 
         * @param a Spatial&Polarization mode
         * @param b Distinguishability mode
         * @param c Occupation number.
         * @param v new amplitude.
         */
        inline void set(Int a, Int b, Int c, const Val& v){Par::insert_or_assign(Key(a, b, c), v);} 

        /**
        * @brief Adds an state-like object to a State. If a Key is already in the state the amplitudes are added.
        * 
        * @param p State-like object to be added.
        */
        inline void add(const Par& p);

        /**
        * @brief Adds an state-like object to a State. If a Key is already in the state the amplitudes are added.
        * 
        * @param p State-like object to be added.
        */
        inline void add(Par&& p);

        /**
        * @brief Adds an scaled state-like object to a State. If a Key is already in the state the amplitudes are added.
        *  
        * @param p State-like object to be added.
        * @param v p is firsted scaled with v before added to this.
        */
        inline void add(const Par& p, const Val& v);

        /**
        * @brief Returns the norm of a state
        * 
        * @return Real The norm of the state
        */
        inline Real norm() const;

        /**
         * @brief Get the Par, i.e. boost::container::flat_map<Key, Val>,  object
         * 
         * @return Par 
         */
        inline Par getPar(){return (*this);}

        /**
         * @brief Get the Par, i.e. boost::container::flat_map<Key, Val>, object moved.
         * 
         * @return Par&& 
         */
        inline Par&& get_parMoved(){return std::move((*this));}

        /**
        * @brief Uses Gram-Schmidt procedure to add a basis element to the orthogonal basis. Used to get orthogonal Distinguishability modes.
        * 
        * @param wf Wave-function of the new entry represented as vector
        * @return std::vector<Val> Returns the inner products of wf with all basis elements after G.-S. procedure.
        */
        OWF addBasisElem(const WF&);

        /**
        * @brief Adds num new photon to the state with the wave function wf in the Spatial&Polarization mode m.
        * 
        * @param wf Wave function of the new photons represented as vector
        * @param mode Spatial&Polarization mode of the new photons
        * @param num Number of new photons
        */
        void addPhoton(const WF&, Int, Int);

        /**
        * @brief Applies a unitary U in second-quantization on the Spatial&Polarization modes in modes (max 2).
        * 
        * @param U Unitary represented as line
        * @param modes modes affected by the operation
        */
        inline void apply(const std::vector<Val>& U, const std::vector<Int>& modes);

        /**
        * @brief Applies a unitary U in second-quantization on the Spatial&Polarization modes in modes (max 1).
        *  
        * @param U Unitary, i.e. here a number with abs value 1
        * @param mode mode affected by the operation
        */
        inline void apply(const Val& U, const Int& mode);

        /**
        * @brief Swaps two Spatial&Polarization modes, used to implement a PBS.
        * 
        * @param a First Spatial&Polarization mode
        * @param b Second Spatial&Polarization mode
        */
        inline void swap(Int a, Int b);

        /**
         * @brief Removes all Key-Val pairs where the abs of the amplitude is lower than tol.
         * 
         */
        inline void clean();

        /**
        * @brief Normalizes the state.
        * 
        */
        inline void normalise();

        /**
        * @brief Multiplies the state with a scalar n
        * 
        * @param n Scalar to multiply the state with
        */
        inline void mul(Val);

        /**
        * @brief Computes the overlap to a measurement pattern and filters according to a fidelity reference.
        * 
        * Filters the state such that the measurement patter is overlapping and the fidelity reference is overlapping. 
        * Returns the part of the state that is accepted by the measurement and post-selection but is not overlapping with the fidelity refence.
        * 
        * @param ref The measurement pattern {S&P mode: total occupation-numer}. All S&P modes not in ref are accepted either way.
        * @param allModes vector of vector of S&P modes of fidelity reference used to pre-filter. At least for one of them all modes have to be occupied. This implements post-selections.
        * @param fref The fidelity reference. Used to further filter the measurement result, keeps that parts that could contribute to the fidelity.
        * @return State<Key, Val, Real> 
        */
        inline State<Key, Val, Real> overlapWithFilter(const boost::container::flat_map<Int, Int>&, const std::vector<std::vector<Int>>&, const std::vector<boost::container::flat_map<Int, Int>>&); 

        /**
        * @brief Keeps the part of the state that is either rejected because of the measurement result or because of post-selection.
        * 
        * @param MVec Measurement pattern {S&P mode: total occupation-numer}. Filters out the Keys that overlap with one of these.
        * @param allModes Modes for post-selection. Filters our those that are accepted by post-selection, i.e. for all of the vectors at least one is not occupied
        */
        inline void overlapCompl(const std::vector<boost::container::flat_map<Int,Int>>&, const std::vector<std::vector<Int>>& = {});

        /**
        * @brief Performes loss of one photons on the modes in modes
        *  
        * @param modes Spatial&Polarization modes to perform loss on
        */
        inline void loss(const std::vector<Int>&);

        /**
        * @brief Maps the current distinguishability conf to a different one - Only use for mapping to less distinguishable conf.
        * 
        * @param K Key that encodes the target distinguishability configuration.
        */
        inline void collapse(const Key&);

        /**
        * @brief Filters out these key that haven't at least one photon in on of the modes in allModes.
        * 
        * @param allModes Spatial&Distinguishability modes to check. 
        */
        inline void notEmpty(const std::vector<std::vector<Int>>&);

        /**
        * @brief Keeps the keys where all Distinguishability Modes are the same for one of the collection of modes in modes.
        * 
        * @param modes For one of the vectors in modes, all Distinguishability Modes have to be the same for the key to be included in the output.
        * @param facs Factors associated to the vectors in modes. These are used for the amplitudes correspoding to the keys, in which all the modes in the correspoding vector have the same Distinguishability Modes. 
        */
        inline void sameDModeDel(const std::vector<std::vector<Int>>&, const std::vector<Val>&);
};

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::add(const Par& p){
    bool nfound;
    typename Par::iterator rit;
    std::pair<typename Par::iterator, bool> pbi;
    for (typename Par::const_iterator it=p.cbegin(); it != p.cend(); it++){
        pbi = Par::emplace(std::make_pair(it->first, it->second));
        rit = pbi.first;
        nfound = pbi.second;
        if (!nfound){
            rit->second += it->second;
        }
    }
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::add(Par&& p){
    bool nfound;
    typename Par::iterator rit;
    std::pair<typename Par::iterator, bool> pbi;
    for (typename Par::iterator it=p.begin(); it != p.end(); it++){
        pbi = Par::emplace(std::make_pair(it->first, it->second));
        rit = pbi.first;
        nfound = pbi.second;
        if (!nfound){
            rit->second += it.second;
        }
    }
    p.clear();
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::add(const Par& p, const Val& v){
    bool nfound;
    typename Par::iterator rit;
    std::pair<typename Par::iterator, bool> pbi;
    for (typename Par::const_iterator it=p.cbegin(); it != p.cend(); it++){
        pbi = Par::emplace(it->first, v*(it->second));
        rit = pbi.first;
        nfound = pbi.second;
        if (!nfound){
            rit->second += (v * it->second);
        }
    }
}
 
template<class Key, class Val, class Real>
inline Real State<Key, Val, Real>::norm() const {
    Real r = 0;
    for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
        r += std::pow(std::abs(it->second), 2);
    }
    return std::sqrt(r);
}

template<class Key, class Val, class Real>
std::vector<Val> State<Key, Val, Real>::addBasisElem(const WF& wf){
    OWF bNew(basis.size(), 0.0);
    Val ov;
    OWF decomp;
    for (int i = 0; i<basis.size();i++){
        ov = -ovlpH<Val, Real>(basis[i], wf, get_ovlp, waves);
        decomp.push_back(-ov);
        for (int j = 0; j<basis[i].size();j++)
            bNew[j] += basis[i][j]*ov;
    }
    bNew.push_back((Val) 1.0);
    waves.push_back(wf);
    Val normC = ovlp<Val, Real>(bNew, bNew, get_ovlp, waves);
    Real norm = std::abs(normC);
    for (int i = 0;i<bNew.size();i++){
        bNew[i] = bNew[i] / std::sqrt(norm);
    }
    basis.push_back(bNew);
    decomp.push_back(ovlpH<Val, Real>(bNew, wf, get_ovlp, waves));
    Real norm2 = 0;
    for (Val a : decomp){
        norm2 += std::pow(std::abs(a), 2);
    }
    if (norm2!=0){
    for (int i = 0;i<bNew.size();i++){
        decomp[i] = decomp[i] / std::sqrt(norm2);
    }}
    return decomp;
}

template<class Key, class Val, class Real>
void State<Key, Val, Real>::addPhoton(const WF& wf, Int mode, Int num){
    Int s = Par::size();
    if (mode >= lossMode) lossMode = mode+1;
    if (s==0){
        set(mode, 0, num, (Val) (1.0));
        waves.push_back(wf);
        basis.push_back({(Val) 1.0});
        return;
    }
    Int index = -1;
    for (Int i = 0; i<waves.size();i++){
        if (std::abs(get_ovlp(wf, waves[i])) == (Real) 1.0) {index = i; break;}
    }
    if (index!= -1){
        Par SNew;
        std::pair<Key, Val> p=std::make_pair(Key(), (Val) 0.0);
        for (typename Par::const_iterator it=Par::cbegin(); it!=Par::cend(); it++){
            p = (*it);
            p.first.incr(mode, index, num);
            SNew.insert(SNew.cend(), std::move(p));
        }
        set(std::move(SNew));
        return;
    }
    Par SNew;
    OWF decomp = addBasisElem(wf);
    std::pair<Key, Val> p, p2;
    for (typename Par::const_iterator it=Par::cbegin(); it!=Par::cend(); it++){
        p = (*it);
        for (Int i = 0; i<decomp.size();i++){
            p2 = p;
            p2.first.addEnd(mode, i, num); 
            p2.second *= decomp[i];
            SNew.insert(SNew.cend(), std::move(p2));
        }
    }
    set(std::move(SNew));
    return;
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::apply(const std::vector<Val>& U, const std::vector<Int>& modes){
    Par S = get_parMoved(), S2;
    std::pair<Key, Val> p;
    for (typename Par::iterator it = S.begin(); it!=S.end(); it++){
        p = (*it);
        S2 = p.first.apply(U, modes, tol);
        add(S2, p.second);
    }
    clean();
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::apply(const Val& U, const Int& mode){
    std::pair<Key, Val> p;
    Val v;
    for (typename Par::iterator it = Par::begin(); it!=Par::end(); it++){
        p = (*it);
        v = p.second;
        it->second = v*p.first.apply(U, mode);
    }
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::swap(Int a, Int b){
    std::pair<Key, Val> p;
    Par par;
    for (typename Par::iterator it = Par::begin(); it!=Par::end(); it++){
        p = (*it);
        p.first.swap(a, b);
        par.insert(par.cend(), std::move(p));
        }
    set(std::move(par));
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::clean(){
    std::pair<Key, Val> p;
    Par par;
    for (typename Par::iterator it = Par::begin(); it!=Par::end(); it++){
        p = (*it);
        if (std::abs(p.second)>tol){
        par.insert(par.cend(), std::move(p));}
        }
    set(std::move(par));
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::normalise(){
    std::pair<Key, Val> p;
    Val n = norm();
    if (n == (Val) 0.0) n = 1.0;
    Par par;
    for (typename Par::iterator it = Par::begin(); it!=Par::end(); it++){
        p = (*it);
        p.second /= n;
        if (std::abs(p.second)>tol){
            par.insert(par.cend(), std::move(p));}
        }
    set(std::move(par));
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::mul(Val n){
    std::pair<Key, Val> p;
    Par par;
    for (typename Par::iterator it = Par::begin(); it!=Par::end(); it++){
        p = (*it);
        p.second *= n;
        if (std::abs(p.second)>tol){
            par.insert(par.cend(), std::move(p));}
        }
    set(std::move(par));
}
 
template<class Key, class Val, class Real>
State<Key, Val, Real> State<Key, Val, Real>::overlapWithFilter(const boost::container::flat_map<Int, Int>& ref, const std::vector<std::vector<Int>>& allModes, const std::vector<boost::container::flat_map<Int, Int>>& fref){
    Par p;
    State<Key, Val, Real> S;
    bool found;
    for (typename Par::iterator it=Par::begin(); it!=Par::end();it++){
        found = false;
        if (it->first.overlapping(ref) && it->first.notEmpty(allModes)){
            for (boost::container::flat_map<Int, Int> m : fref){
                if (it->first.overlapping(m)){
                    p.insert(p.cend(), (*it));
                    found = true;
                    break;
                }
            }
            if (!found){
                S.insert(S.cend(), (*it));}
            }
    }
    Par::operator=(p);
    return S;
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::overlapCompl(const std::vector<boost::container::flat_map<Int,Int>>& MVec, const std::vector<std::vector<Int>>& allModes){
    State<Key,Val,Real> p;
    bool hit = false;
    for (typename Par::iterator it=Par::begin(); it!=Par::end();it++){
        hit = false;
        for (boost::container::flat_map<Int,Int> ref : MVec){
            if (it->first.overlapping(ref)){hit = true;}
        }
        if (!hit || !(it->first.notEmpty(allModes)))
            p.insert(p.cend(), (*it));
    }
    Par::clear();
    set(std::move(p));
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::loss(const std::vector<Int>& modes){
    Par p=get_parMoved();
    Int maxLM=0;
    for (typename Par::iterator it = p.begin(); it!=p.end(); it++){
        add(it->first.template loss<Val>(modes, lossMode, maxLM), it->second);
    }
    if (maxLM>lossMode)
        lossMode = maxLM;
    clean();
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::collapse(const Key& K){
    boost::container::flat_map<Int, Int> f;
    Int j = 0;
    for (typename Key::const_iterator it=K.cbegin(); it!=K.cend(); it++){
        f[j] = it->first.second;
        j++;}
    Key K2;
    Par p;
    Val v;
    std::pair<typename Par::iterator, bool> pib;
    for (typename Par::iterator it = Par::begin(); it != Par::end(); it++){
        K2 = it->first;
        v = it->second;
        K2.collapse(f, v);
        pib = p.emplace(std::make_pair(K2, v));
        if (!pib.second)
            pib.first->second += v;
    }
    set(std::move(p));
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::notEmpty(const std::vector<std::vector<Int>>& allModes){
    Par p;
    for (typename Par::iterator it = Par::begin(); it != Par::end(); it++){
        if (it->first.notEmpty(allModes))
            p.insert(std::make_pair(it->first, it->second));
    }
    set(std::move(p));
    normalise();
}

template<class Key, class Val, class Real>
inline void State<Key, Val, Real>::sameDModeDel(const std::vector<std::vector<Int>>& modes, const std::vector<Val>& facs){
    Par p;
    int i;
    Key K;
    std::pair<typename Par::iterator, bool> pib;
    for (typename Par::iterator it = Par::begin(); it != Par::end(); it++){
        i = 0;
        for (std::vector<Int> m: modes){
            K = it->first;
            if (K.sameDModeDel(m)){
                pib = p.emplace(std::make_pair(K, it->second*facs[i]));
                if (!pib.second)
                    pib.first->second += it->second*facs[i];
                break;
            }
        i++;
        }
    }
    set(std::move(p));
}
#endif