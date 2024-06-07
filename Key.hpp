/**
 * @file Key.hpp
 * @author Fabian Wiesner (fabian.wiesner97@gmail.com)
 * @brief Defines the class Key used for the State sturcture.
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

#ifndef KEY_HPP
#define KEY_HPP

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
#include "KeyAux.hpp"

/**
 * @brief Key class, that implements the functions needed for State. Inherits from boost::container::flat_map<std::pair<Int, Int>, Int>.
 * 
 * @tparam Int Class used for integer numbers.
 */
template<class Int>
class Key : public boost::container::flat_map<std::pair<Int, Int>, Int>{
    public:

    /**
     * @brief Short handle for the parent type that provide the actual data structure.
     * 
     */
    using Par = boost::container::flat_map<std::pair<Int, Int>, Int>;

    /**
     * @brief Handle for the intereger number type to provide it for State.
     * 
     */
    using basetype = Int;

    /**
     * @brief Handle for State-like (actually State-parent) data structure.
     * 
     * @tparam Val Value-type, cf. State.
     */
    template<class Val>
    using SD = boost::container::flat_map<Key<Int>, Val>;

    /**
     * @brief Construct a new Key object.
     * 
     */
    Key(){Par();}

    /**
     * @brief Construct a new Key object from an entry.
     * 
     * @param a spatial & polarization mode
     * @param b distinguishability mode
     * @param c occupation number
     */
    Key(Int a, Int b, Int c){Par::insert_or_assign(std::make_pair(a, b), c);}

    /**
     * @brief Adds a boost::container::flat_map<std::pair<Int, Int>, Int>, if entry exist, occupation numbers are added.
     * 
     * @param p boost::container::flat_map<std::pair<Int, Int>, Int> to add.
     */
    inline void add(const Par& p){
        bool nfound;
        typename Par::iterator rit;
        std::pair<typename Par::iterator, bool> pbi;
        for (typename Par::const_iterator it= p.cbegin(); it != p.cend(); it++){
            pbi = Par::emplace(std::make_pair(it->first, it->second));
            rit = pbi.first;
            nfound = pbi.second;
            if (!nfound){
                rit->second += it->second;
            }
        }
    }

    /**
     * @brief Returns this as boost::container::flat_map<std::pair<Int, Int>, Int>.
     * 
     * @return Par, i.e. boost::container::flat_map<std::pair<Int, Int>, Int> 
     */
    operator Par() const {return (*this);}

    /**
     * @brief Removes entries with occupation number 0.
     * 
     */
    inline void clean(){
        Par p;
        for (typename Par::iterator it = Par::begin(); it!= Par::end(); it++){
            if (it->second>0)
                p.insert_or_assign(p.cend(), it->first, it->second);
        }
        Par::operator= (std::move(p));
    }

    /**
     * @brief Normalization factor for second-quantization for Spatial&Polarization modes a and b.
     * 
     * @tparam Val Value-type, cf. State
     * @param a Spatial&Polarization mode
     * @param b Spatial&Polarization mode
     * @return Val The normalization factor 
     */
    template<class Val>
    Val factor(Int a, Int b) const {
        Val r = 1.0;
        for (typename Par::const_iterator it = Par::cbegin(); it!= Par::cend(); it++){
            if (it->first.first== a || it->first.first== b) r*= (Val) facut(it->second);
        }
        return std::sqrt(r);
    }

    /**
     * @brief Normalization factor for second-quantization for all modes.
     * 
     * @tparam Val Value-type, cf. State
     * @return Val The normalization factor
     */
    template<class Val>
    Val factor() const {
        Val r = 1.0;
        for (typename Par::const_iterator it = Par::cbegin(); it!= Par::cend(); it++){
            r*= (Val) facut(it->second);
        }
        return std::sqrt(r);
    }

    /**
     * @brief Adds entry and the end.
     * 
     * @param a Spatial&Polarization mode
     * @param b Distinguishability mode
     * @param c Occupation number
     */
    void addEnd(const Int& a, const Int& b, const Int& c){Par::insert_or_assign(Par::cend(), std::make_pair(a, b), c);}

    /**
     * @brief Applies a Unitary in second quantization.
     * 
     * @tparam Val Amplitude type that is used in the State 
     * @tparam Real Real number type that should be used, e.g. float
     * @param U Unitary in line format
     * @param modes Affected modes, max 2.
     * @param tol Tolerance for the aplitude, amplitude with absolute value lower than tol are discarded
     * @return SD<Val> State similar data structure one gets, if one applies U on the implicit state (Key : 1.0).
     */
    template<class Val, class Real>
    inline SD<Val> apply(const std::vector<Val>& U, const std::vector<Int>& modes, const Real& tol) const {
        SD<Val> Ret, RetIter, RetIter2;
        Key K;
        Ret.insert(std::make_pair(K, (Val) 1.0));
        Int m, d, n, i;
        Int a = modes[0], b = modes[1];
        Val amp = 1.0;
        typename SD<Val>::const_iterator hint;
        std::pair<Key, Val> pr;
        std::pair<typename SD<Val>::iterator, bool> pib;
        for (typename Par::const_iterator it = Par::cbegin(); it!= Par::cend(); it++){
            m = it->first.first;
            d = it->first.second;
            n = it->second;
            if (m != a && m != b){
                for (typename SD<Val>::iterator it = Ret.begin(); it!= Ret.end(); it++){
                    K = it->first;
                    K.insert_or_assign(K.cend(), std::make_pair(m, d), n);
                    RetIter.insert_or_assign(RetIter.cend(), K, it->second);
                }
                K.clear();
                Ret = std::move(RetIter);
                continue;
            }
            i = (m == a) ? 0 : 1;
            
            for (Int j= 0; j<= n; j++){
                amp *= (Val) std::pow(U[i], j) * binomialCoeff<Val, Int>(n, j);
                amp *= (Val) std::pow(U[i+2], n-j)/((Val) std::sqrt(facut(n)));
                K = Key(modes[0], d, j);
                K.insert_or_assign(std::make_pair(modes[1], d), n-j);
                K.clean();
                RetIter.insert_or_assign(RetIter.cend(), K, amp);
                amp = (Val) (1.0);
            }
            hint = RetIter2.cend();
            for (typename SD<Val>::iterator riit = RetIter.begin(); riit != RetIter.end(); riit++){
                for (typename SD<Val>::iterator rit = Ret.begin(); rit != Ret.end(); rit++){
                    K = riit->first;
                    K.add(rit->first);
                    pr = std::make_pair(K, (riit->second) * (rit->second));
                    pib = RetIter2.emplace(pr);
                    if (!pib.second)
                        pib.first->second += pr.second;
                }
            }
            Ret = std::move(RetIter2);
            RetIter.clear();
        }
        for (typename SD<Val>::iterator it = Ret.begin(); it != Ret.end(); it++){
            if (std::abs(it->second) > tol)
                RetIter.insert_or_assign(RetIter.cend(), it->first, it->second * it->first.template factor<Val>(a, b));}
        return RetIter;
    }

    /**
     * @brief Applies a Unitary in second quantization for a single mode.
     * @tparam Val Amplitude type that is used in the State 
     * @param U Unitary in line format
     * @param mode affected mode
     * @return Val The factor for the amplitude
     */
    template<class Val>
    inline Val apply(const Val& U, const Int& mode) const {
        Int a = 0;
        std::pair<Int, Int> p;
        for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
            p = (it->first);
            if (p.first == mode){
                a += (it->second);
            }
        }
        return std::pow(U, a);
    }

    /**
     * @brief Swaps Spatial&Polarization modes a and b.
     * 
     * @param a Spatial&Polarization mode
     * @param b Spatial&Polarization mode
     */
    inline void swap(Int a, Int b){
        std::pair<std::pair<Int, Int>, Int> p;
        Par par;
        for (typename Par::iterator it= Par::begin(); it!= Par::end(); it++){
            p = (*it);
            if (p.first.first == a){
                p.first.first = b;
            }
            else if (p.first.first == b)
                p.first.first = a;
            par.insert(par.cend(), p);
        }
        Par::operator= (par);
    }

    /**
     * @brief Checks if the Key containes data of a boost::container::flat_map<Int, Int>.
     * 
     * @param ref The data to check if contained
     * @return Bool, if contained or not 
     */
    inline bool overlapping(const boost::container::flat_map<Int, Int>& ref) const {
        boost::container::flat_map<Int, Int> col;
        std::pair<typename boost::container::flat_map<Int, Int>::iterator, bool> pib;
        for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
            if (ref.find(it->first.first)!= ref.cend()){
                pib = col.emplace(std::make_pair(it->first.first, it->second));
                if (!pib.second){
                    pib.first->second += it->second;
                }
            }
        }
        typename boost::container::flat_map<Int, Int>::const_iterator it2;
        for (typename boost::container::flat_map<Int, Int>::const_iterator it = ref.cbegin(); it != ref.cend(); it++){
            it2 = col.find(it->first);
            if (it2 == col.cend()){
                if (it->second != 0) return false;
            }
            else
                if (it2->second != it->second)
                    return false;

                
        }
        return true;
    }
    

/*     inline bool overlapping(const Par& ref) const {
        typename Par::iterator it2;
        for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
            it2 = ref.find(it->first);
            if (it2!= ref.cend()){
                if (it2->second()!= it->second) return false;
            }
        }
        return true;
    } */

    /**
     * @brief Applies loss with uniform probability per photon on the modes.
     * 
     * @tparam Val Amplitude type that is used in the State 
     * @param modes Spatial&Polarization Modes affected by loss
     * @param lossMode Spatial&Polarization Mode the photon goes to
     * @param maxLM Saves the highest mode used for loss.
     * @return SD<Val>, i.e. boost::container::flat_map<Key<Int>, Val>, that is obtained by loss on (Key:1.0)
     */
    template<class Val>
    SD<Val> loss(const std::vector<Int>& modes, const Int& lossMode, Int& maxLM) const {
        SD<Val> R;
        Int total = 0;
        Key k = (*this), k2;
        Int j;
        Int s = modes.size();
        for (typename Par::iterator it = k.begin(); it != k.end(); it++){
            j = s - (modes.cend() - std::find(modes.cbegin(), modes.cend(), it->first.first));
            if (j != s && it->second>0){
                it->second -= 1;
                k2 = k;
                k2.insert_or_assign(k2.cend(), std::make_pair(lossMode+j, it->first.second), 1);
                R.insert_or_assign(R.cend(), k2, (Val) std::sqrt((it->second)+1));
                it->second += 1;
                total += it->second;
                if (lossMode+j>= maxLM) maxLM = lossMode + j + 1;
            }
        }
        if (total== 0){
            R.insert_or_assign(R.cend(), k, 1.0);
            return R;
        }
        Val n = std::sqrt(total);
        for (typename SD<Val>::iterator it = R.begin(); it!= R.end(); it++)
            it->second /= n;
        return R;
    }

    /**
     * @brief Maps the current distinguishability conf to a different one - Only use for mapping to less distinguishable conf.
     * 
     * @tparam Val Amplitude type that is used in the State 
     * @param f Encodes the mapping {old D.mode : new Dmode}
     * @param amp Amplitude for this key before the mapping
     */
    template<class Val>
    inline void collapse(const boost::container::flat_map<Int, Int>& f, Val& amp){
        Val pre = factor<Val>();
        Par p;
        Int a, b, n;
        std::pair<typename Par::iterator, bool> pib;
        for (typename Par::iterator it = Par::begin(); it!= Par::end(); it++){
            a = it->first.first;
            b = f.at(it->first.second);
            n = it->second;
            pib = p.emplace(std::make_pair(std::make_pair(a, b), n));
            if (!pib.second)
                pib.first->second += n;
        }
        Par::operator= (std::move(p));
        amp *= factor<Val>()/pre;
    }    

    /**
     * @brief Checks if the key has at least one photon in on of the modes in allModes.
     * 
     * @param allModes Spatial&Distinguishability modes to check
     * @return true, if there is at least one photon in one of the mode
     * @return false, otherwise
     */
    inline bool notEmpty(const std::vector<std::vector<Int>>& allModes) const {
        std::vector<Int> Msub;
        typename std::vector<Int>::iterator it2;
        for (std::vector<Int> M: allModes){
            Msub = M;
            for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
                it2 = std::find(Msub.begin(), Msub.end(), it->first.first);
                if (it2!= Msub.end() && (it->second)>0)
                    Msub.erase(it2);
                if (Msub.size()== 0)
                    return true;
            }
        }
        return allModes.size()== 0;
    }

    /**
     * @brief Checks if all S&P modes in m have the same distinguishability mode and deletes all modes in m.
     * 
     * @param m Spatial&Polarization modes to check
     * @return true, if all modes in m have the same distinguishability mode
     * @return false, otherwise
     */
    inline bool sameDModeDel(const std::vector<Int>& m){
        typename std::vector<Int>::const_iterator it2;
        Par p;
        Int dmode = -1;
        int founds = 0;
        for (typename Par::const_iterator it = Par::cbegin(); it != Par::cend(); it++){
            it2 = std::find(m.cbegin(), m.cend(), it->first.first);
            if (it2!= m.cend() && it->second!= 0){
                if (dmode == -1) dmode = it->first.second;
                else if (it->first.second != dmode) return false;
                founds++;
            }
            else if(it->second!= 0)
                p.insert(p.cend(), (*it));
        }
        Par::operator= (std::move(p));
        return founds == m.size();
    }
};


#endif