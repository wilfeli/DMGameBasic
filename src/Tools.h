/*
 * Tools.h
 *
 *  Created on: Apr 22, 2014
 *      Author: wilfeli
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <Eigen/Dense>

namespace Tools{

extern void mgrid(Eigen::MatrixXd, Eigen::MatrixXd*);
extern void mgrid_test(Eigen::MatrixXd, Eigen::MatrixXd*);
extern void print_vector(std::vector<double>);
extern void print_vector(std::vector<std::vector<double>>);
    
class MyRNG{
public:
    MyRNG(double);
    
    double GetUniform();
    uint GetUint();
    
    double state = 2013;
    bool hasSpare = false;
    double rn1 = 0.0;
    double rn2 = 0.0;
    
    uint m_w = 521288629;
    uint m_z = 362436069;
};

extern double get_normal(double, double, MyRNG&);
extern int get_int(int, int, MyRNG&);

template <class _RandomAccessIterator>
void
random_shuffle(_RandomAccessIterator __first, _RandomAccessIterator __last, MyRNG& rng){
    typedef typename std::iterator_traits<_RandomAccessIterator>::difference_type difference_type;
    difference_type __d = __last - __first;
    if (__d > 1){
        for (--__last, --__d; __first < __last; ++__first, --__d){
            //get random number from 0 to __d
            difference_type __i = get_int(0, __d, rng);
            if (__i != difference_type(0))
                swap(*__first, *(__first + __i));
        };
    };
};

};


#endif /* TOOLS_H_ */
