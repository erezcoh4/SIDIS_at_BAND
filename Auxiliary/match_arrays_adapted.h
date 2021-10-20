#ifndef __MATCH_ARRAYS_H__
#define __MATCH_ARRAYS_H__


#include <algorithm>
#include <iterator>
#include <utility>


class match_arrays {
    
    // based on [https://stackoverflow.com/questions/36669297/most-efficient-way-to-find-index-of-matching-values-in-two-sorted-arrays-using-c#comment60933333_36669297]
    // adapted to work on the ifarm
    
public:
    match_arrays(){};
    ~match_arrays(){};
    
    // helper structure for the search
    template<class Range, class Out>
    struct search_data {
        // is any there clearer way to get iterator that might be either
        // a Range::const_iterator or const T*?
        using iterator = decltype(std::cbegin(std::declval<Range&>()));
        iterator curr;
        const iterator begin, end;
        Out out;
    };
    
    template<class Range, class Out>
    auto init_search_data(const Range& range, Out out) {
        return search_data<Range, Out>{
            std::begin(range),
            std::begin(range),
            std::end(range),
            out,
        };
    }
    
    template<class Range, class Out1, class Out2>
    void match_indices(const Range& in1, const Range& in2, Out1 out1, Out2 out2) {
        auto search_data1 = init_search_data(in1, out1);
        auto search_data2 = init_search_data(in2, out2);
        
        // initial order is arbitrary
        auto lesser = &search_data1;
        auto greater = &search_data2;
        
        // if either range is exhausted, we are finished
        int Nmatched=0;
        while(lesser->curr != lesser->end
              && greater->curr != greater->end ) {
            // difference of first values in each range
            auto delta = *greater->curr - *lesser->curr;
            
            if(!delta) { // matching value was found
                // store both results and increment the iterators
                *lesser->out++ = std::distance(lesser->begin, lesser->curr++);
                *greater->out++ = std::distance(greater->begin, greater->curr++);
                continue; // then start a new iteraton
            }
            
            if(delta < 0) { // set the order of ranges by their first value
                std::swap(lesser, greater);
                delta = -delta; // delta is always positive after this
            }
            
            // next crossing cannot be farther than the delta
            // this assumption has following pre-requisites:
            // range is sorted, values are integers, values in the range are unique
            auto range_left = std::distance(lesser->curr, lesser->end);
            auto upper_limit =
            std::min(range_left, static_cast<decltype(range_left)>(delta));
            
            // exponential search for a sub range where the value at upper bound
            // is greater than target, and value at lower bound is lesser
            auto target = *greater->curr;
            auto lower = lesser->curr;
            auto upper = std::next(lower, upper_limit);
            for(int i = 1; i < upper_limit; i *= 2) {
                auto guess = std::next(lower, i);
                if(*guess >= target) {
                    upper = guess;
                    break;
                }
                lower = guess;
            }
            
            // skip all values in lesser,
            // that are less than the least value in greater
            lesser->curr = std::lower_bound(lower, upper, target);
        }
    }
    
protected:
    
};
#endif

