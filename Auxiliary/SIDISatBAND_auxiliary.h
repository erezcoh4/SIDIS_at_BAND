// Erez O. C., Oct-6, 2021
#ifndef __DCFID_SIDIS_H__
#define __DCFID_SIDIS_H__


class SIDISatBAND_auxiliary {
public:
    SIDISatBAND_auxiliary();
    ~SIDISatBAND_auxiliary();

    Chi2PID_pion_upperBound( Double_t p, Double_t C);
    Chi2PID_pion_lowerBound( Double_t p, Double_t C);
    
protected:

};

#endif
