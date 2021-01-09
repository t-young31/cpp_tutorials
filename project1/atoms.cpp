//
// Created by tom on 01/01/2021.
//

#include "atoms.h"
#include <vector>

double atomic_weight(int atomic_number){
    const std::vector<double> atomic_weights = {1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067,
                                                15.9994, 18.9984032, 2.01797, 22.989770, 24.3050, 26.981538,
                                                28.0855, 30.973761, 32.065, 35.453, 39.948, 39.0983, 40.078,
                                                44.955910, 47.867, 50.9415, 51.9961, 54.938049, 55.845,
                                                58.933200, 58.6934, 63.546, 65.409, 69.723, 72.64};
    // Atomic numbers start from H = 1, so decrement 1
    return atomic_weights[atomic_number-1];
}
