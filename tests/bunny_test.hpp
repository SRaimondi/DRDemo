//
// Created by Simon on 10.09.2017.
//

#ifndef DRDEMO_BUNNY_TEST_HPP
#define DRDEMO_BUNNY_TEST_HPP

namespace drdemo {

    /**
     * Target bunny render scene
     */
    void BunnyTest(int start_resolution, float res_multiplier, int ref_steps);

    /**
     * Bunny test starting from smooth approximation
     */
    void BunnyTestSmooth(float res_multiplier, int ref_steps);

} // drdemo namespace

#endif //DRDEMO_BUNNY_TEST_HPP
