/**
 *
 *
 * Test suite for aaflib to compare with Julia
 */

#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <string>
#include "aa_aaf.h"

using namespace std;
namespace utf = boost::unit_test;

/**
 * Requires Boost.Test. If using Linux then install using your package manager.
 * To compile (and run) use these commands:

g++ -ggdb -frounding-math -DMAXORDER=40 -I. -I${HOME}/lib/filib-3.0.2/include \
    -I/usr/include -I$(pwd)/aaflib-0.1 -fpermissive -std=c++11 -c testsuite.cpp

g++ -L/usr/lib -L$(pwd)/aaflib-0.1 -L${HOME}/lib/filib-3.0.2/lib -o testsuite \
    testsuite.o -laaf -lprim -lgsl -llapack -lblas -lcblas -lstdc++ \
    -lboost_unit_test_framework

./testsuite --log_level=test_suite

 * Be careful when initializing tolerance
 *
 * TODO: use a proper Makefile
 */
BOOST_AUTO_TEST_SUITE(test_aaf, * utf::tolerance(0.00001))

    BOOST_AUTO_TEST_CASE( test_inverse ) {
        unsigned len = 5;
        // AAF 1
        double center1 = 26.10;
        double dev1[] {2.11, -3.03, 4.59, 1.0, -10.};
        unsigned ind1[] {1,3,5,8,10};
        AAF a1(center1, dev1, ind1, len);
        AAF o1 = inv(a1);
        o1.aafprint();
        cout << "coefficients:" << endl;
        cout.precision(17);
        for(int ii = 0; ii < o1.getlength()+1; ii++)
            cout << o1[ii] << endl;
        cout << "-------------" << endl;
    }

    BOOST_AUTO_TEST_CASE( test_inverse_2 ) {
        unsigned len = 4;
        // AAF 1
        double center1 = 4.0;
        double dev1[] {-3.33, 9.0, -1.5, 5.25};
        unsigned ind1[] {2, 3, 5, 6};
        AAF a1(center1, dev1, ind1, len);
        AAF o1 = inv(a1);
        o1.aafprint();
        cout << "coefficients:" << endl;
        cout.precision(17);
        for(int ii = 0; ii < o1.getlength()+1; ii++)
            cout << o1[ii] << endl;
        cout << "-------------" << endl;
    }

    BOOST_AUTO_TEST_CASE( test_power ) {
        unsigned len = 5;
        // AAF 1
        double center1 = 26.10;
        double dev1[] {2.11, -3.03, 4.59, 1.0, -10.};
        unsigned ind1[] {1,3,5,8,10};
        AAF a1(center1, dev1, ind1, len);
        AAF o1 = a1^2;
        o1.aafprint();
        cout << "coefficients:" << endl;
        cout.precision(17);
        for(int ii = 0; ii < o1.getlength()+1; ii++)
            cout << o1[ii] << endl;
        cout << "-------------" << endl;
    }

    BOOST_AUTO_TEST_CASE( test_generic ) {
        unsigned len = 5;
        // AAF 1
        double center1 = 12.10;
        double dev1[] {0.11, -1.03, 0.59, 0.02, -0.42};
        unsigned ind1[] {1,3,5,8,10};
        AAF a1(center1, dev1, ind1, len);
        // AAF 2
        double center2 = -6.0;
        double dev2[] {-0.22, -2.1, 1.0, 0.03, -0.6};
        unsigned ind2[] {2,3,8,9,10};
        AAF a2(center2, dev2, ind2, len);
        BOOST_CHECK_EQUAL( a1.getcenter(), center1 );
        BOOST_CHECK_EQUAL( a1.rad(), 0.11 + 1.03 + 0.59 + 0.02 + 0.42 );
        BOOST_CHECK_EQUAL( a2.rad(), 0.22 + 2.1 + 1.0 + 0.03 + 0.6 );
        AAF a_addn_cst = a1 + 3.0;
        AAF a_subt_cst = a1 - 2.1;
        AAF a_mult_cst = a1 * 2.0;
        AAF a_divi_cst = a1 / 2.0;
        AAF a_addn = a1 + a2;
        AAF a_subt = a1 - a2;
        AAF a_mult = a1 * a2;
        AAF a_divi = a1 / a2;
        unsigned L = 10;
        AAF a_XXs[] = { a1, a2, a_addn_cst, a_subt_cst, a_mult_cst,
                        a_divi_cst, a_addn, a_subt, a_mult, a_divi };
        string desc[] = {
                "a1", "a2", "const addition", "const subtraction",
                "const multiplication", "const division", "addition",
                "substraction", "multiplication", "division"
        };

        for(int i = 0; i < L; i++) {
            cout << desc[i] << endl;
            a_XXs[i].aafprint();
            cout << endl;
        }
    }

BOOST_AUTO_TEST_SUITE_END()












