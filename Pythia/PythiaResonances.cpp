// Author: Claude Pruneau   09/25/2019

/***********************************************************************
 * Copyright (C) 2019, Claude Pruneau.
 * All rights reserved.
 * Based on the ROOT package and environment
 *
 * For the licensing terms see LICENSE.
 **********************************************************************/

/* this file stores the resonances of interest for suppressing from
 * the particle stack when so required by the analyses               */

const int nNoOfResonanceIndexes = 5000;
int resonanceIndexes[nNoOfResonanceIndexes] = {0};    /* a value of 1 will mean suppress the resonance and their decay products */
int toConsiderResonances[nNoOfResonanceIndexes] = {}; /* the list to potentially considered resonances */

void initializeResonanceIndexes()
{
  resonanceIndexes[113] = 0; // \rho(770)^{0}
  resonanceIndexes[213] = 0; // \rho(770)^{+}
  /*resonanceIndexes[-213] = 0; // \rho(770)^{-} */
  resonanceIndexes[221] = 0;  // \eta
  resonanceIndexes[223] = 0;  // \omega(782)
  resonanceIndexes[331] = 0;  // \eta^{\prime}(958)
  resonanceIndexes[333] = 0;  // \phi(1020)
  resonanceIndexes[1114] = 0; // \Delta(1232)^{-}
  /*resonanceIndexes[-1114] = 0; \bar{\Delta}(1232)^{+} */
  resonanceIndexes[2114] = 0; // \Delta(1232)^{0}
  /*resonanceIndexes[-2114] = 0; \bar{\Delta}(1232)^{0} */
  resonanceIndexes[2214] = 0; // \Delta(1232)^{+}
  /*resonanceIndexes[-2214] = 0; \bar{\Delta}(1232)^{-} */
  resonanceIndexes[2224] = 0; // \Delta(1232)^{++}
  /*resonanceIndexes[-2224] = 0; \bar{\Delta}(1232)^{--} */
  resonanceIndexes[3112] = 0; // \Sigma^{-}
  /*resonanceIndexes[-3112] = 0; \bar{\Sigma}^{+} */
  resonanceIndexes[3114] = 0; // \Sigma(1385)^{-}
  /*resonanceIndexes[-3114] = 0; \bar{\Sigma}(1385)^{+} */
  resonanceIndexes[3122] = 0; // \Lambda
  /*resonanceIndexes[-3122] = 0; // \bar{\Lambda} */
  resonanceIndexes[3212] = 0; // \Sigma^{0}
  /*resonanceIndexes[-3212] = 0; \bar{\Sigma}^{0} */
  resonanceIndexes[3214] = 0; // \Sigma(1385)^{0}
  /*resonanceIndexes[-3214] = 0; \bar{\Sigma}(1385)^{0} */
  resonanceIndexes[3222] = 0; // \Sigma^{+}
  /*resonanceIndexes[-3222] = 0; \bar{\Sigma}^{-} */
  resonanceIndexes[3224] = 0; // \Sigma(1385)^{+}
  /*resonanceIndexes[-3224] = 0; \bar{\Sigma}(1385)^{-} */
  resonanceIndexes[3312] = 0; // \Xi^{-}
  /*resonanceIndexes[-3312] = 0; \bar{\Xi}^{+} */
  resonanceIndexes[3314] = 0; // \Xi(1530)^{-}
  /*resonanceIndexes[-3314] = 0; \bar{\Xi}(1530)^{+} */
  resonanceIndexes[3322] = 0; // \Xi^{0}
  /*resonanceIndexes[-3322] = 0; \bar{\Xi}^{0} */
  resonanceIndexes[3324] = 0; // \Xi(1530)^{0}
  resonanceIndexes[3334] = 0; // \Omega^{-}
  /*resonanceIndexes[-3334] = 0; // \bar{\Omega}^{+} */
}

void initializeToConsiderResonances()
{
  toConsiderResonances[113] = 1; // \rho(770)^{0}
  toConsiderResonances[213] = 1; // \rho(770) ^ {+}
  /*toConsiderResonances[-213] = 1; // \rho(770)^{-} */
  toConsiderResonances[221] = 1;  // \eta
  toConsiderResonances[223] = 1;  // \omega(782)
  toConsiderResonances[331] = 1;  // \eta^{\prime}(958)
  toConsiderResonances[333] = 1;  // \phi(1020)
  toConsiderResonances[1114] = 1; // \Delta(1232)^{-}
  /*toConsiderResonances[-1114] = 1; \bar{\Delta}(1232)^{+} */
  toConsiderResonances[2114] = 1; // \Delta(1232)^{0}
  /*toConsiderResonances[-2114] = 1; \bar{\Delta}(1232)^{0} */
  toConsiderResonances[2214] = 1; // \Delta(1232)^{+}
  /*toConsiderResonances[-2214] = 1; \bar{\Delta}(1232)^{-} */
  toConsiderResonances[2224] = 1; // \Delta(1232)^{++}
  /*toConsiderResonances[-2224] = 1; \bar{\Delta}(1232)^{--} */
  toConsiderResonances[3112] = 1; // \Sigma^{-}
  /*toConsiderResonances[-3112] = 1; \bar{\Sigma}^{+} */
  toConsiderResonances[3114] = 1; // \Sigma(1385)^{-}
  /*toConsiderResonances[-3114] = 1; \bar{\Sigma}(1385)^{+} */
  toConsiderResonances[3122] = 1; // \Lambda
  /*toConsiderResonances[-3122] = 1; // \bar{\Lambda} */
  toConsiderResonances[3212] = 1; // \Sigma^{0}
  /*toConsiderResonances[-3212] = 1; \bar{\Sigma}^{0} */
  toConsiderResonances[3214] = 1; // \Sigma(1385)^{0}
  /*toConsiderResonances[-3214] = 1; \bar{\Sigma}(1385)^{0} */
  toConsiderResonances[3222] = 1; // \Sigma^{+}
  /*toConsiderResonances[-3222] = 1; \bar{\Sigma}^{-} */
  toConsiderResonances[3224] = 1; // \Sigma(1385)^{+}
  /*toConsiderResonances[-3224] = 1; \bar{\Sigma}(1385)^{-} */
  toConsiderResonances[3312] = 1; // \Xi^{-}
  /*toConsiderResonances[-3312] = 1; \bar{\Xi}^{+} */
  toConsiderResonances[3314] = 1; // \Xi(1530)^{-}
  /*toConsiderResonances[-3314] = 1; \bar{\Xi}(1530)^{+} */
  toConsiderResonances[3322] = 1; // \Xi^{0}
  /*toConsiderResonances[-3322] = 1; \bar{\Xi}^{0} */
  toConsiderResonances[3324] = 1; // \Xi(1530)^{0}
  toConsiderResonances[3334] = 1; // \Omega^{-}
  /*toConsiderResonances[-3334] = 1; // \bar{\Omega}^{+} */
}
