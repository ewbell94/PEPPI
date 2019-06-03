#define IPERGRP (4)

#define WINL (-7)
#define WINR (7)

#define NUM_IN	((WINR-WINL+1)*IPERGRP+3+1)	/* number of input units */
#define NUM_HID (55)			/* number of hidden units */
#define NUM_OUT (3) 			/* number of output units */

#define TOTAL		(NUM_IN + NUM_HID + NUM_OUT)

#define OFF_VAL	        (0.01)		/* Target OFF value */
#define ON_VAL	        (0.99)		/* Target ON value */

#define IALPHA		(0.9)		/* Initial smoothing factor */
#define ILRATE		(0.001)		/* Initial learning rate */
#define RWMAX		(0.03)		/* Max random weight value */
#define SAMPINTV	(1)		/* Sampling interval */
#define SAMPSIZE	(0)		/* Samples per weight update */

#define CPU_MAX		(10*3600.0)	/* Maximum CPU time (secs) for learning */

#define noEPOCHMODE
