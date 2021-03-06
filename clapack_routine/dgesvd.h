#ifndef __CLAPACK_DGESVD_H
#define __CLAPACK_DGESVD_H

#include "f2c.h"
#include "blaswrap.h"

#ifdef __cplusplus 	
extern "C" {	
#endif	


/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *s, doublereal *u,
	integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, 
        integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif


#endif /* __CLAPACK_DGESVD_H */
