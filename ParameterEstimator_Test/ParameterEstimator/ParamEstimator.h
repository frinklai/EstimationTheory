/*
 *	ParamEstimator.h
 *
 *	Description:
 *		Basic class of parameter estimator
 *
 *
 *		
 * 	History:
 *	 	Author			Date			Modify Reason		
 *		----------------------------------------------------------------
 *		Chi-Yi Tsai		2014/01/25		File Creation    
 *
 */
 
#ifndef PARAM_ESTIMATOR_ENUM_METHODS
#define PARAM_ESTIMATOR_ENUM_METHODS

#include "..\Matrix\fMatrix.h"
#include "..\MatchingPoints\MatchingPoints.h"
#define Mat1PATH "../Debug/Mat1.txt"
#define Vec1PATH "../Debug/Vec1.txt"
//#include <math.h>
enum ParamEstiMethod {LS = 1, WLS, ML};

typedef struct	st_LS_Param
{
	fMatrix*	pMat_H;
	fVector*	pVec_Z;
	// Optional	(for error-variance computing)
	fMatrix*	pMat_Vz;
} LS_Param;

typedef struct	st_WLS_Param
{
	fMatrix*	pMat_H;
	fMatrix*	pMat_W;
	fVector*	pVec_Z;
	// Optional	(for error-variance computing)
	fMatrix*	pMat_Vz;
} WLS_Param;

typedef struct	st_ML_Param
{
	fMatrix*	pMat_H;
	fVector*	pVec_Z;
	// Necessary for state estimation and error-variance computing
	fMatrix*	pMat_Vz;
} ML_Param;

#endif // PARAM_ESTIMATOR_ENUM_METHODS

#ifndef __PARAM_ESTIMATOR_INCLUDED__
#define __PARAM_ESTIMATOR_INCLUDED__

class CParamEstimator
{
public:
	CParamEstimator();
	~CParamEstimator();

	fMatrix*	SolveOptParam(fVector*	pfVecOptParam);

	void		SetParamEstiMethod(ParamEstiMethod Method);
	void		SetMethodParameters(ParamEstiMethod Method, void*	pParam);
	
	ParamEstiMethod	GetParamEstiMethod(void) const;
	void*		GetMethodParameters(ParamEstiMethod Method) const;

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  Add Variables or Functions, by Yu-Cheng, Lai                           *
 *                                                                         *
 *-------------------------------------------------------------------------*/
	// For Estimator
	fMatrix *Vv,   *Var_ls, *W, *Var_wls, *Var_ml;
	fVector *P_ls, *P_wls,  *v, *P_ml;	

	LS_Param  *	LS1 ;
	WLS_Param *	WLS1;
	ML_Param  * ML1 ;

	ParamEstiMethod SolveMethod;

	fMatrix* CParamEstimator::LeastSquare(fVector* P_ls);
	fMatrix* Weighted_LS(fVector* P_wls);
	fMatrix* Maximum_Likelihood(fVector* P_ml);
	void LMV();

	// For InitiHomographyEstimation
	fMatrix* Mat_Ref_X, *Mat_X3; 
	fMatrix* Mat_Tr  ,	*Mat_Ts;

	fVector* Vec_Xr, *Vec_Yr;
	fVector* Vec_Xs, *Vec_Ys;
	fVector* Vec_X1, *Vec_Y1;
	fVector* Vec_X2, *Vec_Y2;

	double Xr_Mean;
	double Yr_Mean;
	double Xs_Mean;
	double Ys_Mean;
	double Sr, Ss;
	
	fMatrix* B;
	fMatrix* InitiHomographyEstimation(double &err);
	fMatrix Generate_tmpB(double* Data);
	fVector Generate_tmpY(double* Data);
};

#endif // __PARAM_ESTIMATOR_INCLUDED__