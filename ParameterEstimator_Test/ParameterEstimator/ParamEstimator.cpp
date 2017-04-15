#include "stdafx.h"
#include "ParamEstimator.h"

CParamEstimator::CParamEstimator()
{
	//Dynamic Allocate LS Parameter Memory 
	LS1 = new LS_Param;
	LS1->pMat_H = new fMatrix;
	LS1->pVec_Z = new fVector;
	P_ls = new fVector;
	Var_ls = new fMatrix;

	//Dynamic Allocate WLS Parameter Memory 
	WLS1 = new WLS_Param;
	WLS1->pMat_H = new fMatrix;
	WLS1->pVec_Z = new fVector;
	P_wls = new fVector;
	Var_wls = new fMatrix;

	//Dynamic Allocate ML Parameter Memory 
	ML1 = new ML_Param;
	ML1->pMat_H = new fMatrix;
	ML1->pVec_Z = new fVector;
	Var_ml  = new fMatrix;
	P_ml	= new fVector;
	
}

CParamEstimator::~CParamEstimator()
{
	delete LS1->pMat_H;
	delete LS1->pVec_Z;
	delete LS1;

	delete WLS1->pMat_H;
	delete WLS1->pVec_Z;
	delete WLS1;

	delete ML1->pMat_H;
	delete ML1->pVec_Z;
	delete ML1;
}

fMatrix* CParamEstimator::SolveOptParam(fVector*	pfVecOptParam)
{
	switch(SolveMethod)
	{
	case LS:
		return LeastSquare(pfVecOptParam);
		break;

	case WLS:
		return Weighted_LS(pfVecOptParam);
		break;

	case ML:
		return Maximum_Likelihood(pfVecOptParam);
		break;

	default:
		cout<<"ERROR TYPE!!\n";
		break;
	}
}

void CParamEstimator::SetParamEstiMethod(ParamEstiMethod Method)
{
	this->SolveMethod = Method;
}

void CParamEstimator::SetMethodParameters(ParamEstiMethod Method, void*	pParam)
{
	switch(Method)
	{
		case LS:
			this->LS1 = (LS_Param*)pParam;
			SetParamEstiMethod(LS);
			break;

		case WLS:
			this->WLS1 = (WLS_Param*)pParam;
			SetParamEstiMethod(WLS);
			break;

		case ML:
			this->ML1 = (ML_Param*)pParam;
			SetParamEstiMethod(ML);
			break;

		default:
			cout<<"ERROR TYPE!!\n";
			break;
	}
}

ParamEstiMethod	CParamEstimator::GetParamEstiMethod(void) const
{
	return this->SolveMethod;
}

void* CParamEstimator::GetMethodParameters(ParamEstiMethod Method) const
{
	switch(Method)
	{
	case LS:
		return (void*)this->LS1;

	case WLS:
		return (void*)this->WLS1;

	case ML:
		return (void*)this->ML1;

	default:
		cout<<"ERROR TYPE!!\n";
		break;
	}
}

fMatrix* CParamEstimator::LeastSquare(fVector* P_ls)
{
	cout<<"\n=========================LS=============================\n";
	
	//=========================NEW=========================
	v	 = new fVector;
	Vv	 = new fMatrix;
	//=========================NEW=========================

	*P_ls = (Inverse(Transp(*LS1->pMat_H) * *LS1->pMat_H) * Transp(*(LS1->pMat_H)) * *LS1->pVec_Z);

	*v    = (*LS1->pVec_Z - *LS1->pMat_H * *P_ls);
	*v	  = (*v - Mean(*v));
	
	*Vv   = (Outer(*v, *v)/(v->GetSize()-1)) ;
	*Vv	  =  Diag(Diag(*Vv));
	*Var_ls = (Inverse( Transp(*LS1->pMat_H) * *LS1->pMat_H ) * Transp(*LS1->pMat_H) 
			 * *Vv * *LS1->pMat_H * Inverse( Transp(*LS1->pMat_H) * *LS1->pMat_H ));

	//=========================delete=========================
	delete v;
	delete Vv;
	//=========================delete=========================
	return Var_ls;

}

fMatrix* CParamEstimator::Weighted_LS(fVector* P_wls)
{
	cout<<"\n======================Weighted_LS======================\n";
	//=========================NEW=========================
	v	 = new fVector;
	Vv	 = new fMatrix;
	W	 = new fMatrix;
	//=========================NEW=========================
	
	//----------------------------------LS------------------------------------------------------
	*P_ls = (Inverse(Transp(*WLS1->pMat_H) * *WLS1->pMat_H) * Transp(*(WLS1->pMat_H)) * *WLS1->pVec_Z);

	*v    = (*WLS1->pVec_Z - *WLS1->pMat_H * *P_ls);
	*v	  = (*v - Mean(*v));

	*Vv   = (Outer(*v, *v)/(v->GetSize()-1)) ;
	*Vv	  =  Diag(Diag(*Vv));
	*Var_ls = (Inverse( Transp(*WLS1->pMat_H) * *WLS1->pMat_H ) * Transp(*WLS1->pMat_H) 
		* *Vv * *WLS1->pMat_H * Inverse( Transp(*WLS1->pMat_H) * *WLS1->pMat_H ));
	//-------------------------------------------------------------------------------------------------

	
	*v = *WLS1->pVec_Z - *WLS1->pMat_H * *this->P_ls;
	*W = Diag( 1/(1+DotMul( *v, *v )) );
	*P_wls   = Inverse( Transp(*WLS1->pMat_H)   * *W * *WLS1->pMat_H ) * Transp(*WLS1->pMat_H) * *W * *WLS1->pVec_Z;

	*Var_wls = Inverse( Transp(*WLS1->pMat_H) * *W * *WLS1->pMat_H ) * Transp(*WLS1->pMat_H) * *W 
			   * *Vv * *W * *WLS1->pMat_H * Inverse( Transp(*WLS1->pMat_H)* *W * *WLS1->pMat_H );
	
	//=========================delete=========================
	delete v;
	delete Vv;
	delete W;
	//=========================delete=========================
	
	return Var_wls;
}

fMatrix* CParamEstimator::Maximum_Likelihood(fVector* P_ml)
{
	cout<<"\n===================Maximum_Likelihood===================\n";
	//=========================NEW=========================
	v	 = new fVector;
	Vv	 = new fMatrix;
	W = new fMatrix;
	//=========================NEW=========================

	//----------------------------------LS------------------------------------------------------
	*P_ls = (Inverse(Transp(*ML1->pMat_H) * *ML1->pMat_H) * Transp(*(ML1->pMat_H)) * *ML1->pVec_Z);

	*v    = (*ML1->pVec_Z - *ML1->pMat_H * *P_ls);
	*v	  = (*v - Mean(*v));

	*Vv   = (Outer(*v, *v)/(v->GetSize()-1)) ;
	*Vv	  =  Diag(Diag(*Vv));
	*Var_ls = (Inverse( Transp(*ML1->pMat_H) * *ML1->pMat_H ) * Transp(*ML1->pMat_H) 
		* *Vv * *ML1->pMat_H * Inverse( Transp(*ML1->pMat_H) * *ML1->pMat_H ));
	//-------------------------------------------------------------------------------------------------

	*W = Inverse(*Vv);
	*P_ml   = Inverse( Transp(*ML1->pMat_H) * *W * *ML1->pMat_H) * Transp(*ML1->pMat_H) * *W * *ML1->pVec_Z;
	*Var_ml = Inverse( Transp(*ML1->pMat_H) * Inverse(*Vv) * *ML1->pMat_H );

	//=========================delete=========================
	delete v;
	delete Vv;
	delete W;
	//=========================delete=========================

	return Var_ml;
}

fMatrix* CParamEstimator::InitiHomographyEstimation(double &Out_err)
{
	//=========================NEW=========================
	fVector* Y = new fVector;
	Mat_Ref_X  = new fMatrix(g_Ref_X, g_nNumPoints, 2);
	Mat_X3     = new fMatrix(g_x3,	 g_nNumPoints, 2);

	Vec_Xr	   = new fVector;		Vec_Yr	  = new fVector;	
	Vec_Xs	   = new fVector;		Vec_Ys	  = new fVector;
	Vec_X1	   = new fVector;		Vec_Y1	  = new fVector;
	Vec_X2	   = new fVector;		Vec_Y2	  = new fVector;
	//=========================NEW=========================

	//=================================Step1================================
	*Vec_Xr    = Mat_Ref_X->GetCol(0);
	*Vec_Yr    = Mat_Ref_X->GetCol(1);

	Xr_Mean = Mean(*Vec_Xr);
	Yr_Mean = Mean(*Vec_Yr);

	double sqrt2 = pow(2, 0.5);
	Sr =  sqrt2/Mean( DotPow(DotPow((*Vec_Xr - Xr_Mean), 2) + DotPow((*Vec_Yr - Yr_Mean), 2), 0.5) );

	double tmpData1[] = {	Sr, 0 , -Sr*Xr_Mean,   
							0 , Sr, -Sr*Yr_Mean,    
							0 , 0 , 1 };
	Mat_Tr	  = new fMatrix(tmpData1, 3, 3);

	*Vec_X2 = Mat_Tr->GetMatVal(0,0)* *Vec_Xr + Mat_Tr->GetMatVal(0,2);
	*Vec_Y2 = Mat_Tr->GetMatVal(1,1)* *Vec_Yr + Mat_Tr->GetMatVal(1,2);

	//=================================Step2================================
	*Vec_Xs = Mat_X3->GetCol(0);
	*Vec_Ys = Mat_X3->GetCol(1);
	Xs_Mean = Mean(*Vec_Xs);
	Ys_Mean = Mean(*Vec_Ys);

	Ss = sqrt2/Mean( DotPow(DotPow((*Vec_Xs - Xs_Mean), 2) + DotPow((*Vec_Ys - Ys_Mean), 2), 0.5) );
	
	double tmpData2[] = {	Ss, 0 , -Ss*Xs_Mean,   
							0 , Ss, -Ss*Ys_Mean,    
							0 , 0 , 1 };
	Mat_Ts	  = new fMatrix(tmpData2, 3, 3);

	*Vec_X1 = Mat_Ts->GetMatVal(0,0)* *Vec_Xs + Mat_Ts->GetMatVal(0,2);
	*Vec_Y1 = Mat_Ts->GetMatVal(1,1)* *Vec_Ys + Mat_Ts->GetMatVal(1,2);

	//===============================Method 2===================================
	double tmp_B_Data[] = 
	{  
		Vec_X1->GetVectorVal(0), Vec_Y1->GetVectorVal(0), 1, 0, 0, 0, -(Vec_X2->GetVectorVal(0))*Vec_X1->GetVectorVal(0), -(Vec_X2->GetVectorVal(0))* Vec_Y1->GetVectorVal(0),
		0, 0, 0, Vec_X1->GetVectorVal(0), Vec_Y1->GetVectorVal(0), 1, -(Vec_Y2->GetVectorVal(0))*Vec_X1->GetVectorVal(0), -(Vec_Y2->GetVectorVal(0))* Vec_Y1->GetVectorVal(0)
	};
	double tmp_Y_Data[] = 
	{  
		Vec_X2->GetVectorVal(0), Vec_Y2->GetVectorVal(0)
	};

	B = new fMatrix(tmp_B_Data, 2, 8);
	fVector *tmp_Y = new fVector(tmp_Y_Data, 2);
	*Y = *tmp_Y;

	for (int i=1; i<Vec_Xr->GetSize() ; i++)
	{
		double tmp_B_Data[] = 
		{  
			Vec_X1->GetVectorVal(i), Vec_Y1->GetVectorVal(i), 1, 0, 0, 0, -(Vec_X2->GetVectorVal(i))*Vec_X1->GetVectorVal(i), -(Vec_X2->GetVectorVal(i))* Vec_Y1->GetVectorVal(i),
			0, 0, 0, Vec_X1->GetVectorVal(i), Vec_Y1->GetVectorVal(i), 1, -(Vec_Y2->GetVectorVal(i))*Vec_X1->GetVectorVal(i), -(Vec_Y2->GetVectorVal(i))* Vec_Y1->GetVectorVal(i)
		};
		double tmp_Y_Data[] = 
		{  
			Vec_X2->GetVectorVal(i), Vec_Y2->GetVectorVal(i)
		};
		*B = V_Cat(*B, Generate_tmpB(tmp_B_Data));
		*Y = Cat(  *Y, Generate_tmpY(tmp_Y_Data));
	}

	//======================LS======================
	CParamEstimator *PE2 = new CParamEstimator;
	fVector* tmp_Vec_LS = new fVector;
	fMatrix* tmp_Mat_LS = new fMatrix;

	// Get the Parameter
	*(PE2->LS1->pMat_H) = *B;	// Set Mat Param
	*(PE2->LS1->pVec_Z) = *Y;	// Set Vec Param 

	// Set the Parameter and Solve Method
	PE2->SetMethodParameters(LS, PE2->LS1);

	// Solve 
	tmp_Mat_LS = PE2->SolveOptParam(tmp_Vec_LS);

	//===============HomographyEstimation===============
	double  tmp_arr[]={1};
	fVector *tmp = new fVector(tmp_arr,1);
	*tmp_Vec_LS = Cat(*tmp_Vec_LS, *tmp);
	delete tmp;

	tmp_Mat_LS->ReShape(*tmp_Vec_LS ,3 ,3);		
	*tmp_Mat_LS = Inverse(*Mat_Tr)* *tmp_Mat_LS * *Mat_Ts;	

	//x_p process
	fMatrix *tmp_x_p1 = new fMatrix;
	tmp_x_p1->ReShape(*Vec_Xs, 1, Vec_Xs->GetSize());
	fMatrix *tmp_x_p2 = new fMatrix;
	tmp_x_p2->ReShape(*Vec_Ys, 1, Vec_Ys->GetSize());

	fMatrix x_p;
	x_p = V_Cat(*tmp_x_p1, *tmp_x_p2);
	x_p = V_Cat(x_p, Ones(1,Vec_Xs->GetSize()));
	x_p = *tmp_Mat_LS * x_p;

	x_p.SetRow(0, x_p.GetRow(0)/x_p.GetRow(2));
	x_p.SetRow(1, x_p.GetRow(1)/x_p.GetRow(2));
	x_p.SetRow(2, x_p.GetRow(2)/x_p.GetRow(2));

	//x_r process
	fMatrix EstiErr;

	fMatrix *tmp_x_r1 = new fMatrix;
	tmp_x_r1->ReShape(*Vec_Xr, 1, Vec_Xr->GetSize());		
	fMatrix *tmp_x_r2 = new fMatrix;
	tmp_x_r2->ReShape(*Vec_Yr, 1, Vec_Yr->GetSize());

	EstiErr = V_Cat(*tmp_x_r1, *tmp_x_r2);
	EstiErr = V_Cat(EstiErr, Ones(1,Vec_Xr->GetSize()));
	EstiErr = (EstiErr-x_p);

	
	double tmp_err_arr[]={0};
	tmp_err_arr[0] = TwoNorm( DotPow(EstiErr.GetCol(0), 2) );
	fVector err(tmp_err_arr, 1);
	fVector err_reg = err;
	
	for (int i=1;i<Vec_Xr->GetSize();i++)
	{
		err_reg = TwoNorm( DotPow(EstiErr.GetCol(i), 2) );
		err = Cat(err, err_reg);
	}
	Out_err = Mean(err);

	//=========================delete=========================
	delete Vec_Xr;		delete Vec_Yr;
	delete Vec_Xs;		delete Vec_Ys;
	delete Vec_X1;		delete Vec_Y1;
	delete Vec_X2;		delete Vec_Y2;
	delete Mat_Tr;		delete Mat_Ts;
	delete Mat_X3;		
	delete Mat_Ref_X;	delete tmp_Y;

	delete tmp_x_r1;
	delete tmp_x_p1;	delete tmp_x_p2;
	delete Y;			delete B;
	delete tmp_Vec_LS;

	//=========================delete=========================
	
	return tmp_Mat_LS;
}
//======================暫用 - 之後可改更彈性=============================
fMatrix CParamEstimator::Generate_tmpB(double* Data)
{
	fMatrix tmp_B(Data, 2, 8);
	return tmp_B;
}
fVector CParamEstimator::Generate_tmpY(double* Data)
{
	fVector tmp_Y(Data, 2);
	return tmp_Y;
}
//======================暫用 - 之後可改更彈性=============================
void CParamEstimator::LMV()
{
	/*cout<<"\n=========================LMV=============================\n";
	P_lmv = P_ls + Var_ls * Transp(Mat_H) * Inverse( Mat_H*Var_ls*Transp(Mat_H)+Vv )*(Vec_Z-Mat_H*P_ls);
	Var_lmv = Inverse( Inverse(Var_ls) + Transp(Mat_H) * Inverse(Vv)*Mat_H );

	cout<<"\nP_lmv = \n";
	P_lmv.Show(ColVec);

	cout<<"\nVar_lmv = \n";
	Var_lmv.Show();*/
}





