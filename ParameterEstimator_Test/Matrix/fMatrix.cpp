#include "stdafx.h"
//#include "fVector.h"
#include "fMatrix.h"

const fMatrix fMatrix::Null(0,0);
int fMatrix::nMatCount = 0;
fMatrix::fMatrix( int n_rows , int n_cols  )
{
	if (n_rows * n_cols > 0)
	{
		rows = n_rows;
		cols = n_cols;
		size = n_rows * n_cols;
		elem = new Float[size];
	}
	nMatCount++;
}

fMatrix::fMatrix( const fMatrix &m1 )
{
	elem = new Float[m1.size];
	size = m1.size;
	cols = m1.cols;
	rows = m1.rows;
	for (int i=0 ; i<m1.size ; i++)
		elem[i] = m1.elem[i];
	nMatCount++;
}

fMatrix::fMatrix( Float *Array, int n_rows , int n_cols  )
{
	rows = n_rows;
	cols = n_cols;
	size = n_rows * n_cols;
	elem = new Float[size];
	for (int i=0 ; i<size ; i++)
		*(elem+i) = *(Array+i);
	nMatCount++;
}

fMatrix::fMatrix( int n_rows , int n_cols , Float *Array )
{
	
	rows = n_rows;
	cols = n_cols;
	size = n_rows * n_cols;
	elem = new Float[size];
	for (int i=0 ; i<size ; i++)
		*(elem+i) = *(Array+i);
	nMatCount++;
}


fMatrix::~fMatrix()
{
	delete[] elem;
	nMatCount--;
}

void fMatrix::Show()const
{
	for (int r = 0, k=0 ; r<rows ; r++ )
	{
		for (int c = 0 ; c<cols ; c++,k++)
		{
			printf("%7.4f ", elem[k]);
		}
		cout<<endl;
	}
}


fMatrix Transp( const fMatrix & m1)
{
	fMatrix tmp(m1.cols, m1.rows);/*(0, 0)*/;
	//fMatrix tmp;
	int w = m1.rows;	int h = m1.cols;		
	//tmp.elem = new Float[m1.rows * m1.cols];

	for(int i=0,k=0 ; i<w ; i++)
	{
		for(int j=0 ; j<h ; j++,k++)
		{
			tmp.elem[j*w + i] = m1.elem[k];
		}
	}
	return tmp;
	//return fMatrix(tmp.elem, m1.cols, m1.rows);
}

void ShowMatrix ( const fMatrix & m1 )
{
	for(int i=0,k=0 ; i<m1.rows ; i++)
	{
		for(int j=0 ; j<m1.cols ; j++,k++)
		{
			printf("%.4f, ",m1.elem[k]);
		}
		cout<<endl;
	}
}


fMatrix  operator +  ( const fMatrix &m1, const fMatrix &m2 )
{
	fMatrix tmp(0, 0);
	if ((m1.rows!=m1.cols)||(m2.rows!=m2.cols))
	{
		fMatrix N_m1(0,0);		fMatrix N_m2(0,0);
		if(m1.rows!=m1.cols)	N_m1 = Trans2Square(m1);
		if(m2.rows!=m2.cols)	N_m2 = Trans2Square(m2);

		tmp.elem = new Float[N_m1.cols*N_m1.rows];
		int k=0;
		for(int i=0; i<N_m1.rows ; i++)
		{
			for(int j=0 ; j<N_m1.cols ; j++)
			{
				tmp.elem[k] = N_m1.elem[k]+N_m2.elem[k];
				k++;
			}
		}
		return fMatrix(tmp.elem, N_m1.rows, N_m1.cols);
	}
	else
	{
		tmp.elem = new Float[m1.cols*m1.rows];
		int k=0;
		for(int i=0; i<m1.rows ; i++)
		{
			for(int j=0 ; j<m1.cols ; j++)
			{
				tmp.elem[k] = m1.elem[k]+m2.elem[k];
				k++;
			}
		}
		return fMatrix(tmp.elem, m1.rows, m1.cols);
	}
}

fMatrix  operator - ( const fMatrix &m1, const fMatrix &m2 )
{
	fMatrix tmp(0, 0);
	if ((m1.rows!=m1.cols)||(m2.rows!=m2.cols))
	{
		fMatrix N_m1(0,0);		fMatrix N_m2(0,0);
		if(m1.rows!=m1.cols)	N_m1 = Trans2Square(m1);
		if(m2.rows!=m2.cols)	N_m2 = Trans2Square(m2);

		tmp.elem = new Float[N_m1.cols*N_m1.rows];
		int k=0;
		for(int i=0; i<N_m1.rows ; i++)
		{
			for(int j=0 ; j<N_m1.cols ; j++)
			{
				tmp.elem[k] = N_m1.elem[k]-N_m2.elem[k];
				k++;
			}
		}
		return fMatrix(tmp.elem, N_m1.rows, N_m1.cols);
	}
	else
	{
		tmp.elem = new Float[m1.cols*m1.rows];
		int k=0;
		for(int i=0; i<m1.rows ; i++)
		{
			for(int j=0 ; j<m1.cols ; j++)
			{
				tmp.elem[k] = m1.elem[k]-m2.elem[k];
				k++;
			}
		}
		return fMatrix(tmp.elem, m1.rows, m1.cols);
	}
}

fMatrix  operator * ( const fMatrix &m1, const fMatrix &m2 )
{
	fMatrix tmp(0, 0);	
	tmp.elem = new Float[m1.rows*m2.cols];

	int m = m1.rows;
	int p = m2.cols ;
	int n = m2.rows;;
	int id=0;
	for(int i=0;i<m;i++)		//列 m
	{
		for(int j=0;j<p;j++)	//行 p
		{
			tmp.elem[id]=0;
			for(int k=0;k<n;k++)//   n
			{
				tmp.elem[id] += (  m1.elem[i*m1.cols+k] * m2.elem[k*m2.cols+j]);
			}
			id++;
		}
	}
	tmp.rows = m;
	tmp.cols = p;
	return fMatrix(tmp.elem, tmp.rows, tmp.cols);
}

fVector  operator *  ( const fMatrix &m1, const fVector &v1 )
{
	fMatrix tmp(0, 0);	
	fVector tmp_V(0);	tmp_V = v1;

	int m = m1.rows;
	int n = tmp_V.GetSize();
	int p = 1;
	tmp.elem = new Float[m1.rows*p];
	
	int id=0;
	for(int i=0;i<m;i++)		//列 m
	{
		for(int j=0;j<p;j++)	//行 p
		{
			tmp.elem[id]=0;
			for(int k=0;k<n;k++)//   n
			{
				tmp.elem[id] += (  m1.elem[i*m1.cols+k] * tmp_V.GetVectorVal(k));
			}
			id++;
		}
	}
	tmp.rows = m;
	tmp.cols = p;
	return fVector(tmp.elem, m);
}

fVector  operator *  ( const fVector &v1, const fMatrix &m1 )
{
	fMatrix tmp(0, 0);	
	fVector tmp_V(0);	tmp_V = v1;

	int m = 1;
	int n = tmp_V.GetSize();//3
	int p = m1.cols;
	tmp.elem = new Float[n];

	int id=0;
	for(int i=0;i<m;i++)		//列 m
	{
		for(int j=0;j<p;j++)	//行 p
		{
			tmp.elem[id]=0;
			for(int k=0;k<n;k++)//   n
			{
				tmp.elem[id] += ( tmp_V.GetVectorVal(k) *  m1.elem[k*n+j]);
			}
			id++;
		}
	}
	tmp.rows = m;
	tmp.cols = p;
	return fVector(tmp.elem, n);
}

fMatrix  operator - ( const fMatrix &m1)
{
	return -1*m1;
}

fMatrix  operator * ( const fMatrix &m1, Float n)
{
	fMatrix tmp(0, 0);	
	tmp.elem = new Float[m1.rows*m1.cols];

	for(int i=0,k=0 ; i<m1.rows ; i++)
	{
		for(int j=0 ; j<m1.cols ; j++,k++)
		{
			tmp.elem[k] = m1.elem[k] * n;
		}
	}
	return fMatrix(tmp.elem, m1.rows, m1.cols);
}

fMatrix  operator * ( Float n, const fMatrix &m1)
{
	return m1*n;
}

fMatrix  operator / (const fMatrix &m1,  Float n)
{
	return m1*(1/n);
}

fMatrix& operator += (       fMatrix &m1, const fMatrix &m2 )
{
	return m1 = m1 + m2;
}

fMatrix& operator -= (       fMatrix &m1, const fMatrix &m2 )
{
	return m1 = m1 - m2;
}

fMatrix& operator *= (       fMatrix &m1, Float n )
{
	return m1 = m1*n;
}
fVector& operator *= (       fVector &v1, const fMatrix &m1 )
{
	v1 = v1*m1;
	return v1;
}
fMatrix& operator /= (       fMatrix &m1, Float n )
{
	return m1 = m1*(1/n);
}

fMatrix &fMatrix::operator=( const fMatrix &M )
{
	elem = NULL;
	elem = new Float[M.size];
	size = M.size;
	rows = M.rows;
	cols = M.cols;
	for(int i=0,k=0 ; i<rows ; i++)
	{
		for(int j=0 ; j<cols ; j++,k++)
		{
			elem[k] = M.elem[k];
		}
	}
	return *this;
}

fMatrix& operator *= (       fMatrix &m1, const fMatrix &m2 )
{
	m1 = m1*m2;
	m1.cols = m2.cols;
	return m1;
}

fMatrix  ATranspA( const fMatrix & m1)	// Computes Transp(A) * A.
{
	return Transp(m1)*m1;
}
fMatrix  AATransp( const fMatrix & m1)	// Computes A * Transp(A)
{
	return m1*Transp(m1);
}
fMatrix V_Cat(const fMatrix &m1, const fMatrix &m2)
{
	fMatrix tmp(0,0);
	int i=0;	int j=0;	int k=0;
	if (m1.rows == m2.rows)
	{
		tmp.cols = m1.cols;
		tmp.rows = m1.rows*2;
		tmp.size = tmp.cols*tmp.rows;
		tmp.elem = new Float[tmp.size];

		for (int r =0 ; r<tmp.rows ;r++)
		{
			for (int c=0 ; c<tmp.cols ; c++)
			{
				if(r >= m1.rows )
					tmp.elem[k] = m2.elem[i++];
				else 
					tmp.elem[k] = m1.elem[j++];
				k++;
			}
		}
	}
	else
	{
		tmp.rows = m1.rows+m2.rows;
		tmp.cols = (m1.cols>m2.cols)?m1.cols :m2.cols ;

		tmp.size = tmp.cols*tmp.rows;
		tmp.elem = new Float[tmp.size];

		for (int r =0 ; r<tmp.rows ;r++)
		{
			for (int c=0 ; c<tmp.cols ; c++)
			{
				if(k < m1.size )
					tmp.elem[k] = m1.elem[i++];
				else 
					tmp.elem[k] = m2.elem[j++];
				k++;
			}
		}
	}
	return fMatrix(tmp.elem, tmp.rows, tmp.cols);
}

fMatrix H_Cat(const fMatrix &m2, const fMatrix &m1 )
{
	fMatrix tmp(0,0);
	int i=0;	int j=0;	int k=0;
	if (m1.cols == m2.cols)
	{
		tmp.cols = m1.cols*2;
		tmp.rows = m1.rows;
		tmp.size = tmp.cols*tmp.rows;
		tmp.elem = new Float[tmp.size];

		for (int r =0 ; r<tmp.rows ;r++)
		{
			for (int c=0 ; c<tmp.cols ; c++)
			{
				if(c >= m1.cols )
					tmp.elem[k] = m1.elem[i++];
				else 
					tmp.elem[k] = m2.elem[j++];
				k++;
			}
		}
	}
	else
	{
		tmp.cols = m1.cols+m2.cols;
		tmp.rows = (m1.rows>m2.rows)?m1.rows :m2.rows ;

		tmp.size = tmp.cols*tmp.rows;
		tmp.elem = new Float[tmp.size];

		for (int r =0 ; r<tmp.rows ;r++)
		{
			for (int c=0 ; c<tmp.cols ; c++)
			{
				if (c >= m1.cols)
					tmp.elem[k] = m2.elem[i++];
				else
					tmp.elem[k] = m1.elem[j++];
				k++;
			}
		}
	}
	
	return fMatrix(tmp.elem, tmp.rows, tmp.cols);
}
fMatrix Trans2Square(const fMatrix &m1)
{
	fMatrix tmp(0,0);

	if(m1.cols > m1.rows)		
	{
		tmp.rows = m1.cols;
		tmp.cols = m1.cols;
	}	
	else
	{
		tmp.cols = m1.rows;
		tmp.rows = m1.rows;
	}
	tmp.size = tmp.cols * tmp.rows;
	tmp.elem = new Float[tmp.size];

	int delta_size = tmp.size - m1.size;
	int delta_row, delta_col;
	Float *delta_arr = new Float[delta_size]; 

	if (m1.rows > m1.cols)	//H_Cat
	{
		delta_row = m1.rows;
		delta_col = tmp.cols - m1.cols;
		for(int i=0;i<delta_size;i++)
			delta_arr[i] = 0;
		fMatrix delta_Mat(delta_arr, delta_row, delta_col);
		return H_Cat(m1, delta_Mat);
	}
	else					//V_Cat
	{
		delta_col = m1.cols;
		delta_row = tmp.rows - m1.rows;
		for(int i=0;i<delta_size;i++)
			delta_arr[i] = 0;
		fMatrix delta_Mat(delta_arr, delta_row, delta_col);
		return V_Cat(m1, delta_Mat);
	}

}

fVector fMatrix::GetCol( int col ) const
{
	Float* arr = new Float[rows];
	int k=0;
	for (int r = 0;r<rows;r++)
	{
		for(int c=0;c<cols;c++)
		{
			if (c == col)
				arr[k++] = elem[r*cols+c];
		}
	}
	return fVector(arr, k);
}

void fMatrix::SetCol( int col, const fVector &v1 )
{
	fVector tmp_V(0);
	tmp_V = v1;
	Float* arr = new Float[rows];
	int k=0;
	for (int r = 0;r<rows;r++)
	{
		for(int c=0;c<cols;c++)
		{
			if (c == col)
				elem[r*cols+c] = tmp_V.GetVectorVal(k++);
		}
	}
}
void fMatrix::SetSize( int rows, int cols )
{
	this->rows = rows;
	this->cols = cols;
}
void fMatrix::SetRow( int row, const fVector &v1 )
{
	fVector tmp_V(0);
	tmp_V = v1;
	Float* arr = new Float[rows];
	int k=0;
	for (int r = 0;r<rows;r++)
	{
		for(int c=0;c<cols;c++)
		{
			if (r == row)
				elem[r*cols+c] = tmp_V.GetVectorVal(k++);
		}
	}
}

fVector fMatrix::GetRow( int row ) const
{
	Float* arr = new Float[cols];
	int k=0;
	for (int r = 0;r<rows;r++)
	{
		for(int c=0;c<cols;c++)
		{
			if (r == row)
			{
				arr[k++] = elem[r*cols+c];
			}
		}
	}
	return fVector(arr, k);
}

//fMatrix fMatrix::GetBlock( int imin, int imax, int jmin, int jmax )const
//{
//	Float* arr = new Float[(imax-imin+1)*(jmax-jmin+1)];
//	int k = 0;
//	for( int i=0 ;i <=cols;i++)
//	{
//		for( int j=0;j<=rows;j++)
//		{
//			if((j>=jmin && j<=jmax) && (i>=imin && i<=imax))
//				arr[k++] = elem[i*cols+j];
//		}
//	}
//	return fMatrix(arr, imax-imin+1 ,jmax-jmin+1);
//}

fMatrix fMatrix::GetBlock( int imin, int imax, int jmin, int jmax )const
{
	Float* arr = new Float[(imax-imin+1)*(jmax-jmin+1)];
	int k = 0;
	for( int i=0;i<=rows;i++)
	{
		for( int j=0 ;j <=cols;j++)
		{
			if((j>=jmin && j<=jmax) && (i>=imin && i<=imax))
				arr[k++] = elem[i*cols+j];
		}
	}
	return fMatrix(arr, imax-imin+1 ,jmax-jmin+1);
}

void fMatrix::SetBlock( int imin, int imax, int jmin, int jmax, const fMatrix & m1)
{
	int k = 0;
	for( int i=0 ;i <=cols;i++)
	{
		for( int j=0;j<=rows;j++)
		{
			if((j>=jmin && j<=jmax) && (i>=imin && i<=imax))
				this->elem[i*cols+j] = m1.elem[k++];
		}
	}
}

void fMatrix::SetBlock( int imin, int imax, int jmin, int jmax, const fVector & v1)
{
	int k = 0;
	fVector tmp_V(0);
	tmp_V = v1;
	for( int i=0 ;i <=cols;i++)
	{
		for( int j=0;j<=rows;j++)
		{
			if((j>=jmin && j<=jmax) && (i>=imin && i<=imax))
				this->elem[i*cols+j] = tmp_V.GetVectorVal(k++);
		}
	}
}

fVector Diag( const fMatrix & m1)
{
	int tmp_size = (m1.rows>m1.cols)? m1.cols : m1.rows;
	int k=0;
	Float *arr = new Float[tmp_size];
	for(int i=0 ; i<m1.rows ; i++)
	{
		for(int j=0 ; j<m1.cols ; j++)
		{
			if (i==j)
				arr[k++] = m1.elem[i*m1.cols+j];
		}
	}
	return fVector(arr, k);
}

fMatrix Diag( Float x, Float y, Float z )
{
 Float *arr = new Float[9];
 int k=0;
 for (int i=0;i<3;i++)
 {
	 for (int j=0;j<3;j++)
	 {
		 if (i==j)
		 {
			 if(i==0)		arr[k++] = x;
			 else if(i==1)	arr[k++] = y;
			 else if(i==2)	arr[k++] = z;
		 }
		 else
			 arr[k++] = 0.0;
	 }
 }
 return fMatrix(arr, 3, 3);
}

fMatrix  Diag(const fVector & v1)
{
	fVector tmp_V(0);
	tmp_V = v1;
	int tmp_size = tmp_V.GetSize();

	Float *arr = new Float[tmp_size*tmp_size];
	int k=0;
	for (int i=0;i<tmp_size;i++)
	{
		for (int j=0;j<tmp_size;j++)
		{
			if (i==j)
				arr[i*tmp_size+j] = tmp_V.GetVectorVal(k++);
			else
				arr[i*tmp_size+j] = 0.0;
		}
	}
	return fMatrix(arr, tmp_size, tmp_size);
}

double Trace( const fMatrix & m1)
{
	double ans=0;
	for(int i=0 ; i<m1.rows ; i++)
		ans += m1.elem[i*m1.cols+i];
	return ans;
}

fMatrix  &fMatrix::operator=( Float s )
{
	int k=0;
	for (int r = 0;r<rows;r++)
	{
		for(int c=0;c<cols;c++)
		{
			elem[k++] = s;
		}
	}
	return *this;
}

fMatrix &fMatrix::SwapRows( int i1, int i2 )
{
	fVector v1(0);	fVector v2(0);	fVector tmp_V(0);
	v1 = this->GetRow(i1);
	v2 = this->GetRow(i2);

	tmp_V = v1;
	v1	  = v2;
	v2	  = tmp_V;

	this->SetRow(i1, v1);
	this->SetRow(i2, v2);
	
	return *this;
}

fMatrix &fMatrix::SwapCols( int j1, int j2 )
{
	fVector v1(0);	fVector v2(0);	fVector tmp_V(0);
	v1 = this->GetCol(j1);
	v2 = this->GetCol(j2);

	tmp_V = v1;
	v1	  = v2;
	v2	  = tmp_V;

	this->SetCol(j1, v1);
	this->SetCol(j2, v2);

	return *this;
}

fVector Mean( const fMatrix & m1)//Col_Mean
{
	Float *arr;
	arr = new Float[m1.cols];
	int i=0;
	for(i=0 ; i<m1.cols ;)
		arr[i++] = Mean(m1.GetCol(i));
	return fVector(arr, i);
}

fVector Row_Mean( const fMatrix & m1)
{
	Float *arr;
	arr = new Float[m1.rows];
	int i=0;
	for(i=0 ; i<m1.rows ;)
		arr[i++] = Mean(m1.GetRow(i));
	return fVector(arr, i);
}

fMatrix Outer( const fVector &v1, const fVector &v2 )
{
	fVector tmp_V1(0);	fVector tmp_V2(0);
	tmp_V1 = v1;		tmp_V2 = v2;

	int V1_Size = tmp_V1.GetVectorSize();
	int V2_Size = tmp_V2.GetVectorSize();
	if (V1_Size!=V2_Size)	
	{
		cout<<"\n\n\nERROR, V1 V2長度需相等!!!!!!!!\n\n\n";
		system("pause");
	}
	Float *arr = new Float[V1_Size*V2_Size];

	int k=0;
	for(int r = 0; r<V1_Size ; r++)
	{
		for (int c = 0; c<V2_Size ;c++)
		{
			arr[k++] = tmp_V1.GetVectorVal(r) *tmp_V2.GetVectorVal(c);
		}
	}
	return fMatrix(arr, V1_Size, V2_Size);
}

fMatrix  Identity( int nSize )
{
	Float *arr = new Float[nSize*nSize];
	int k=0;
	for (int r = 0; r<nSize ; r++)
	{
		for (int c=0; c<nSize ; c++)
		{
			if(r==c)	arr[k++] = 1;
			else		arr[k++] = 0;
		}
	}
	return fMatrix(arr, nSize, nSize);
}

fMatrix  Cholesky( const fMatrix & m1)
{
	Float* L = new Float[m1.rows*m1.cols];
	double s;

	for (int i = 0; i < m1.cols; i++)
	{
		for (int j = 0; j < m1.cols; j++) 
		{
			if(i<j)	
			{
				L[i * m1.cols + j] = 0;
			}
			else
			{
				s = 0;
				for (int k = 0; k < j; k++)
				{
					s += L[i * m1.cols + k] * L[j * m1.cols + k];
				}
				L[i * m1.cols + j] = (i == j) ? sqrt(m1.elem[i * m1.cols + i] - s)
											  :(1.0 / L[j * m1.cols + j] * (m1.elem[i * m1.cols + j] - s));
			}	
		}
	}
	return Transp(fMatrix(L, m1.cols, m1.cols));
	//return fMatrix(L, m1.cols, m1.cols);
}

double   OneNorm( const fMatrix & m1)
{
	double Max_Col_sum	= 0;
	double Col_Sum		= 0;

	for (int i = 0;i<m1.cols;i++)
	{
		Col_Sum = fabs_OneNorm(m1.GetCol(i));
		if (Col_Sum > Max_Col_sum)	
			Max_Col_sum	 = Col_Sum;
	}
	return Max_Col_sum;
}

double   InfNorm( const fMatrix & m1)
{
	double Max_Row_sum	= 0;
	double Row_Sum		= 0;

	for (int i = 0;i<m1.rows;i++)
	{
		Row_Sum = fabs_OneNorm(m1.GetRow(i));
		if (Row_Sum > Max_Row_sum)	
			Max_Row_sum	 = Row_Sum;
	}
	return Max_Row_sum;
}

double fMatrix::GetMatVal( int i, int j )const
{
	return elem[i*rows+j];
}

void fMatrix::SetMatVal( int i, int j ,Float val)
{
	elem[i*rows+j] = val;
}

double fMatrix::GetProdTrace( int i, int j )//計算det時會用到的，由左乘到右或是由右乘到左
{
	if(rows!=cols)	cout<<"ERROR";
	double product = 1;
	
	if ((i>j)||((i==0)&&(j==0)))
	{
		product = 1;
		for(int r=0 ; r<rows ;r++)
		{
			if(i>=rows)	i=0;
			if(j>=rows)	j=0;
			product *= this->GetMatVal(i++,j++);
		}
	}
	else
	{
		product = 1;
		for(int r=0 ; r<rows ;r++)
		{
			if(i>=rows)	i=0;
			if(j<0)		j=rows-1;
			int a= this->GetMatVal(i++,j--);
			product *= a;
		}
	}
	
	return product;
}

fMatrix  Cov( const fMatrix &m1)
{
	int k=0;int a=0;
	Float* arr;
	arr = new Float[m1.rows*m1.cols];

	for (int r=0;r<m1.rows;r++)
	{
		for (int c=0;c<m1.cols;c++)
		{
			if (r==c)
				arr[k++] = Var(m1.GetCol(c));
			else
				arr[k++] = Var(m1.GetCol(c), m1.GetCol(r));
		}
	}
	return fMatrix(arr, m1.rows, m1.cols);
}

fMatrix  Cov( const fVector &v1)
{
	fVector tmp_V(0);	
	tmp_V = v1;
	tmp_V = tmp_V - Mean(tmp_V);
	return Outer(tmp_V, tmp_V)/(tmp_V.GetSize()-1);
}

double  Determinant ( const fMatrix & m1)
{
	double DetValue = 0.0;
	fMatrix tmp(0,0);
	tmp.elem = new Float[m1.size];
	tmp = m1;

	if ( tmp.size > 0 && tmp.rows == tmp.cols )
	{
		if ( tmp.rows > 2 )
		{
			for (int r=0; r<tmp.rows; r++)
			{
				DetValue += pow( -1.0 , (double)r+2 )* tmp.GetMatVal(r,0) *( Determinant((tmp.Minor(r,0))) );
			}
		} 
		else if ( tmp.rows == 2 )
		{
			DetValue = tmp.elem[0]*tmp.elem[3] - tmp.elem[1]*tmp.elem[2];
		} 
		else
		{
			DetValue = tmp.elem[0];
		}	
	}
	return DetValue;
}

fMatrix fMatrix::Minor(int i,int j)
{
	Float* arr;
	int tmp_size = (rows-1)*(cols-1);
	arr = new Float[tmp_size];

	int k=0;
	for (int r=0;r<rows;r++)
	{
		for (int c=0;c<cols;c++)
		{
			if (!((r==i)||(c==j)))
				arr[k++] = elem[r*cols+c];
		}
	}
	return fMatrix(arr, rows-1, cols-1);
}

fMatrix&  fMatrix::Row_Elim(int i, Float val, int j )
{
	this->SetRow(j, ((val * this->GetRow(i)) + this->GetRow(j)) ) ;
	return *this;
}

fMatrix&  fMatrix::Row_Elim(int i, Float val )
{
	this->SetRow(i, val * this->GetRow(i)  ) ;
	return *this;
}

fMatrix  Inverse  ( const fMatrix & m1)
{
	if (m1.cols!=m1.rows)
	{
		cout<<"ERROR!!! 需為方陣\n";
		system("pause");
	}
	fMatrix I(m1.rows, m1.cols);	// 單位方陣矩陣
	fMatrix Aug(0,0);				// 增廣矩陣

	I	= Identity(m1.rows);
	Aug = H_Cat(m1, I);				// 初始化增廣矩陣

	//GJE

	//若第一個元素為0，跟不為0的元素換列
	if(Aug.elem[0]==0)
	{
		for(int k=1  ; k < Aug.rows ; k++)
		{
			if (Aug.elem[k*Aug.cols+0] != 0)
			{
				Aug.SwapRows(0,k);
			}
		}
	}

	//高斯-喬丹消去法
	for (int r=0,c=0;r<Aug.rows;r++,c++)
	{
		//檢查對角線元素
		if(r==c)
		{
			if (Aug.elem[r*Aug.cols+r] != 1)
				Aug.Row_Elim(r, 1/Aug.elem[r*Aug.cols+r] );	//化成 Leading 1

			for(int k=0  ; k < Aug.rows ; k++)
			{
				if ((Aug.elem[k*Aug.cols+r] != 0)&&(k!=r))
				{
					Float tmp = -1*(Aug.elem[k*Aug.cols+r]);
					Aug.Row_Elim(r, tmp, k);
				}
			}
		}
	}
	return Aug.GetBlock(0, Aug.rows-1, m1.cols, Aug.cols-1);
}

fMatrix &fMatrix::Inv(void)
{
	*this = Inverse(*this);
	return *this;
}

fMatrix fMatrix::LoadMatrix(char* tmp_FilePath)
{
	fstream file;
	char *str;
	int id;
	double *data;
	int  size[2];
	char *tmp = new char[2];

	file.open(tmp_FilePath, ios::in); //將檔案開啟為輸入狀態

	if(!file) 
	{
		cerr << "Can't open file!" << endl;
		system("PAUSE");
		exit(1); //在不正常情形下，中斷程式的執行
	}

	//Determine Size
	for(int i=0;i<2;i++)
	{
		if(file >>tmp)
			size[i] = atoi(tmp);
		else
		{
			cout<<"error";
			break;
		}
			
		if(i==1)
		{
			str  = new char[size[0]*size[1]];
			data = new double[size[0]*size[1]];
		}
	}
	//Get Data
	for(int i=0 ; file >> str ; i++)//讀取記錄，若讀取至檔案結尾則傳回0
	{
		data[i] = atof(str);
	}

	delete[] str;
	file.close();
	return fMatrix(data, size[0], size[1]);
	
}

void fMatrix::ReShape(fVector &v1, int row, int col)
{
	this->elem=NULL;
	this->elem = new Float[row*col];
	memcpy(this->elem, v1.GetVector(), (row*col)*sizeof(double));
	this->rows = row;
	this->cols = col;
}
fMatrix  Ones( int row,int col )
{
	fMatrix one(row, col);
	
	for (int r = 0; r<row ; r++)
	{
		for (int c=0; c<col ; c++)
		{
			one.SetMatVal(r,c,1);
		}
	}
	return one;
}