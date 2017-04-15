/*
 *	fMatrix.h
 *
 *	Description:
 *		Basic matrix class with some associated methods.
 *
 *
 *		
 * 	History:
 *	 	Author			Date			Modify Reason		
 *		----------------------------------------------------------------
 *		Chi-Yi Tsai		2015/02/26		File Creation    
 *
 */

#ifndef __MATRIX_INCLUDED__
#define __MATRIX_INCLUDED__
#include "..\Vector\fVector.h"
class fMatrix 
{
/*-------------------------------------------------------------------------*
 *                                                                         *
 *  FRIEND OPERATORS                                                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
// 1. A+B
/*OK*/friend fMatrix  operator +  ( const fMatrix &, const fMatrix & );
// 2. A-B
/*OK*/friend fMatrix  operator -  ( const fMatrix &                 );
/*OK*/friend fMatrix  operator -  ( const fMatrix &, const fMatrix & );
// 3. c*A or A*c
/*OK*/friend fMatrix  operator *  ( const fMatrix &,       Float    );
/*OK*/friend fMatrix  operator *  (       Float  ,  const fMatrix & );
// 4. A/c
/*OK*/friend fMatrix  operator /  ( const fMatrix &,       Float    );
// 5. A*B
/*OK*/friend fMatrix  operator *  ( const fMatrix &, const fMatrix & );
/*OK*/friend fVector  operator *  ( const fMatrix &, const fVector & );
/*OK*/friend fVector  operator *  ( const fVector &, const fMatrix & );

/*OK*/friend fMatrix& operator += (       fMatrix &, const fMatrix & );
/*OK*/friend fMatrix& operator -= (       fMatrix &, const fMatrix & );
/*OK*/friend fMatrix& operator *= (       fMatrix &,       Float    );
/*OK*/friend fMatrix& operator *= (       fMatrix &, const fMatrix & );
/*OK*/friend fVector& operator *= (       fVector &, const fMatrix & );
/*OK*/friend fMatrix& operator /= (       fMatrix &,       Float    );

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  FRIEND FUNCTIONS                                                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
/*OK*/friend fMatrix  Transp      ( const fMatrix & );	// Transpose of a matrix
/*OK*/friend fMatrix  AATransp    ( const fMatrix & );	// Computes A * Transp(A).
/*OK*/friend fMatrix  ATranspA    ( const fMatrix & );	// Computes Transp(A) * A.
/*OK*/friend fMatrix  Outer       ( const fVector &, const fVector & );  // Computes the outer product of two vectors.
/*OK*/friend fMatrix  Identity	( int nSize ); // Returns an nSizexnSize identity matrix.

/*OK*/friend fMatrix  Diag        ( const fVector & );// Returns the square matrix with the elements of the vector d along its diagonal.		把取出的對角線方陣化
/*OK*/friend fVector  Diag        ( const fMatrix & );// Returns the vector consisting of the diagonal elements of the matrix M				取出矩陣的對角線
/*OK*/friend fMatrix  Diag        ( Float, Float, Float );// Returns the 3 x 3 diagonal matrix with x, y, and z as its diagonal elements.		指定矩陣對角線的3個元素

/*OK*/friend double   Determinant ( const fMatrix & );// Computes the determinant of a square matrix
/*OK*/friend double   Trace       ( const fMatrix & );// Computes the trace of a square matrix											//對角線之合
/*OK*/friend double   OneNorm     ( const fMatrix & );// Computes the L1-norm of the matrix A, which is the maximum absolute column sum.
/*OK*/friend double   InfNorm     ( const fMatrix & );// Computes the Inf-norm of the matrix A, which is the maximum absolute row sum.

/*OK*/friend fMatrix  Inverse  ( const fMatrix & );// Computes the inverse of a square matrix.
/*OK*/friend fMatrix  Cholesky	( const fMatrix & );// Computes Cholesky decomposition of a square matrix.	
/*OK*/friend fVector  Mean		( const fMatrix & );// Computes column mean value of a matrix.	
/*OK*/friend fMatrix  Cov			( const fMatrix & );// Returns a covariance matrix of a square matrix.
/*OK*/friend fMatrix  Cov			( const fVector & );// Returns a covariance matrix of a vector, using outer product.		http://zh.wikipedia.org/wiki/%E5%8D%8F%E6%96%B9%E5%B7%AE%E7%9F%A9%E9%98%B5
friend void     SVDcmp		( fMatrix &AU, fVector &W, fMatrix &V); // Computes SVD decomposition of a matrix.

/*OK*/friend void     ShowMatrix  ( const fMatrix & );// Print a matrix on screen.

public:
/*-------------------------------------------------------------------------*
 *                                                                         *
 *  C O N S T R U C T O R S  & D E S T R U C T O R S                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
	// Initinalize constructor.
    /*OK*/fMatrix( int n_rows = 0, int n_cols = 0 );
	// Assign constructor.
	/*OK*/fMatrix( Float *Array, int n_rows , int n_cols  );
	/*OK*/fMatrix( int n_rows , int n_cols , Float *Array );
	// Copy constructor.
	/*OK*/fMatrix( const fMatrix &/*A = Null */);  
	// Destructor
   /*OK*/~fMatrix( void );

    static  const fMatrix Null;

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  A S S I G N M E N T    O P E R A T O R S                               *
 *                                                                         *
 *-------------------------------------------------------------------------*/
	// 6. A=B
/*OK*/fMatrix  &operator=( const fMatrix &M );
/*OK*/fMatrix  &operator=( Float s );

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  MATRIX OPERATION FUNCTIONS                                             *
 *                                                                         *
 *-------------------------------------------------------------------------*/
	// 7. Swap
/*OK*/fMatrix &SwapRows( int i1, int i2 );
/*OK*/fMatrix &SwapCols( int j1, int j2 );
	// 8. Inverse
/*OK*/fMatrix &Inv(void);

/*OK*/void     SetCol( int col, const fVector & );
/*OK*/void     SetRow( int row, const fVector & );
/*OK*/void     SetBlock( int imin, int imax, int jmin, int jmax, const fMatrix & );
/*OK*/void     SetBlock( int imin, int imax, int jmin, int jmax, const fVector & );
/*OK*/void     SetSize( int rows, int cols = 0 );

/*OK*/fVector  GetCol( int col ) const;
/*OK*/fVector  GetRow( int row ) const;
/*OK*/fMatrix  GetBlock( int imin, int imax, int jmin, int jmax ) const;

/*OK*/void	   Show() const;

private:
    int    rows; // Number of rows in the matrix.
    int    cols; // Number of columns in the matrix.
    Float *elem; // Pointer to the actual data.

	static int nMatCount;

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  Add Variables or Functions, by Yu-Cheng, Lai                           *
 *                                                                         *
 *-------------------------------------------------------------------------*/
public:
	/*void Show();*/
	inline void ShowErrorMessage(int Errcode){cout<<"Error code = " << Errcode;}
	friend fMatrix V_Cat(const fMatrix &, const fMatrix & );
	friend fMatrix H_Cat(const fMatrix &, const fMatrix & );
	friend fMatrix Trans2Square(const fMatrix &);
	friend fVector  Row_Mean( const fMatrix & );
	fMatrix&  Row_Elim(int i, Float val, int j );//列消去 GJE用
	fMatrix&  Row_Elim(int i, Float val        );//列消去 GJE用
	double GetMatVal( int i, int j)const;
	void   SetMatVal( int i, int j ,Float val);
	double GetProdTrace( int i, int j );
	fMatrix Minor(int iRow,int iCol);
	fMatrix LoadMatrix(char*);
	void ReShape(fVector &v1, int row, int col);
	friend fMatrix  Ones( int row,int col );
private:
	int size;

};

#endif // __MATRIX_INCLUDED__
