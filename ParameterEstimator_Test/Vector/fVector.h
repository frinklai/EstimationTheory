/*
 *	fVector.h
 *
 *	Description:
 *		Basic vector class with some associated methods.
 *
 *
 *		
 * 	History:
 *	 	Author			Date			Modify Reason		
 *		----------------------------------------------------------------
 *		Chi-Yi Tsai		2015/02/26		File Creation    
 *
 */
#ifndef __VECTOR_INCLUDED__
#define __VECTOR_INCLUDED__

#ifndef DOUBLE_PRECISION
#define DOUBLE_PRECISION
#endif

#ifndef Float
#ifdef DOUBLE_PRECISION
#define Float double
#else
#define Float float
#endif
#endif

enum VecType {ColVec = 1, RowVec};
 
class fVector 
{
/*-------------------------------------------------------------------------*
 *                                                                         *
 *  FRIEND OPERATORS                                                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
/*ok*/friend fVector  operator +  ( const fVector &, const fVector & );
/*ok*/friend fVector  operator -  ( const fVector &, const fVector & ); // Binary minus.
/*ok*/friend fVector  operator -  ( const fVector &                 );  // Unary minus.
/*ok*/friend fVector  operator -  ( const fVector &, Float ); 
/*ok*/friend fVector  operator -  (		Float	, const fVector & );
/*ok*/friend fVector  operator *  ( const fVector &,        Float   );
/*ok*/friend fVector  operator *  (       Float   , const fVector & );
/*ok*/friend fVector  operator /  ( const fVector &,        Float   );
/*ok*/friend fVector  operator /  ( const fVector &, const fVector & ); // Element-wise division
/*ok*/friend double   operator *  ( const fVector &, const fVector & ); // Inner-product between two vectors
/*ok*/friend fVector  operator ^  ( const fVector &, const fVector & ); // Cross-product between two vectors
/*ok*/friend fVector& operator += (       fVector &, const fVector & );
/*ok*/friend fVector& operator -= (       fVector &, const fVector & );
/*ok*/friend fVector& operator *= (       fVector &,        Float   );
/*ok*/friend fVector& operator /= (       fVector &,        Float   );

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  FRIEND FUNCTIONS                                                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
/*ok*/friend fVector  Min         ( const fVector &, const fVector & ); // Element-wise minimum-element extraction between two vectors
/*ok*/friend fVector  Max         ( const fVector &, const fVector & ); // Element-wise maximum-element extraction between two vectors
/*ok*/friend double   Dist        ( const fVector &, const fVector & ); // Returns two norm distance between two vectors
/*ok*/friend fVector  Normalize   ( const fVector & ); // Normalizes a vector into an unit vector
/*ok*/friend double   OneNorm     ( const fVector & ); // Returns one norm value of a vector
/*ok*/friend double   TwoNorm     ( const fVector & ); // Returns two norm value of a vector
/*ok*/friend double   TwoNormSqr  ( const fVector & ); // Returns square of the two norm value of a vector

/*ok*/friend fVector  Sqrt		( const fVector & ); // Element-wise square root of a vector
/*ok*/friend double   Mean		( const fVector & ); // Mean value of a vector.
/*ok*/friend double   Var			( const fVector & ); // Variance of a vector. 
/*ok*/friend double   Std			( const fVector & ); // Standard derivation of a vector.    	

/*ok*/friend void     ShowVector  ( const fVector &, VecType Type = ColVec );
 
public:
/*-------------------------------------------------------------------------*
 *                                                                         *
 *  C O N S T R U C T O R S  & D E S T R U C T O R S                       *
 *                                                                         *
 *-------------------------------------------------------------------------*/
	
	 // Initinalize constructor.
/*ok*/fVector( int size = 0  );			
	 // Copy constructor.
/*ok*/fVector( const fVector &/*A = Null*/);			
	 // Assign constructor.
/*ok*/fVector( const Float *x, int n );	
/*ok*/fVector( int n, const Float *x );	
/*ok*/fVector( Float, Float );			
/*ok*/fVector( Float, Float, Float );		
/*ok*/~fVector();

    static  const fVector Null;

public:
/*ok*/fVector &operator=( const fVector & );
/*OK*/void    operator=( Float );
/*OK*/void    SetSize( int );

/*ok*/fVector &Swap( int i, int j );
/*ok*/fVector GetBlock( int i, int j ) const;
/*ok*/void    SetBlock( int i, int j, const fVector & );
/*ok*/void	  Show(VecType Type = RowVec) const;
	
private:
    int    size;
    Float* elem;

	static int nVecCount;

/*-------------------------------------------------------------------------*
 *                                                                         *
 *  Add Variables or Functions, by Yu-Cheng, Lai                           *
 *                                                                         *
 *-------------------------------------------------------------------------*/
public:
	inline void   SetVectorID(int id, Float value){*(elem+id) = value;}	
	inline double  GetVectorVal(int id){return *(elem+id);}
	inline Float*  GetVector(){return (elem);}
	inline int	  GetVectorSize(){return size;}
	friend double Angle(const fVector &v1, const fVector &v2);
	inline friend double det(Float, Float, Float, Float);
	inline double operator()(const int n ){return *(elem+n);}	
	friend fVector Cat( const fVector &, const fVector &  );		//Horizontal cat 2Vector
	friend fVector operator<<=(const fVector &, const fVector & );
	friend int GetnVecCount();
	inline int GetSize(){return size;}
	friend double   fabs_OneNorm( const fVector & ); // Returns one norm value of a vector
	friend double   Var			( const fVector &, const fVector &  ); 
	friend fVector DotMul(const fVector &, const fVector & );
	friend fVector DotPow(const fVector &, double n );
	fVector LoadVector(char*);
	friend fVector  operator +  ( const fVector &, Float );
	friend fVector  operator +  (		Float	, const fVector & );
	friend fVector  operator /  (		Float	, const fVector & );

private:
	friend fVector  operator^( const fVector &, int n );     //vector bit-wise power		
};

#endif // __VECTOR_INCLUDED__
