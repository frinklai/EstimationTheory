#include "stdafx.h"
#include "fVector.h"
int fVector::nVecCount = 0;
const fVector fVector::Null(0);

fVector::fVector( const fVector & v1)
{
	elem = new Float[v1.size];
	size = v1.size;
	for (int i=0 ; i<v1.size ; i++)
		elem[i] = v1.elem[i];
	nVecCount++;
}

fVector::fVector( const Float *x, int n )
{
	nVecCount++;
	size = n;
	elem = new Float[size];
	for (int i=0 ; i<size ; i++)
	{
		*(elem+i) = *(x+i);
	}
}

fVector::fVector( int n, const Float *x )
{
	nVecCount++;
	size = n;
	elem = new Float[size];
	for (int i=0 ; i<size ; i++)
	{
		*(elem+i) = *(x+i);
	}
}
fVector::fVector( int n )
{
	nVecCount++;
	if (n>0)
	{
		size = n;
		elem = new Float[size];
	}
}

fVector::~fVector()
{
	nVecCount--;
	delete[] elem;
}

fVector::fVector( Float x, Float y, Float z)
{
	nVecCount++;
	elem = new Float[3];
	*(elem+0) = x;
	*(elem+1) = y;
	*(elem+2) = z;
}

fVector::fVector( Float x, Float y)
{
	nVecCount++;
	elem = new Float[2];
	*(elem+0) = x;
	*(elem+1) = y;
}

void fVector::Show(VecType Type) const
{
	if (Type == RowVec)
		for (int i=0 ; i< size ; i++)
			printf("%8.4f, ",elem[i]);
	else
		for (int i=0 ; i< size ; i++)
			printf("%8.4f\n",elem[i]);
	printf("\n");
	cout<<"size = "<<size;
}

void fVector::SetSize( int size)
{
	this->size = size;
}

void fVector::SetBlock( int i, int j, const fVector &v1 )
{
	for (int k=i,m=0;k<=j;k++,m++)
		elem[k] = v1.elem[m];
}
fVector& fVector::Swap( int i, int j )
{
	Float tmp = 0;
	tmp = this->elem[i];
	this->elem[i] = this->elem[j];
	this->elem[j] = tmp;
	return *this;
}
fVector fVector::GetBlock( int i, int j ) const
{
	fVector tmp(0);
	tmp.elem = new Float[(j-i)+1];
	for (int k=i,m=0;k<=j;k++,m++)
		tmp.elem[m] = this->elem[k];

	return fVector(tmp.elem, (j-i)+1);
}


void ShowVector ( const fVector & v1, VecType Type)
{
	if (Type == RowVec)
		for (int i=0 ; i< v1.size ; i++)
			printf("%.2f, ",v1.elem[i]);
	else
		for (int i=0 ; i< v1.size ; i++)
			printf("%.2f\n",v1.elem[i]);

		printf("\n");
}
fVector &fVector::operator=( const fVector & v1)
{
	elem = NULL;
	elem = new Float[v1.size];
	size = v1.size;
	for (int i=0 ; i<v1.size ; i++)
		elem[i] = v1.elem[i];

	return *this;
}
//======================================
//===========friend functions===========
//======================================
fVector operator + (const fVector &v1, const fVector &v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i] + v2.elem[i];

	return fVector(tmp.elem, v1.size);
}

fVector operator + (const fVector &v1, Float n)
{
	fVector tmp(1);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i] + n;
	return fVector(tmp.elem, v1.size);
}

fVector operator + (Float n, const fVector &v1)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = n +  v1.elem[i];
	return fVector(tmp.elem, v1.size);
}

fVector operator - (const fVector &v1, const fVector &v2)
{
	return v1+(-v2);
}

fVector operator - (const fVector &v1)
{
	return -1*v1;
}

fVector operator - (const fVector &v1, Float n)
{
	fVector tmp(1);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i]-n;
	return fVector(tmp.elem, v1.size);
}

fVector operator - (Float n, const fVector &v1)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = n -  v1.elem[i];
	return fVector(tmp.elem, v1.size);
}

fVector operator * (const fVector &v1, Float n)
{
	fVector tmp(1);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i]*n;
	return fVector(tmp.elem, v1.size);
}

fVector operator * (Float n, const fVector &v1)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = n*v1.elem[i];
	return fVector(tmp.elem, v1.size);
}

fVector operator / (const fVector &v1, Float n)
{
	return v1 * (1/n);
}

fVector operator / (const fVector &v1, const fVector &v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i] / v2.elem[i];

	return fVector(tmp.elem, v1.size);
}

fVector operator / (Float n, const fVector &v1)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = n /  v1.elem[i];
	return fVector(tmp.elem, v1.size);
}

double operator * (const fVector &v1, const fVector &v2)
{
	Float ans = 0;
	for (int i=1 ; i<=v1.size ; i++)
		ans += (v1.elem[i-1] * v2.elem[i-1]);

	return ans;
}

double det(Float a1, Float a2, Float b1, Float b2)
{
	return (a1*b2)-(a2*b1);
}

fVector operator ^ (const fVector &u, const fVector &v)
{
	fVector tmp(0);
	tmp.elem = new Float[v.size];
	for (int i=0, j=1, k=2 ; i<v.size ; i++, j++, k++)
	{
		if (j>(v.size-1))	j = 0;
		if (k>(v.size-1))	k = 0;
		fVector a(0);
		fVector b(0);
		a = /*const_cast<fVector&>*/(u);
		b = /*const_cast<fVector&>*/(v);
		tmp.elem[i] = det(a(j), a(k), b(j), b(k));
	}
	return fVector(tmp.elem, v.size);
}

fVector  operator^( const fVector &v1, int n )
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i = 0; i<v1.size; i++)
		tmp.elem[i] = pow(v1.elem[i], n);

	return fVector(tmp.elem, v1.size);
}

fVector& operator += ( fVector &v1, const fVector &v2)
{
	return v1 = v1+v2;
}

fVector& operator -= ( fVector &v1, const fVector &v2)
{
	return v1 = v1-v2;
}

fVector& operator *= ( fVector &v1, Float n)
{
	return v1 = v1*n;
}

fVector& operator /= ( fVector &v1, Float n)
{
	return v1 = v1/n;
}

fVector Min(const fVector &v1, const fVector &v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = (v1.elem[i] <= v2.elem[i]) ?v1.elem[i] :v2.elem[i] ;
	return fVector(tmp.elem, v1.size);
}

fVector Max(const fVector &v1, const fVector &v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = (v1.elem[i] >= v2.elem[i]) ?v1.elem[i] :v2.elem[i] ;
	return fVector(tmp.elem, v1.size);
}

double OneNorm(const fVector &v1)
{
	double ans=0;
	for (int i=0 ; i<v1.size ; i++)
		ans = ans + v1.elem[i];
	return ans;
}

double fabs_OneNorm(const fVector &v1)
{
	double ans=0;
	for (int i=0 ; i<v1.size ; i++)
		ans = ans + fabs(v1.elem[i]);
	return ans;
}

double Mean(const fVector &v1)
{
	return OneNorm(v1)/(v1.size);
}

double Var( const fVector & v1)
{
	return OneNorm(((v1 - Mean(v1))^2))/(v1.size-1);
}

double Var( const fVector & v1, const fVector & v2)
{
	return OneNorm(DotMul((v1 - Mean(v1)),(v2 - Mean(v2))))/(v2.size-1);
}

double Std( const fVector & v1)
{
	return sqrt(Var(v1));
}

fVector Sqrt( const fVector & v1)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
	{
		tmp.elem[i] = (Float)sqrt(v1.elem[i]);
	}
	return fVector(tmp.elem, v1.size);
}

double Dist( const fVector & v1, const fVector & v2)
{
	Float ans = 0;
	return sqrt(ans += (v1-v2)*(v1-v2));
}

double TwoNorm(const fVector &v1)
{
	Float ans = 0;
	return sqrt(ans += (v1*v1));
}

double TwoNormSqr(const fVector &v1)
{
	Float ans = 0;
	return sqrt(sqrt(ans += (v1*v1)));
}

double Angle(const fVector &v1, const fVector &v2)
{
	return acos((v1*v2)/( TwoNorm(v1)*TwoNorm(v2) ));
}

fVector Normalize( const fVector & v1)
{
	double len = TwoNorm(v1);
	return fVector((v1/len).elem, v1.size);
}

fVector Cat( const fVector &v1, const fVector &v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size+v2.size];
	tmp.SetBlock(0, v1.size-1, v1);
	tmp.SetBlock(v1.size, v1.size+(v2.size-1), v2);
	return fVector( tmp.elem, v1.size+v2.size );
}

fVector operator<<=(const fVector &v1, const fVector &v2 )//Combine 2 Vectors to 1 Vector
{
	return Cat(v1,v2);
}


int GetnVecCount()
{
	return fVector::nVecCount;
}

void fVector::operator=( Float n)
{
	for(int i=0;i<size;i++)
		elem[i] = n;
}

fVector DotMul(const fVector &v1, const fVector & v2)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = v1.elem[i] * v2.elem[i];

	return fVector(tmp.elem, v1.size);
}
fVector DotPow(const fVector &v1, double n)
{
	fVector tmp(0);
	tmp.elem = new Float[v1.size];
	for (int i=0 ; i<v1.size ; i++)
		tmp.elem[i] = pow(v1.elem[i], n);

	return fVector(tmp.elem, v1.size);
}
fVector fVector::LoadVector(char* tmp_FilePath)
{
	fstream file;
	char *str;
	int id;
	double *data;
	int  size;
	char *tmp = new char[1];

	file.open(tmp_FilePath, ios::in); //將檔案開啟為輸入狀態

	if(!file) 
	{
		cerr << "Can't open file!" << endl;
		system("PAUSE");
		exit(1); //在不正常情形下，中斷程式的執行
	}

	//Determine Size
	if(file >>tmp)
		size = atoi(tmp);
	else
		cout<<"error";

	str  = new char[size];
	data = new double[size];
	
	//Get Data
	for(int i=0 ; file >> str ; i++)//讀取記錄，若讀取至檔案結尾則傳回0
	{
		data[i] = atof(str);
	}

	delete[] str;
	file.close();

	return fVector(data, size);

}