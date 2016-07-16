
#ifndef _SCALAR_VECTOR_H_
#define _SCALAR_VECTOR_H_

class ScalarVector
{

private:
	float*		m_vector;
	int			m_length;

public:
	ScalarVector(int length);
	ScalarVector(const ScalarVector& v);
	~ScalarVector();

	void CopyFrom(const ScalarVector& v);

	int GetLength() const { return m_length; }
	void Set(int n, float val);
	float Get(int n) const;

	float Norm() const;
	float Norm2() const;

	void Add(const ScalarVector& v);
	void Multiply(float val);
};

void VectorMult(ScalarVector* result, const ScalarVector& v, const float m);
void VectorAdd(ScalarVector* result, const ScalarVector& a, const ScalarVector& b);
void VectorSubtract(ScalarVector* result, const ScalarVector& a, const ScalarVector& b);

float DotProd(const ScalarVector& a, const ScalarVector& b);
float Norm2(const EVector& a);


#endif
