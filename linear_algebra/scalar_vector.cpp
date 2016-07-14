#include "../common/common_shared.h"
#include "scalar_vector.h"

//---------------------------------------------------------------------//
ScalarVector::ScalarVector(int length)
{
	xassert(length > 0);
	m_vector = new float[length];
	m_length = length;
}

ScalarVector::ScalarVector(const ScalarVector& v)
{
	xassert(v.m_length > 0)
	m_length = v.m_length;
	m_vector = new float[m_length];
	memcpy(m_vector, v.m_vector, sizeof(float) * m_length);
}

ScalarVector::~ScalarVector()
{
	delete[] m_vector;
}

void ScalarVector::CopyFrom(const ScalarVector& v)
{
	xassert(m_length == v.m_length);
	memcpy(m_vector, v.m_vector, sizeof(float) * m_length);
}

//---------------------------------------------------------------------//
void ScalarVector::Set(int n, float val)
{
	xassert(n >= 0 && n < m_length);
	m_vector[n] = val;
}

float ScalarVector::Get(int n) const
{
	xassert(n >= 0 && n < m_length);
	return m_vector[n];
}

float ScalarVector::Norm() const
{
	return sqrtf(Norm2());
}

float ScalarVector::Norm2() const
{
	float sum = 0;
	for (int i = 0; i < m_length; i++)
	{
		sum += m_vector[i] * m_vector[i];
	}
	return sum;
}

//---------------------------------------------------------------------//
void ScalarVector::Add(const ScalarVector& v)
{
	xassert(m_length == v.GetLength());
	for (int ii = 0; ii < m_length; ii++)
	{
		m_vector[ii] += v.Get(ii);
	}
}

void ScalarVector::Multiply(float val)
{
	for (int ii = 0; ii < m_length; ii++)
	{
		m_vector[ii] *= val;
	}
}

void VectorMult(ScalarVector* result, const ScalarVector& v, const float m)
{
	xassert(result->GetLength() == v.GetLength());
	for (int ii = 0; ii < v.GetLength(); ii++)
	{
		result->Set(ii, v.Get(ii) * m);
	}
}

void VectorAdd(ScalarVector* result, const ScalarVector& a, const ScalarVector& b)
{
	xassert(result->GetLength() == a.GetLength());
	xassert(result->GetLength() == b.GetLength());

	for (int ii = 0; ii < a.GetLength(); ii++)
	{
		result->Set(ii, a.Get(ii) + b.Get(ii));
	}
}

void VectorSubtract(ScalarVector* result, const ScalarVector& a, const ScalarVector& b)
{
	xassert(result->GetLength() == a.GetLength());
	xassert(result->GetLength() == b.GetLength());

	for (int ii = 0; ii < a.GetLength(); ii++)
	{
		result->Set(ii, a.Get(ii) - b.Get(ii));
	}
}

float DotProd(const ScalarVector& a, const ScalarVector& b)
{
	xassert(a.GetLength() == b.GetLength());
	float sum = 0.f;

	for (int ii = 0; ii < a.GetLength(); ii++)
	{
		sum += a.Get(ii) * b.Get(ii);
	}

	return sum;
}

float Norm2(const EVector& a)
{
	float sum = 0;
	for (int i = 0; i < a.rows(); i++)
	{
		sum += a(i) * a(i);
	}
	return sum;
}
