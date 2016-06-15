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
