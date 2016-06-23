
#include "../common/common_shared.h"
#include "gradient.h"

//----------------------------------------------------------------------------------------//
Gradient::Gradient()
	: m_length(0)
	, m_funcs(nullptr)
{}

Gradient::Gradient(int n)
{
	xassert(n > 0);

	m_length = n;
	m_funcs = new ScalarF[m_length];
}

Gradient::~Gradient()
{
	delete[] m_funcs;
}

//----------------------------------------------------------------------------------------//
void Gradient::Set(int i, ScalarF f)
{
	xassert(i >= 0 && i < m_length);
	m_funcs[i] = f;
}

void Gradient::Evaluate(const ScalarVector& input, ScalarVector* output) const
{
	xassert(input.GetLength() == output->GetLength());

	for (int ii = 0; ii < input.GetLength(); ii++)
	{
		float res = (m_funcs[ii])(input);
		output->Set(ii, res);
	}
}

//----------------------------------------------------------------------------------------//
EGradient::EGradient()
	: m_length(0)
	, m_funcs(nullptr)
{}

EGradient::EGradient(int n)
{
	xassert(n > 0);

	m_length = n;
	m_funcs = new ScalarFunc[m_length];
}

EGradient::~EGradient()
{
	delete[] m_funcs;
}

//----------------------------------------------------------------------------------------//
void EGradient::Set(int i, ScalarFunc f)
{
	xassert(i >= 0 && i < m_length);
	m_funcs[i] = f;
}

void EGradient::Evaluate(const EVector& input, EVector* output) const
{
	xassert(input.rows() == output->rows());

	for (int ii = 0; ii < input.rows(); ii++)
	{
		float res = (m_funcs[ii])(input);
		(*output)(ii) = res;
	}
}