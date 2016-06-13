
class ScalarVector
{

private:
	float*		m_vector;
	int			m_length;

public:
	ScalarVector(int length);
	ScalarVector(const ScalarVector& v);
	~ScalarVector();


};