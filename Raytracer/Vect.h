#ifndef VECT_H
#define VECT_H

template<typename T>
class Vect{
	public:
		Vect<T>();
		Vect<T>& operator=(const Vect<T>& a);
		void operator/=(const float rhs);
		void operator+=(const Vect<T> rhs);
		void operator-=(const Vect<T> rhs);
		const T& operator[] (int i) const;
		Vect<T> operator- () const;

		T& operator[] (int i);
		
		T norm();
		void normalize();
	private:
		T elm[3];

};

template<typename T>
inline Vect<T>::Vect() {

}

template<typename T>
inline const T& Vect<T>::operator[] (int i) const {
	assert(i >= 0);
	assert(i < 3);
	return elm[i];
}

template<typename T>
inline T& Vect<T>::operator[] (int i) {
	assert(i >= 0);
	assert(i < 3);
	return elm[i];
}

template<typename T>
inline Vect<T> Vect<T>::operator- () const {
	Vect<T> c;
	c[0] = -elm[0];
	c[1] = -elm[1];
	c[2] = -elm[2];
	return c;
}
template<typename T>
T Vect<T>::norm() {
	T temp = dot(*this, *this);
	return sqrt(temp);
}

template<typename T>
void Vect<T>::normalize() {
	T temp = this->norm();
	(*this)[0] /= temp;
	(*this)[1] /= temp;
	(*this)[2] /= temp;
}

template<typename T>
Vect<T> cross(Vect<T>& a, Vect<T>& b) {
	Vect<T> c;
	/*c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = -(a[0] * b[2] - a[2] * b[0]);
	c[2] = a[0] * b[1] - a[1] * b[0];*/
	//left hand rule
	c[0] = b[1] * a[2] - b[2] * a[1];
	c[1] = -(b[0] * a[2] - b[2] * a[0]);
	c[2] = b[0] * a[1] - b[1] * a[0];
	return c;
}

template<typename T>
T dot(Vect<T>& a, Vect<T>& b) {
	T temp;
	temp = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	return temp;
}

template<typename T>
Vect<T>& Vect<T>::operator=(const Vect<T>& a) {
	if (this != &a) {
		(*this)[0] = a[0];
		(*this)[1] = a[1];
		(*this)[2] = a[2];
	}
	return *this;
}

template<typename T>
inline void Vect<T>::operator/=(const float rhs) {
	elm[0] /= rhs;
	elm[1] /= rhs;
	elm[2] /= rhs;
}

template<typename T>
inline void Vect<T>::operator+=(const Vect<T> rhs) {
	elm[0] += rhs.elm[0];
	elm[1] += rhs.elm[1];
	elm[2] += rhs.elm[2];
}

template<typename T>
inline void Vect<T>::operator-=(const Vect<T> rhs) {
	elm[0] -= rhs.elm[0];
	elm[1] -= rhs.elm[1];
	elm[2] -= rhs.elm[2];
}

template<typename T>
Vect<T> operator*(double a, Vect<T>& b) {
	Vect<T> c;
	c[0] = b[0] * a;
	c[1] = b[1] * a;
	c[2] = b[2] *a;
	return c;

}

template<typename T>
Vect<T> operator*(float a, Vect<T>& b)
{
	Vect<T> c;
	c[0] = b[0] * a;
	c[1] = b[1] * a;
	c[2] = b[2] * a;
	return c;
}

template<typename T>
Vect<T> operator*(Vect<T>& b, float a) {
	return operator*(a, b);
}

//pointwise multiplicatin
template<typename T>
Vect<T> operator*(Vect<T>& a, Vect<T>& b) {
	Vect<T> c;
	c[0] = a[0] * b[0];
	c[1] = a[1] * b[1];
	c[2] = a[2] * b[2];
	return c;
}

template<typename T>
Vect<T> operator+(Vect<T>& a, Vect<T>& b)
{
	Vect<T> c;
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
	return c;
}

template<typename T>
Vect<T> operator-(Vect<T>& a, Vect<T>& b) {
	Vect<T> c;
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
	return c;
}


#endif