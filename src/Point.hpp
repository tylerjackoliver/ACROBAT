#ifndef _POINT_H_
#define _POINT_H_

#include <ostream>
#include <stdexcept>

template <typename Type>
struct Point
{
    std::vector<Type> state;

    // Class constructors
    Point() : state(6)
    {
        for (unsigned i = 0; i < state.size(); ++i) state[i] = 0.0;
    };

    Point(Type x_, Type y_, Type z_) : state(6)
    {
        this->state[0] = x_;
        this->state[1] = y_;
        this->state[2] = z_;
        this->state[3] = 0.0;
        this->state[4] = 0.0;
        this->state[5] = 0.0;
    }
    
    Point(Type x_, Type y_, Type z_, Type vx_, Type vy_, Type vz_) : state(6)
    {
        this->state[0] = x_;
        this->state[1] = y_;
        this->state[2] = z_;
        this->state[3] = vx_;
        this->state[4] = vy_;
        this->state[5] = vz_;
    }

    Point(std::vector<Type> &in) : state(6)
    {
        this->state = in;
    }

    Point(const Point& in)
    {
        this->state = in.state;
    }

    // Operator overloads
    Point& operator+=(Point const& rhs)
    {
        for (unsigned i = 0; i < this->state.size(); ++i) this->state[i] += rhs.state[i];

        return *this;
    }

    Point& operator-=(Point const& rhs)
    {
        for (unsigned i = 0; i < this->state.size(); ++i) this->state[i] -= rhs.state[i];

        return *this;
    }

    template <typename multType>
    Point& operator*=(multType a)
    {
        for (unsigned i = 0; i < this->state.size(); ++i) this->state[i] *= a;

        return *this;        
    }

    template <typename multType>
    Point& operator/=(multType a)
    {
        for (unsigned i = 0; i < this->state.size(); ++i) this->state[i] /= a;
    }

    Point operator+(Point const &a) const
    {
        std::vector<Type> ret(state.size());
        for (unsigned i = 0; i < this->state.size(); ++i) ret[i] = this->state[i] + a[i];

        return Point(ret);
    }

    Point operator-(Point const &a) const
    {
        std::vector<Type> ret(state.size());
        for (unsigned i = 0; i < this->state.size(); ++i) ret[i] = this->state[i] - a[i];

        return Point(ret);
    }

    bool operator==(const Point& a) const
    {
        return this->state == a.state;
    }

    void vectorToPoint(const std::vector<Type> &in)
    {
        this->state = in;
    }

    Type& operator[](int a)
    {
        return this->state[a];
    }

    const Type& operator[](int a) const
    {
        return this->state[a];
    }
};

template <typename Type>
std::ostream& operator<<(std::ostream& os, const Point<Type>& in)
{
    os << "(r, v): (";
    for (unsigned i = 0; i < 5; ++i) os << in[i] << ", ";
    return os << in[5] << ")";
};

#endif // _POINT_H_