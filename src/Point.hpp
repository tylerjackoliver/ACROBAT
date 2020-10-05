#ifndef _POINT_H_
#define _POINT_H_

#include <ostream>
#include <stdexcept>

template <typename Type>
struct Point
{
    Type x;
    Type y;
    Type z;
    Type vx;
    Type vy;
    Type vz;

    // Class constructors

    Point(){};  // Default initialiser;leave x, y, z to be created later
    Point(Type x_, Type y_, Type z_)
    {
        this->x = x_;
        this->y = y_;
        this->z = z_;
        this->vx = 0.0;
        this->vy = 0.0;
        this->vz = 0.0;
    }
    
    Point(Type x_, Type y_, Type z_, Type vx_, Type vy_, Type vz_)
    {
        this->x = x_;
        this->y = y_;
        this->z = z_;
        this->vx = vx_;
        this->vy = vy_;
        this->vz = vz_;
    }

    Point(const Point& in)
    {
        this->x = in.x;
        this->y = in.y;
        this->z = in.z;
        this->vx = in.vx;
        this->vy = in.vy;
        this->vz = in.vz;
    }

    // Operator overloads
    Point& operator+=(Point const& rhs)
    {
        this->x += rhs.x;
        this->y += rhs.y;
        this->z += rhs.z;
        this->vx += rhs.vx;
        this->vy += rhs.vy;
        this->vz += rhs.vz;

        return *this;
    }

    Point& operator-=(Point const& rhs)
    {
        this->x -= rhs.x;
        this->y -= rhs.y;
        this->z -= rhs.z;
        this->vx -= rhs.vx;
        this->vy -= rhs.vy;
        this->vz -= rhs.vz;

        return *this;
    }

    template <typename multType>
    Point& operator*=(multType a)
    {
        this->x *= rhs.x;
        this->y *= rhs.y;
        this->z *= rhs.z;
        this->vx *= rhs.vx;
        this->vy *= rhs.vy;
        this->vz *= rhs.vz;

        return *this;        
    }

    template <typename multType>
    Point& operator/=(multType a)
    {
        this->x /= rhs.x;
        this->y /= rhs.y;
        this->z /= rhs.z;
        this->vx /= rhs.vx;
        this->vy /= rhs.vy;
        this->vz /= rhs.vz;
    }

    Point operator+(Point const &a) const
    {
        Type xx = this->x + a.x;
        Type yy = this->y + a.y;
        Type zz = this->z + a.z;
        Type vxx = this->vx + a.vx;
        Type vyy = this->vy + a.vy;
        Type vzz = this->vz + a.vz;

        return Point(xx, yy, zz, vxx, vyy, vzz);
    }

    Point operator-(Point const &a) const
    {
        Type xx = this->x - a.x;
        Type yy = this->y - a.y;
        Type zz = this->z - a.z;
        Type vxx = this->vx - a.vx;
        Type vyy = this->vy - a.vy;
        Type vzz = this->vz - a.vz;

        return Point(xx, yy, zz, vxx, vyy, vzz);
    }

    template <typename multType>
    Point operator*(multType a)
    {
        Type xx = this->x * a.x;
        Type yy = this->y * a.y;
        Type zz = this->z * a.z;
        Type vxx = this->vx * a.vx;
        Type vyy = this->vy * a.vy;
        Type vzz = this->vz * a.vz;

        return Point(xx, yy, zz, vxx, vyy, vzz);
    }

    bool operator==(const Point& a) const
    {
        bool xx, yy, zz, vxx, vyy, vzz, ret;
        xx = this->x == a.x;
        yy = this->y == a.y;
        zz = this->z == a.z;
        vxx = this->vx == a.vx;
        vyy = this->vy == a.vy;
        vzz = this->vz == a.vz;

        ret = xx && yy && zz && vxx && vyy && vzz;

        return ret;
    }

    Point& operator=(const Point& a)
    {
        this->x = a.x;
        this->y = a.y;
        this->z = a.z;
        this->vx = a.vx;
        this->vy = a.vy;
        this->vz = a.vz;

        return *this;
    }

    void vectorToPoint(const std::vector<Type> &in)
    {
        this->x = in[0];
        this->y = in[1];
        this->z = in[2];
        this->vx = in[3];
        this->vy = in[4];
        this->vz = in[5];
    }

    Type& operator[](int a)
    {
        Type ret;
        switch(a)
        {
            case 0:
                return this->x; // I'm so sorry, this is disgusting
                break;
            case 1:
                return this->y;
                break;
            case 2:
                return this->z;
                break;
            case 3:
                return this->vx;
                break;
            case 4:
                return this->vy;
                break;
            case 5:
                return this->vz;
                break;
            default:
                throw std::out_of_range("Point index out of range."); // Throw error in bad case
        }
    }

    const Type& operator[](int a) const
    {
        Type ret;
        switch(a)
        {
            case 0:
                return this->x;
                break;
            case 1:
                return this->y;
                break;
            case 2:
                return this->z;
                break;
            case 3:
                return this->vx;
                break;
            case 4:
                return this->vy;
                break;
            case 5:
                return this->vz;
                break;
            default:
                throw std::out_of_range("Point index out of range."); // Throw error in bad case
        }
    }

};

template <typename Type>
std::ostream& operator<<(std::ostream& os, const Point<Type>& in)
{
    return (os << "(r, v): (" << in.x << ", " << in.y << ", " << in.z << ", " << in.vx << ", " << in.vy << ", " << in.vz << ")");    
}

#endif // _POINT_H_
