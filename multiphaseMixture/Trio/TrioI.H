/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Swap.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<typename T1, typename T2, typename T3>
using trio = std::tuple<T1, T2, T3>;

template<class T>
inline int Foam::Trio<T>::compare(const Trio<T>& a, const Trio<T>& b)
{
    if (a.first() == b.first() && a.second() == b.second() && a.third() == b.third())
    {
        return 1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::Trio<T>::Trio(const T& f, const T& s, const T& t)
{
    first()  = f;
    second() = s;
	third()  = t;
}


template<class T>
inline Foam::Trio<T>::Trio(T&& f, T&& s, T&& t)
{
    first()  = std::move(f);
    second() = std::move(s);
    third()  = std::move(t);
}


template<class T>
inline Foam::Trio<T>::Trio(const trio<T,T,T>& vals)
{
    first()  = vals.first;
    second() = vals.second;
    third()  = vals.third;
}


template<class T>
inline Foam::Trio<T>::Trio(trio<T,T,T>&& vals)
{
    first()  = std::move(vals.first);
    second() = std::move(vals.second);
    third()  = std::move(vals.third);
}


template<class T>
inline Foam::Trio<T>::Trio(const FixedList<T, 3>& list)
:
    FixedList<T, 3>(list)
{}


/*template<class T>
inline Foam::Trio<T>::Trio(const T& f, const T& s, const T& t, const bool doSort)
{
	first()  = f;
	second() = s;
	third()  = t;
}


template<class T>
inline Foam::Trio<T>::Trio(const FixedList<T, 3>& list, const bool doSort)
:
    Trio<T>(list.first(), list.middle(), list.last(), doSort)
{}*/


template<class T>
inline Foam::Trio<T>::Trio(Istream& is)
:
    FixedList<T, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const T& Foam::Trio<T>::second() const noexcept
{
	const T* TT = cdata();//this;
    return TT[1];
}


template<class T>
inline T& Foam::Trio<T>::second() noexcept
{
	T* TT = data();//this;
    return TT[1];
}


template<class T>
inline const T& Foam::Trio<T>::third() const noexcept
{
    return last();
}


template<class T>
inline T& Foam::Trio<T>::third() noexcept
{
    return last();
}


template<class T>
inline const T& Foam::Trio<T>::previous(const T& a) const
{
    if (a == first())
    {
        return third();
    }
	else if (a == second())
    {
        return first();
    }
	else if (a == third())
    {
        return second();
    }
    else
    {
        FatalErrorInFunction
            << "Trio " << *this
            << " does not contain " << a << abort(FatalError);
    }

    return first();
}


template<class T>
inline const T& Foam::Trio<T>::next(const T& a) const
{
    if (a == first())
    {
        return second();
    }
	else if (a == second())
    {
        return third();
    }
	else if (a == third())
    {
        return first();
    }
    else
    {
        FatalErrorInFunction
            << "Trio " << *this
            << " does not contain " << a << abort(FatalError);
    }

    return first();
}


/*template<class T>
inline void Foam::Trio<T>::flip()
{
    Foam::Swap(first(), second());
}


template<class T>
inline bool Foam::Trio<T>::sorted() const
{
    return !(second() < first());
}


template<class T>
inline void Foam::Trio<T>::sort()
{
    if (second() < first())
    {
        flip();
    }
}*/


// ************************************************************************* //
